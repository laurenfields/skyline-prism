using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// The arm resolution both the Differential pane and <c>prism differential</c> run through. It
/// decides which samples land on which side of a test, so the rules it refuses on matter as much as
/// the ones it accepts.
/// </summary>
public class ContrastArmsTests
{
    private static readonly string?[] Types =
        { "Control", "Mild", "Severe", "Control", "Severe", null, "Mild" };

    [Fact]
    public void SingleLevelPerArm_SelectsThoseColumns()
    {
        var arms = ContrastArms.Resolve(Types, new[] { "Control" }, new[] { "Severe" });

        Assert.True(arms.Ok);
        Assert.Null(arms.Error);
        Assert.Equal(new[] { 0, 3 }, arms.A);
        Assert.Equal(new[] { 2, 4 }, arms.B);
    }

    [Fact]
    public void SeveralLevels_AreTheUnion_InColumnOrder()
    {
        // The whole point of multi-select: Control + Mild is ONE arm, not two contrasts.
        var arms = ContrastArms.Resolve(Types, new[] { "Control", "Mild" }, new[] { "Severe" });

        Assert.True(arms.Ok);
        Assert.Equal(new[] { 0, 1, 3, 6 }, arms.A);
        Assert.Equal(new[] { 2, 4 }, arms.B);
    }

    [Fact]
    public void NullValuedSamples_JoinNeitherArm()
    {
        // A sample with no value for the group-by column is not evidence for either side. Column 5
        // is null above and must appear in neither list.
        var arms = ContrastArms.Resolve(Types, new[] { "Control", "Mild", "Severe" }, new[] { "Absent" });

        Assert.False(arms.Ok);
        Assert.DoesNotContain(5, arms.A);
    }

    [Fact]
    public void ALevelInBothArms_IsRefused_NotSilentlyDroppedFromOneSide()
    {
        // Dropping it from A and dropping it from B give different answers, and the output would
        // record neither choice - so neither is taken.
        var arms = ContrastArms.Resolve(Types, new[] { "Control", "Mild" }, new[] { "Mild", "Severe" });

        Assert.False(arms.Ok);
        Assert.Contains("Mild", arms.Error);
        Assert.Contains("both arms", arms.Error);
        Assert.Empty(arms.A);
        Assert.Empty(arms.B);
    }

    [Fact]
    public void EveryOverlappingLevel_IsNamed_NotJustTheFirst()
    {
        var arms = ContrastArms.Resolve(
            Types, new[] { "Control", "Mild", "Severe" }, new[] { "Mild", "Severe" });

        Assert.False(arms.Ok);
        Assert.Contains("Mild", arms.Error);
        Assert.Contains("Severe", arms.Error);
    }

    [Theory]
    [InlineData("Absent", "Severe", "Arm A matched no samples.")]
    [InlineData("Control", "Absent", "Arm B matched no samples.")]
    [InlineData("Absent", "AlsoAbsent", "Neither arm matched any sample.")]
    public void AnEmptyArm_SaysWhichOne(string a, string b, string expected)
    {
        // Which arm is empty is the difference between fixing one typo and re-reading the command.
        var arms = ContrastArms.Resolve(Types, new[] { a }, new[] { b });

        Assert.False(arms.Ok);
        Assert.Equal(expected, arms.Error);
    }

    [Theory]
    [InlineData(true, false)]
    [InlineData(false, true)]
    [InlineData(true, true)]
    public void NoLevelsPicked_AsksForSome(bool emptyA, bool emptyB)
    {
        var arms = ContrastArms.Resolve(
            Types,
            emptyA ? Array.Empty<string>() : new[] { "Control" },
            emptyB ? Array.Empty<string>() : new[] { "Severe" });

        Assert.False(arms.Ok);
        Assert.Contains("at least one value", arms.Error);
    }

    [Fact]
    public void Matching_IsCaseSensitive()
    {
        // Metadata levels come from one document's own annotation, where "control" and "Control" are
        // two distinct values. Folding them here would merge two groups a user kept apart.
        var arms = ContrastArms.Resolve(Types, new[] { "control" }, new[] { "Severe" });

        Assert.False(arms.Ok);
        Assert.Equal("Arm A matched no samples.", arms.Error);
    }

    [Fact]
    public void Describe_JoinsLevelsTheWayTheStatusLineReadsThem()
        => Assert.Equal("Control + Mild", ContrastArms.Describe(new[] { "Control", "Mild" }));
}

/// <summary>
/// The control-sample vocabulary, shared by the QC pane's group filter and the variance prior's
/// "fit on controls" option.
/// </summary>
public class ControlSampleTypesTests
{
    [Theory]
    [InlineData("qc", true)]
    [InlineData("QC", true)]
    [InlineData("reference", true)]
    [InlineData("Standard", true)]
    [InlineData("Quality Control", true)]
    [InlineData("STANDARD", true)]      // spellings vary in case between sources
    [InlineData("experimental", false)]
    [InlineData("Unknown", false)]
    [InlineData(null, false)]
    public void IsControl_CoversBothVocabularies(string? value, bool expected)
        => Assert.Equal(expected, ControlSampleTypes.IsControl(value));

    [Fact]
    public void PriorGroups_KeepsControlTypesApart()
    {
        // One group per control TYPE: pooling QC and reference would count the systematic gap
        // between two different materials as measurement noise.
        var types = new string?[] { "experimental", "qc", "reference", "qc", "reference", "experimental" };

        var groups = ControlSampleTypes.PriorGroups(types);

        Assert.NotNull(groups);
        Assert.Equal(2, groups!.Count);
        Assert.Contains(groups, g => g.SequenceEqual(new[] { 1, 3 }));
        Assert.Contains(groups, g => g.SequenceEqual(new[] { 2, 4 }));
    }

    [Fact]
    public void PriorGroups_DropsAControlTypeWithOnlyOneReplicate()
    {
        // A group of one has no variance to contribute at all.
        var types = new string?[] { "experimental", "qc", "reference", "qc" };

        var groups = ControlSampleTypes.PriorGroups(types);

        Assert.NotNull(groups);
        Assert.Single(groups!);
        Assert.Equal(new[] { 1, 3 }, groups![0]);
    }

    [Fact]
    public void PriorGroups_IsNullWhenNothingCanBeFit()
    {
        // Null, not empty: the caller turns the option off rather than fitting a prior on nothing.
        Assert.Null(ControlSampleTypes.PriorGroups(
            new string?[] { "experimental", "experimental", "qc" }));
    }
}

/// <summary>
/// The linear-trend design is declared but not implemented. It must refuse, not quietly run the
/// two-arm contrast under its name.
/// </summary>
public class UnimplementedDesignTests
{
    [Fact]
    public void LinearTrend_Throws_RatherThanRunningAnUnpairedContrast()
    {
        var expr = new double[,] { { 1, 2, 3, 4 }, { 2, 3, 4, 5 } };
        var options = new DifferentialOptions { Design = DifferentialDesign.LinearTrend };

        var ex = Assert.Throws<NotImplementedException>(() => Differential.Run(
            expr, new[] { "F1", "F2" }, new[] { 0, 1 }, new[] { 2, 3 }, options));

        // The message has to separate it from the variance prior of nearly the same name, which is
        // implemented and is the default - that collision is why the guard reads the way it does.
        Assert.Contains("linear-trend design", ex.Message);
        Assert.Contains("IntensityTrend", ex.Message);
    }

    [Theory]
    [InlineData(DifferentialDesign.Unpaired)]
    [InlineData(DifferentialDesign.Paired)]
    public void TheImplementedDesigns_AreNotCaughtByTheGuard(DifferentialDesign design)
    {
        var expr = new double[,] { { 1, 2, 3, 4 }, { 2, 3, 4, 5 } };
        var options = new DifferentialOptions
        {
            Design = design,
            SubjectLabels = design == DifferentialDesign.Paired
                ? new string?[] { "s1", "s2", "s1", "s2" }
                : null,
        };

        // Whatever else they do, they do not throw NotImplementedException.
        var ex = Record.Exception(() => Differential.Run(
            expr, new[] { "F1", "F2" }, new[] { 0, 1 }, new[] { 2, 3 }, options));

        Assert.IsNotType<NotImplementedException>(ex);
    }
}
