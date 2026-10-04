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
/// A trend design has no A and B arms, so the two-arm entry point must refuse it and say where to
/// go - not fall through and run an ordinary contrast under a trend's name.
/// </summary>
public class TrendDesignRoutingTests
{
    private static readonly double[,] Expr = { { 1, 2, 3, 4 }, { 2, 3, 4, 5 } };
    private static readonly string[] Ids = { "F1", "F2" };

    [Theory]
    [InlineData(DifferentialDesign.LinearTrend)]
    [InlineData(DifferentialDesign.LinearTrendWithinSubject)]
    public void TheTwoArmEntryPoint_RefusesATrend_AndNamesTheRightOne(DifferentialDesign design)
    {
        var ex = Assert.Throws<ArgumentException>(() => Differential.Run(
            Expr, Ids, new[] { 0, 1 }, new[] { 2, 3 }, new DifferentialOptions { Design = design }));

        Assert.Contains("no A and B arms", ex.Message);
        Assert.Contains("RunTrend", ex.Message);
    }

    [Theory]
    [InlineData(DifferentialDesign.Unpaired)]
    [InlineData(DifferentialDesign.Paired)]
    public void TheTwoArmDesigns_AreNotCaughtByThatCheck(DifferentialDesign design)
    {
        var options = new DifferentialOptions
        {
            Design = design,
            SubjectLabels = design == DifferentialDesign.Paired
                ? new string?[] { "s1", "s2", "s1", "s2" }
                : null,
        };

        var ex = Record.Exception(() => Differential.Run(
            Expr, Ids, new[] { 0, 1 }, new[] { 2, 3 }, options));

        Assert.False(ex is ArgumentException { Message: var m } && m.Contains("RunTrend"));
    }

    [Fact]
    public void RunTrend_RefusesATwoArmDesign()
    {
        var ex = Assert.Throws<ArgumentException>(() => Differential.RunTrend(
            Expr, Ids, new[] { 0, 1, 2, 3 }, new[] { 0.0, 1.0, 2.0, 3.0 },
            new DifferentialOptions { Design = DifferentialDesign.Unpaired }));

        Assert.Contains("linear-trend design", ex.Message);
    }

    [Fact]
    public void ATrendColumnWithOneValue_IsRefused_RatherThanFittingARankDeficientDesign()
    {
        var ex = Assert.Throws<ArgumentException>(() => Differential.RunTrend(
            Expr, Ids, new[] { 0, 1, 2, 3 }, new[] { 7.0, 7.0, 7.0, 7.0 },
            new DifferentialOptions { Design = DifferentialDesign.LinearTrend }));

        Assert.Contains("only one value", ex.Message);
    }

    [Fact]
    public void AWithinSubjectTrend_WithoutASubjectColumn_SaysSo()
    {
        var ex = Assert.Throws<ArgumentException>(() => Differential.RunTrend(
            Expr, Ids, new[] { 0, 1, 2, 3 }, new[] { 0.0, 1.0, 2.0, 3.0 },
            new DifferentialOptions { Design = DifferentialDesign.LinearTrendWithinSubject }));

        Assert.Contains("subject column", ex.Message);
    }
}

/// <summary>
/// Where the intensity-trend prior takes its groups from, and why it matters.
/// </summary>
/// <remarks>
/// The prior's per-feature scale should describe MEASUREMENT variance. Fitted on the design groups
/// of a real study it describes measurement variance plus the biology those groups contain, and the
/// moderation then shrinks the very effects the analysis is looking for. These pin that the two
/// sources genuinely produce different numbers, and that the result records which one ran - two
/// results are not comparable unless they used the same source.
/// </remarks>
public class VariancePriorSourceTests
{
    // Two arms whose WITHIN-arm spread is large (biological heterogeneity), plus a set of control
    // replicates whose spread is small (technical only). The two prior sources must disagree.
    private static (double[,] Expr, string[] Ids, int[] A, int[] B, IReadOnlyList<IReadOnlyList<int>> Controls)
        Cohort()
    {
        const int nFeatures = 40, nPerArm = 8, nControls = 6;
        var rng = new Random(4242);
        var total = 2 * nPerArm + nControls;
        var expr = new double[nFeatures, total];
        for (var f = 0; f < nFeatures; f++)
        {
            var level = 10 + f * 0.25;
            for (var s = 0; s < nPerArm; s++)
            {
                // Wide within-arm spread: this is the biology the prior must not absorb.
                expr[f, s] = level + 1.2 * (rng.NextDouble() - 0.5) * 2;
                expr[f, nPerArm + s] = level + 0.8 + 1.2 * (rng.NextDouble() - 0.5) * 2;
            }

            for (var s = 0; s < nControls; s++)
                // Tight replicate spread: this is the measurement variance.
                expr[f, 2 * nPerArm + s] = level + 0.12 * (rng.NextDouble() - 0.5) * 2;
        }

        var controls = new IReadOnlyList<int>[]
        {
            Enumerable.Range(2 * nPerArm, nControls).ToArray(),
        };
        return (expr, Enumerable.Range(0, nFeatures).Select(i => $"f{i}").ToArray(),
            Enumerable.Range(0, nPerArm).ToArray(), Enumerable.Range(nPerArm, nPerArm).ToArray(),
            controls);
    }

    private static DifferentialResult Run(IReadOnlyList<IReadOnlyList<int>>? priorGroups)
    {
        var (expr, ids, a, b, _) = Cohort();
        var (_, _, _, _, controls) = Cohort();
        return Differential.Run(expr, ids, a, b, new DifferentialOptions
        {
            Prior = VariancePrior.IntensityTrend,
            Correction = MultipleTesting.None,
            PriorGroupColumns = priorGroups is null ? null : controls,
        });
    }

    [Fact]
    public void TheResultRecordsWhichSourceTheScaleCameFrom()
    {
        var (_, _, _, _, controls) = Cohort();

        Assert.Equal("intensity-trend from design groups", Run(null).VariancePrior);
        Assert.Equal("intensity-trend from controls", Run(controls).VariancePrior);
    }

    [Fact]
    public void ControlsAndDesignGroups_GiveDifferentNumbers()
    {
        var (_, _, _, _, controls) = Cohort();
        var fromGroups = Run(null);
        var fromControls = Run(controls);

        // The prior degrees of freedom are estimated from the study residuals either way, so the
        // AMOUNT of shrinkage is unchanged - only the target moves.
        Assert.Equal(fromGroups.DfResidual, fromControls.DfResidual, 12);

        var g = fromGroups.Rows.ToDictionary(r => r.FeatureId);
        var c = fromControls.Rows.ToDictionary(r => r.FeatureId);

        // The effect sizes come from the same OLS fit and must not move at all...
        foreach (var id in g.Keys)
            Assert.Equal(g[id].LogFc, c[id].LogFc, 12);

        // ...while the moderated statistics must, or the option is doing nothing.
        Assert.Contains(g.Keys, id => Math.Abs(g[id].T - c[id].T) > 1e-6);
    }

    [Fact]
    public void ATighterControlPrior_ShrinksTowardSmallerVariance_SoTIsLarger()
    {
        // The direction the lab's practice exists to produce: a prior taken from nominal replicates
        // is smaller than one taken from groups carrying biology, so genuine effects survive
        // moderation instead of being shrunk toward nothing.
        var (_, _, _, _, controls) = Cohort();
        var g = Run(null).Rows.ToDictionary(r => r.FeatureId);
        var c = Run(controls).Rows.ToDictionary(r => r.FeatureId);

        var larger = g.Keys.Count(id => Math.Abs(c[id].T) > Math.Abs(g[id].T));
        Assert.True(larger > g.Count / 2,
            $"expected the control prior to raise |t| for most features; it raised {larger}/{g.Count}");
    }
}
