using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Which samples a trend design can use. The exclusions matter as much as the inclusions: a sample
/// or subject that cannot carry slope information still costs a parameter if it is kept.
/// </summary>
public class TrendSamplesTests
{
    private static readonly int[] AllColumns = { 0, 1, 2, 3, 4, 5 };

    [Fact]
    public void Independent_KeepsEverySampleWithAValue()
    {
        var x = new[] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };

        var sel = TrendSamples.Resolve(AllColumns, x, null, withinSubject: false);

        Assert.Equal(AllColumns, sel.Columns);
        Assert.Equal(x, sel.X);
        Assert.Equal(0, sel.SubjectCount);
        Assert.Empty(sel.Messages);
    }

    [Fact]
    public void ASampleWithNoTrendValue_IsLeftOut_AndCounted()
    {
        // NaN is "this sample has no value in that column", which is not evidence about the slope.
        var x = new[] { 0.0, double.NaN, 2.0, 3.0, double.NaN, 5.0 };

        var sel = TrendSamples.Resolve(AllColumns, x, null, withinSubject: false);

        Assert.Equal(new[] { 0, 2, 3, 5 }, sel.Columns);
        Assert.Contains(sel.Messages, m => m.Contains("2 sample(s) have no value"));
    }

    [Fact]
    public void WithinSubject_GroupsByLabel_NotBySampleOrder()
    {
        // Interleaved, the way a real metadata column arrives.
        var x = new[] { 0.0, 0.0, 1.0, 1.0, 2.0, 2.0 };
        var subjects = new string?[] { "A", "B", "A", "B", "A", "B" };

        var sel = TrendSamples.Resolve(AllColumns, x, subjects, withinSubject: true);

        Assert.Equal(2, sel.SubjectCount);
        Assert.Equal(6, sel.Columns.Count);
        // Each subject's three samples share one index, whatever order they arrived in.
        var bySubject = sel.Columns
            .Select((c, i) => (c, s: sel.SubjectOf[i]))
            .GroupBy(t => t.s)
            .ToDictionary(g => g.Key, g => g.Select(t => t.c).OrderBy(c => c).ToArray());
        Assert.Equal(new[] { 0, 2, 4 }, bySubject[bySubject.Keys.First(k => bySubject[k][0] == 0)]);
        Assert.Equal(new[] { 1, 3, 5 }, bySubject[bySubject.Keys.First(k => bySubject[k][0] == 1)]);
    }

    [Fact]
    public void ASubjectWithOneSample_CarriesNoSlope_AndIsLeftOut()
    {
        // Its dummy would absorb its single point entirely: no slope information, one parameter
        // spent. Excluded out loud.
        var x = new[] { 0.0, 1.0, 2.0, 0.0, 1.0, 2.0 };
        var subjects = new string?[] { "A", "A", "A", "B", "C", "D" };

        var sel = TrendSamples.Resolve(AllColumns, x, subjects, withinSubject: true);

        Assert.Equal(1, sel.SubjectCount);
        Assert.Equal(new[] { 0, 1, 2 }, sel.Columns);
        Assert.Contains(sel.Messages, m => m.Contains("3 subject(s)") && m.Contains("two distinct"));
    }

    [Fact]
    public void ASubjectMeasuredTwiceAtTheSameTime_AlsoCarriesNoSlope()
    {
        // Two samples is not enough - they have to be at two different values of x.
        var x = new[] { 0.0, 0.0, 0.0, 0.0, 1.0, 2.0 };
        var subjects = new string?[] { "A", "A", "A", "B", "B", "B" };

        var sel = TrendSamples.Resolve(AllColumns, x, subjects, withinSubject: true);

        Assert.Equal(1, sel.SubjectCount);
        Assert.Equal(new[] { 3, 4, 5 }, sel.Columns);
        Assert.Contains(sel.Messages, m => m.Contains("A"));
    }

    [Fact]
    public void ASampleWithNoSubject_IsLeftOut_AndCounted()
    {
        var x = new[] { 0.0, 1.0, 2.0, 0.0, 1.0, 2.0 };
        var subjects = new string?[] { "A", "A", "A", null, "", "B" };

        var sel = TrendSamples.Resolve(AllColumns, x, subjects, withinSubject: true);

        Assert.Equal(new[] { 0, 1, 2 }, sel.Columns);
        Assert.Contains(sel.Messages, m => m.Contains("no subject"));
    }

    [Fact]
    public void WithinSubject_WithoutLabels_Throws()
        => Assert.Throws<ArgumentException>(() => TrendSamples.Resolve(
            AllColumns, new[] { 0.0, 1, 2, 3, 4, 5 }, null, withinSubject: true));

    [Fact]
    public void SubjectIndices_AreContiguousFromZero_SoTheDummyBlockIsWellFormed()
    {
        // Three subjects survive out of five; their indices must be 0,1,2 - not the original
        // positions, which would leave gaps and an over-wide design.
        var x = new[] { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
        var subjects = new string?[] { "gap1", "gap1", "gap2", "gap2", "gap3", "gap3" };

        var sel = TrendSamples.Resolve(AllColumns, x, subjects, withinSubject: true);

        Assert.Equal(3, sel.SubjectCount);
        Assert.Equal(new[] { 0, 1, 2 }, sel.SubjectOf.Distinct().OrderBy(i => i));
    }
}

/// <summary>
/// The trend design is a strict GENERALIZATION of the two-arm contrast, and this is what makes one
/// effect-size threshold serve both.
/// </summary>
/// <remarks>
/// A trend reports the modeled change across the observed range of x. Code a two-level grouping as
/// x = 0 and x = 1 and that range is 1, so the reported change IS the log2 fold change and the
/// moderated t is the same test. If this ever stops holding, "|log2FC| >= 1" and
/// "|log2 change| >= 1" have quietly become two different thresholds wearing one number.
/// </remarks>
public class TrendGeneralizesTheTwoArmContrastTests
{
    private static (double[,] Expr, string[] Ids) Data(int nFeatures, int nPerArm, int seed)
    {
        var rng = new Random(seed);
        var expr = new double[nFeatures, 2 * nPerArm];
        for (var f = 0; f < nFeatures; f++)
        {
            var level = 12 + f * 0.37;
            var shift = f % 3 == 0 ? 1.4 : 0.0;
            for (var s = 0; s < nPerArm; s++)
            {
                expr[f, s] = level + 0.4 * (rng.NextDouble() - 0.5);
                expr[f, nPerArm + s] = level + shift + 0.4 * (rng.NextDouble() - 0.5);
            }
        }

        return (expr, Enumerable.Range(0, nFeatures).Select(i => $"f{i}").ToArray());
    }

    [Fact]
    public void ATwoLevelTrendCodedZeroOne_ReproducesTheTwoArmContrast()
    {
        const int nFeatures = 12, nPerArm = 9;
        var (expr, ids) = Data(nFeatures, nPerArm, seed: 20260922);
        var armA = Enumerable.Range(0, nPerArm).ToArray();
        var armB = Enumerable.Range(nPerArm, nPerArm).ToArray();

        // x = 0 for arm A, 1 for arm B: the same grouping, expressed as a number.
        var x = new double[2 * nPerArm];
        for (var s = nPerArm; s < 2 * nPerArm; s++)
            x[s] = 1.0;

        // The global prior on both sides - a trend has no groups for the intensity trend, so
        // asking for it would compare two different estimators rather than two designs.
        var options = new DifferentialOptions
        {
            Prior = VariancePrior.Global,
            Correction = MultipleTesting.BenjaminiHochberg,
        };

        var arms = Differential.Run(expr, ids, armA, armB, options);
        var trend = Differential.RunTrend(
            expr, ids, Enumerable.Range(0, 2 * nPerArm).ToArray(), x,
            options with { Design = DifferentialDesign.LinearTrend });

        Assert.Equal(1.0, trend.TrendRange, 12);
        Assert.Equal(arms.DfResidual, trend.DfResidual, 12);
        Assert.Equal(arms.DfPrior, trend.DfPrior, 9);

        var byArms = arms.Rows.ToDictionary(r => r.FeatureId);
        var byTrend = trend.Rows.ToDictionary(r => r.FeatureId);
        foreach (var id in ids)
        {
            Assert.Equal(byArms[id].LogFc, byTrend[id].LogFc, 10);
            Assert.Equal(byArms[id].T, byTrend[id].T, 10);
            Assert.Equal(byArms[id].PValue, byTrend[id].PValue, 10);
            Assert.Equal(byArms[id].AdjPValue, byTrend[id].AdjPValue, 10);
        }
    }

    [Fact]
    public void RescalingTheTrendColumn_LeavesTheReportedEffectUnchanged()
    {
        // The whole reason the reported effect is a change across the range rather than a slope:
        // measuring the same experiment in days instead of weeks must not move the hit list.
        const int nFeatures = 10, nPerArm = 8;
        var (expr, ids) = Data(nFeatures, nPerArm, seed: 7);
        var columns = Enumerable.Range(0, 2 * nPerArm).ToArray();
        var options = new DifferentialOptions
        {
            Design = DifferentialDesign.LinearTrend,
            Prior = VariancePrior.Global,
            Correction = MultipleTesting.None,
        };

        var weeks = Enumerable.Range(0, 2 * nPerArm).Select(i => (double)(i % 4)).ToArray();
        var days = weeks.Select(w => w * 7.0).ToArray();

        var inWeeks = Differential.RunTrend(expr, ids, columns, weeks, options);
        var inDays = Differential.RunTrend(expr, ids, columns, days, options);

        // The SPAN differs by exactly the unit change...
        Assert.Equal(inWeeks.TrendRange * 7.0, inDays.TrendRange, 10);

        // ...and the reported effect, and every p-value, do not move at all.
        var a = inWeeks.Rows.ToDictionary(r => r.FeatureId);
        var b = inDays.Rows.ToDictionary(r => r.FeatureId);
        foreach (var id in ids)
        {
            Assert.Equal(a[id].LogFc, b[id].LogFc, 10);
            Assert.Equal(a[id].PValue, b[id].PValue, 10);
        }
    }
}
