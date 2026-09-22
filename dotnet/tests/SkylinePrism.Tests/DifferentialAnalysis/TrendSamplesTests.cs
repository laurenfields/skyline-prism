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
