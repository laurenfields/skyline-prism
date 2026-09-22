using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// The hit rule: which p, what cut, how big an effect. Read in four places, so the thing worth
/// pinning is that one object decides all of them the same way.
/// </summary>
public class SignificanceRuleTests
{
    private static DifferentialRow Row(double logFc, double p, double adjP) =>
        new("F1", logFc, 0, 0, 0, p, adjP, 0, 0);

    [Fact]
    public void TheDefault_IsTheConventionalRule()
    {
        var rule = SignificanceRule.Default;

        Assert.Equal(0.05, rule.PThreshold);
        Assert.True(rule.UseAdjusted);
        Assert.Equal(1.0, rule.Log2FcThreshold);
        Assert.Equal("adj.P < 0.05, |log2FC| >= 1", rule.Describe());
    }

    [Fact]
    public void BothConditionsMustHold()
    {
        var rule = SignificanceRule.Default;

        Assert.True(rule.IsSignificant(Row(logFc: 1.5, p: 1e-6, adjP: 0.01)));
        Assert.False(rule.IsSignificant(Row(logFc: 0.5, p: 1e-6, adjP: 0.01)));  // effect too small
        Assert.False(rule.IsSignificant(Row(logFc: 1.5, p: 1e-6, adjP: 0.20)));  // p too large
    }

    [Fact]
    public void RawVersusAdjusted_ChangesWhichColumnIsJudged()
    {
        // The case the option exists for: nothing survives correction, but the raw ranking is not
        // empty. Same feature, opposite verdicts.
        var row = Row(logFc: 1.5, p: 0.001, adjP: 0.40);

        Assert.False(new SignificanceRule().IsSignificant(row));
        Assert.True(new SignificanceRule { UseAdjusted = false }.IsSignificant(row));
    }

    [Fact]
    public void AZeroEffectCut_LeavesThePValueAlone_AndSaysSo()
    {
        var rule = new SignificanceRule { Log2FcThreshold = 0 };

        Assert.True(rule.IsSignificant(Row(logFc: 0.01, p: 1e-9, adjP: 1e-6)));
        // No "|log2FC| >= 0" clause, which would read as a constraint rather than its absence.
        Assert.Equal("adj.P < 0.05", rule.Describe());
    }

    [Fact]
    public void ANaNInEitherQuantity_IsNeverAHit()
    {
        var rule = SignificanceRule.Default;

        Assert.False(rule.IsSignificant(Row(logFc: double.NaN, p: 1e-9, adjP: 1e-9)));
        Assert.False(rule.IsSignificant(Row(logFc: 2.0, p: 1e-9, adjP: double.NaN)));
    }

    [Theory]
    // Adjusted p asked for AND a correction ran: the axis is the adjusted value.
    [InlineData(true, true, "-log10(adjusted p-value)")]
    // Correct = None, so the adjusted column just holds the raw p - calling it adjusted would lie.
    [InlineData(true, false, "-log10(p-value)")]
    // The reader asked for raw explicitly.
    [InlineData(false, true, "-log10(p-value)")]
    [InlineData(false, false, "-log10(p-value)")]
    public void TheAxisNamesWhicheverPIsActuallyPlotted(bool useAdjusted, bool corrected, string expected)
        => Assert.Equal(expected, new SignificanceRule { UseAdjusted = useAdjusted }.YAxisLabel(corrected));

    [Fact]
    public void Describe_TakesTheEffectName_SoATrendDoesNotSayFoldChange()
    {
        var rule = new SignificanceRule { PThreshold = 0.01, UseAdjusted = false, Log2FcThreshold = 0.585 };

        Assert.Equal("P < 0.01, |log2 change across week| >= 0.585",
            rule.Describe("log2 change across week"));
    }

    [Fact]
    public void PValueOf_ReportsTheColumnTheRuleJudgesBy()
    {
        var row = Row(logFc: 1, p: 0.002, adjP: 0.04);

        Assert.Equal(0.04, new SignificanceRule().PValueOf(row));
        Assert.Equal(0.002, new SignificanceRule { UseAdjusted = false }.PValueOf(row));
    }
}
