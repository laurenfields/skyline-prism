using System;
using MathNet.Numerics;

namespace SkylinePrism.Core.DifferentialAnalysis.Detection;

/// <summary>
/// McNemar's test of marginal homogeneity for paired binary outcomes.
/// </summary>
/// <remarks>
/// <para>The paired analogue of Fisher's exact test on detection. Only the DISCORDANT pairs carry
/// information: a peptide detected in both halves of a pair, or in neither, says nothing about
/// whether detection differs between the two conditions, because the subject is its own control.
/// What is tested is whether the two kinds of discordance are equally likely.</para>
/// <para>Exact, not the chi-square approximation - matching
/// <c>statsmodels.stats.contingency_tables.mcnemar</c>, whose default is <c>exact=True</c>. With
/// b + c discordant pairs and no difference, the number falling one way is Binomial(b + c, 1/2), so
/// the two-sided p is twice the lower tail. This is a genuine exact test at any size, so unlike the
/// rank tests in <see cref="SimpleTests"/> there is no small-sample caveat to carry.</para>
/// </remarks>
internal static class McNemar
{
    /// <summary>
    /// Two-sided exact p for <paramref name="b"/> pairs discordant one way and <paramref name="c"/>
    /// the other.
    /// </summary>
    /// <remarks>
    /// Returns 1 when there are no discordant pairs at all: every subject agreed with itself, which
    /// is the least possible evidence of a difference, not the most. Computed through log-gamma
    /// rather than factorials so a cohort with hundreds of discordant pairs does not overflow.
    /// </remarks>
    public static double TwoSidedP(int b, int c)
    {
        if (b < 0 || c < 0)
            throw new ArgumentOutOfRangeException(nameof(b), "discordant counts cannot be negative");

        var n = b + c;
        if (n == 0)
            return 1.0;

        var k = Math.Min(b, c);

        // P(X <= k) for X ~ Binomial(n, 1/2), summed in log space.
        var logHalfPow = -n * Math.Log(2.0);
        var lower = 0.0;
        for (var i = 0; i <= k; i++)
            lower += Math.Exp(LogChoose(n, i) + logHalfPow);

        var p = 2.0 * lower;
        return p > 1.0 ? 1.0 : p;
    }

    private static double LogChoose(int n, int k) =>
        SpecialFunctions.GammaLn(n + 1) - SpecialFunctions.GammaLn(k + 1) - SpecialFunctions.GammaLn(n - k + 1);
}
