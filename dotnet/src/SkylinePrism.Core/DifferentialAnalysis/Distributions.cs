using System;
using MathNet.Numerics.Distributions;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// The distribution tails the tests in this namespace refer their statistics to.
/// </summary>
/// <remarks>
/// <para>Each one existed already as a private line inside the single caller that needed it -
/// <c>Differential.ModeratedPValue</c> and <c>DetectionGlm.Chi2SurvivalDf1</c> - which was fine while
/// there was one test. With a menu of tests there is no longer a single caller, and a second private
/// copy is how two tests end up disagreeing in the last digits for no stated reason.</para>
/// <para>Every one is written as an UPPER tail evaluated on the negative side rather than as
/// <c>1 - CDF</c>. Far out in the tail the CDF is within an ulp of 1, so subtracting it from 1 throws
/// away every significant digit of the answer - which is exactly the region a p-value is read in. The
/// forms here are the ones scipy uses for the same reason.</para>
/// </remarks>
internal static class Distributions
{
    /// <summary>
    /// Two-sided p from a t statistic: <c>P(|T| &gt;= |t|)</c> on <paramref name="df"/> degrees of
    /// freedom, or the normal limit when <paramref name="df"/> is infinite (which is what an
    /// empirical-Bayes prior with no estimable excess variance produces).
    /// </summary>
    public static double TwoSidedT(double t, double df)
    {
        if (double.IsNaN(t) || double.IsNaN(df))
            return double.NaN;
        var absT = Math.Abs(t);
        return double.IsInfinity(df)
            ? 2.0 * Normal.CDF(0.0, 1.0, -absT)
            : 2.0 * new StudentT(0.0, 1.0, df).CumulativeDistribution(-absT);
    }

    /// <summary>
    /// chi-square survival function for df = 1: <c>P(X &gt; x) = erfc(sqrt(x/2))</c>, computed
    /// tail-accurately as <c>2 * Phi(-sqrt(x))</c> to match <c>scipy.stats.chi2.sf(x, 1)</c>.
    /// </summary>
    public static double Chi2SurvivalDf1(double x)
    {
        if (double.IsNaN(x))
            return double.NaN;
        if (x <= 0.0)
            return 1.0;
        return 2.0 * Normal.CDF(0.0, 1.0, -Math.Sqrt(x));
    }

    /// <summary>Standard normal two-sided p: <c>P(|Z| &gt;= |z|)</c>.</summary>
    public static double TwoSidedNormal(double z)
    {
        if (double.IsNaN(z))
            return double.NaN;
        return 2.0 * Normal.CDF(0.0, 1.0, -Math.Abs(z));
    }

    /// <summary>
    /// Welch-Satterthwaite degrees of freedom for two samples with unequal variances:
    /// <c>(v1/n1 + v2/n2)^2 / ((v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1))</c>.
    /// </summary>
    /// <remarks>
    /// Deliberately NOT rounded to an integer - it is a fractional quantity and scipy's
    /// <c>ttest_ind(equal_var=False)</c> refers the statistic to exactly this fractional df. Rounding
    /// it is a classic way to disagree with every other tool in the third decimal of a p-value.
    /// Returns NaN when either arm has no spread to estimate from, which the caller reports rather
    /// than papering over.
    /// </remarks>
    public static double WelchSatterthwaiteDf(double var1, int n1, double var2, int n2)
    {
        if (n1 < 2 || n2 < 2)
            return double.NaN;
        var a = var1 / n1;
        var b = var2 / n2;
        var denom = a * a / (n1 - 1) + b * b / (n2 - 1);
        if (!(denom > 0))
            return double.NaN; // both arms constant: the statistic is not defined either
        return (a + b) * (a + b) / denom;
    }
}
