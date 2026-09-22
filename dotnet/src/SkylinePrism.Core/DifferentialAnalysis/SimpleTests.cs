using System;
using System.Collections.Generic;
using System.Linq;
using SkylinePrism.Core.Numerics;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// The per-feature two-sample tests that are not empirical-Bayes: Welch's t, Student's t and
/// Mann-Whitney U.
/// </summary>
/// <remarks>
/// <para>Each one produces the same <see cref="DifferentialResult"/> shape as the moderated t, so the
/// Volcano, the hit table and the per-feature boxplot work unchanged - the only thing that differs is
/// how <c>t</c> and <c>PValue</c> were arrived at.</para>
/// <para><b>These do not borrow strength across features.</b> That is the entire point of choosing
/// one: an independent test per feature, with no shared variance prior, which is what a reader
/// expecting "a t-test" means. It is also why they are less powerful than the moderated t on the
/// small-n designs proteomics usually has, and why the moderated t stays the default.</para>
/// <para>Covariates are not accepted here. A t-test has no design matrix to put them in; adjusting
/// for a covariate means fitting a linear model, which is the moderated path. The caller reports
/// that rather than silently ignoring a ticked covariate.</para>
/// </remarks>
internal static class SimpleTests
{
    /// <summary>
    /// Run <paramref name="test"/> on every feature complete across the selected samples.
    /// </summary>
    public static DifferentialResult Run(
        double[,] exprLog2FeaturesBySamples,
        IReadOnlyList<string> featureIds,
        IReadOnlyList<int> groupAColumns,
        IReadOnlyList<int> groupBColumns,
        DifferentialTest test,
        MultipleTesting correction,
        int minPerGroup,
        IReadOnlyList<string> messages)
    {
        var nFeatures = exprLog2FeaturesBySamples.GetLength(0);
        if (featureIds.Count != nFeatures)
            throw new ArgumentException("featureIds length must match the matrix rows", nameof(featureIds));
        if (groupAColumns.Count < minPerGroup || groupBColumns.Count < minPerGroup)
            throw new ArgumentException($"each group needs at least {minPerGroup} samples");
        if (groupAColumns.Intersect(groupBColumns).Any())
            throw new ArgumentException("the two groups share a sample column");

        var nA = groupAColumns.Count;
        var nB = groupBColumns.Count;

        // Complete-case, exactly as the moderated path does: a feature is tested only where it is
        // present in every selected sample, so every feature is tested on the same samples.
        var tested = new List<int>(nFeatures);
        for (var f = 0; f < nFeatures; f++)
        {
            var ok = true;
            foreach (var c in groupAColumns.Concat(groupBColumns))
                if (!double.IsFinite(exprLog2FeaturesBySamples[f, c]))
                {
                    ok = false;
                    break;
                }

            if (ok)
                tested.Add(f);
        }

        if (tested.Count == 0)
            throw new InvalidOperationException("No feature is observed across every selected sample.");

        var a = new double[nA];
        var b = new double[nB];
        var rows = new DifferentialRow[tested.Count];
        var pValues = new double[tested.Count];

        for (var i = 0; i < tested.Count; i++)
        {
            var f = tested[i];
            for (var j = 0; j < nA; j++)
                a[j] = exprLog2FeaturesBySamples[f, groupAColumns[j]];
            for (var j = 0; j < nB; j++)
                b[j] = exprLog2FeaturesBySamples[f, groupBColumns[j]];

            var (stat, p, logFc) = test switch
            {
                DifferentialTest.WelchT => TTest(a, b, pooled: false),
                DifferentialTest.StudentT => TTest(a, b, pooled: true),
                DifferentialTest.MannWhitney => MannWhitneyU(a, b),
                _ => throw new ArgumentOutOfRangeException(nameof(test), test, "not a simple test"),
            };

            pValues[i] = p;
            var meanA = NumpyMath.Mean(a);
            var meanB = NumpyMath.Mean(b);
            rows[i] = new DifferentialRow(
                featureIds[f], logFc, Math.Pow(2.0, logFc), (meanA * nA + meanB * nB) / (nA + nB),
                stat, p, double.NaN, meanA, meanB);
        }

        var adjusted = Fdr.Adjust(pValues, correction);
        for (var i = 0; i < rows.Length; i++)
            rows[i] = rows[i] with { AdjPValue = adjusted[i] };

        return new DifferentialResult(
            rows.OrderBy(r => r.PValue).ToList(), nA, nB, nFeatures, tested.Count,
            // No linear model, so no residual df shared across features (Welch's differs per
            // feature) and no prior at all. NaN rather than a number that would be read as one.
            double.NaN, double.NaN, "none", Array.Empty<string>(), messages, Array.Empty<string>());
    }

    /// <summary>
    /// Two-sample t. <paramref name="pooled"/> false is Welch (unequal variances, Satterthwaite df),
    /// true is Student (pooled variance, <c>nA + nB - 2</c> df). The contrast is B - A, so a positive
    /// value means higher in B, matching the moderated path and the axis label.
    /// </summary>
    private static (double T, double P, double LogFc) TTest(double[] a, double[] b, bool pooled)
    {
        int nA = a.Length, nB = b.Length;
        var meanA = NumpyMath.Mean(a);
        var meanB = NumpyMath.Mean(b);
        var varA = NumpyMath.Var(a, ddof: 1);
        var varB = NumpyMath.Var(b, ddof: 1);
        var logFc = meanB - meanA;

        double se, df;
        if (pooled)
        {
            df = nA + nB - 2;
            var pooledVar = ((nA - 1) * varA + (nB - 1) * varB) / df;
            se = Math.Sqrt(pooledVar * (1.0 / nA + 1.0 / nB));
        }
        else
        {
            se = Math.Sqrt(varA / nA + varB / nB);
            df = Distributions.WelchSatterthwaiteDf(varA, nA, varB, nB);
        }

        // Both arms constant: the statistic is 0/0. NaN rather than an infinity that would sort to
        // the top of the hit table as the most significant thing in the run.
        if (!(se > 0) || !double.IsFinite(df))
            return (double.NaN, double.NaN, logFc);

        var t = logFc / se;
        return (t, Distributions.TwoSidedT(t, df), logFc);
    }

    /// <summary>
    /// Mann-Whitney U, two-sided, by the normal approximation with a tie correction and a continuity
    /// correction - <c>scipy.stats.mannwhitneyu(..., method="asymptotic")</c>.
    /// </summary>
    /// <remarks>
    /// <para><b>Always asymptotic, never the exact permutation distribution.</b> scipy's
    /// <c>method="auto"</c> switches to the exact test when the larger sample is 8 or fewer AND there
    /// are no ties; above that - which is every realistic cohort - it uses this. PRISM does not
    /// implement the exact branch, so on a very small group its p-value will differ from a scipy
    /// default-argument run. The caller says so rather than leaving the reader to discover it.</para>
    /// <para>The reported effect is the difference in MEDIANS, not means: a rank test makes no claim
    /// about means, and reporting one beside a rank p-value invites reading the pair as a t-test.</para>
    /// </remarks>
    private static (double U, double P, double LogFc) MannWhitneyU(double[] a, double[] b)
    {
        int nA = a.Length, nB = b.Length;
        var all = new double[nA + nB];
        Array.Copy(a, 0, all, 0, nA);
        Array.Copy(b, 0, all, nA, nB);

        var ranks = Stats.RankAverage(all);
        var rankSumB = 0.0;
        for (var i = 0; i < nB; i++)
            rankSumB += ranks[nA + i];

        // U for B against A, so the sign convention matches logFc = median(B) - median(A).
        var uB = rankSumB - nB * (nB + 1) / 2.0;
        var mean = nA * (double)nB / 2.0;

        // Tie correction: sum over tied groups of (t^3 - t), which shrinks the variance.
        var n = nA + nB;
        var tieTerm = 0.0;
        foreach (var run in all.OrderBy(v => v).GroupBy(v => v))
        {
            double t = run.Count();
            if (t > 1)
                tieTerm += t * t * t - t;
        }

        var variance = nA * (double)nB / 12.0 * (n + 1 - tieTerm / (n * (double)(n - 1)));
        var logFc = Median(b) - Median(a);
        if (!(variance > 0))
            return (uB, double.NaN, logFc); // every value tied: no rank information at all

        // Continuity correction of 0.5 toward the mean, as scipy's use_continuity default does.
        var z = (Math.Abs(uB - mean) - 0.5) / Math.Sqrt(variance);
        if (z < 0)
            z = 0;
        return (uB, Distributions.TwoSidedNormal(z), logFc);
    }

    private static double Median(double[] values)
    {
        var copy = (double[])values.Clone();
        Array.Sort(copy);
        var n = copy.Length;
        return n % 2 == 1 ? copy[n / 2] : 0.5 * (copy[n / 2 - 1] + copy[n / 2]);
    }
}
