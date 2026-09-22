using System;
using System.Collections.Generic;
using System.Linq;
using SkylinePrism.Core.Numerics;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// Per-feature prior variances for the moderated t, in the estimator the lab's
/// <c>proteomics-toolkit</c> uses rather than limma's.
/// </summary>
/// <remarks>
/// <para><b>Why this exists beside <see cref="EmpiricalBayes.SqueezeVarTrend"/>.</b> Both are called
/// an "intensity trend" and they are different estimators. limma's <c>trend=TRUE</c> fits a natural
/// cubic spline of log(residual variance) against mean <b>log2</b> expression and re-estimates the
/// prior degrees of freedom from that fit. The toolkit fits a LOWESS of log(within-group variance)
/// against log(within-group mean) on <b>raw, pre-log</b> intensities, one point per (feature, group),
/// converts back to log space by the delta method, and leaves the prior degrees of freedom at the
/// global value. They differ in the smoother, the space, what contributes a point, and whether the
/// prior df moves - so they disagree by a median 3-7% on p-values and up to 178% on individual
/// features (<c>docs/differential-analysis.md</c>). Neither is wrong; they are not interchangeable,
/// and a result has to say which one produced it.</para>
/// <para><b>Why raw intensities are free here.</b> The toolkit has to carry a separate pre-log copy
/// of the matrix because its dispatcher log-transforms in place. PRISM's matrix is log2 of a LINEAR
/// parquet, so the raw intensity is just <c>2^x</c> - exact, since that is the transform that
/// produced it.</para>
/// </remarks>
internal static class VariancePriors
{
    /// <summary>The toolkit's minimum: fewer points than this and the trend is not a fit, it is noise.</summary>
    private const int MinTrendPoints = 5;

    /// <summary>One (feature, group) observation: the raw-scale mean, sd and sample count.</summary>
    private readonly record struct GroupStat(int Feature, double Mean, double Sd, int N);

    /// <summary>
    /// The toolkit's <c>intensity_trend</c> prior: per-feature prior variance in LOG2 space, or null
    /// when there are too few usable points to fit a trend (the caller then falls back to the global
    /// prior and says so).
    /// </summary>
    /// <param name="exprLog2Tested">
    /// The LOG2 matrix actually fitted, <c>[feature, sample]</c>, already reduced to the tested
    /// features and the selected sample columns - so its column order matches
    /// <paramref name="groups"/>.
    /// </param>
    /// <param name="groups">
    /// Column indices INTO <paramref name="exprLog2Tested"/>, one list per prior group. Normally the
    /// two contrast arms; with the QC/reference override, whichever replicates were nominated.
    /// </param>
    public static double[]? IntensityTrend(double[,] exprLog2Tested, IReadOnlyList<IReadOnlyList<int>> groups)
    {
        var nFeatures = exprLog2Tested.GetLength(0);
        var stats = CollectGroupStats(exprLog2Tested, groups);

        // The trend is fitted only on points that can carry one; every point is then PREDICTED from
        // it, including the ones excluded from the fit, which is what the reference does.
        var fitPoints = stats
            .Where(s => s.N >= 2 && double.IsFinite(s.Mean) && double.IsFinite(s.Sd) && s.Mean > 0 && s.Sd > 0)
            .ToList();
        if (fitPoints.Count < MinTrendPoints)
            return null;

        // LOWESS wants x ascending, and the interpolation below needs the sorted curve anyway.
        var ordered = fitPoints.OrderBy(s => Math.Log(s.Mean)).ToList();
        var x = ordered.Select(s => Math.Log(s.Mean)).ToArray();
        var y = ordered.Select(s => 2.0 * Math.Log(s.Sd)).ToArray(); // log(variance) = 2*log(sd)
        var yhat = Lowess.Fit(x, y, frac: 0.5, iterations: 3);

        // Sample-size-weighted mean of the per-group predicted log-space variance, per feature.
        var sum = new double[nFeatures];
        var weight = new double[nFeatures];
        foreach (var s in stats)
        {
            if (!(s.N > 0) || !double.IsFinite(s.Mean) || !(s.Mean > 0))
                continue;
            // Edge-clamped interpolation onto the fitted curve, matching numpy.interp with
            // left=ys[0], right=ys[-1]: outside the fitted range the trend is held flat rather than
            // extrapolated, which a LOWESS has no basis to do.
            var logVarRaw = Stats.Interp(Math.Log(s.Mean), x, yhat);
            // Delta method for x -> log2(x): var_log2 = var_raw / mean^2 / (ln 2)^2.
            var varLog2 = Math.Exp(logVarRaw) / (s.Mean * s.Mean) / (Ln2 * Ln2);
            if (!double.IsFinite(varLog2) || !(varLog2 > 0))
                continue;
            sum[s.Feature] += varLog2 * s.N;
            weight[s.Feature] += s.N;
        }

        var prior = new double[nFeatures];
        for (var f = 0; f < nFeatures; f++)
            prior[f] = weight[f] > 0 ? sum[f] / weight[f] : double.NaN;

        // A feature with no usable group takes the mean of those that had one. Leaving it NaN would
        // silently drop the feature from the result at the p-value step.
        var known = prior.Where(double.IsFinite).ToArray();
        if (known.Length == 0)
            return null;
        var fallback = known.Average();
        for (var f = 0; f < nFeatures; f++)
            if (!double.IsFinite(prior[f]) || !(prior[f] > 0))
                prior[f] = fallback;

        return prior;
    }

    private const double Ln2 = 0.6931471805599453;

    /// <summary>
    /// Per (feature, group) raw-scale mean, sd (ddof = 1) and count, from the log2 matrix. Groups with
    /// fewer than two samples contribute nothing - an sd needs two.
    /// </summary>
    private static List<GroupStat> CollectGroupStats(
        double[,] exprLog2Tested, IReadOnlyList<IReadOnlyList<int>> groups)
    {
        var nFeatures = exprLog2Tested.GetLength(0);
        var stats = new List<GroupStat>(nFeatures * Math.Max(groups.Count, 1));
        var buffer = new List<double>();

        foreach (var cols in groups)
        {
            if (cols.Count < 2)
                continue;
            for (var f = 0; f < nFeatures; f++)
            {
                buffer.Clear();
                foreach (var c in cols)
                {
                    var v = exprLog2Tested[f, c];
                    if (double.IsFinite(v))
                        buffer.Add(Math.Pow(2.0, v)); // back to the linear scale the parquet held
                }

                if (buffer.Count < 2)
                {
                    stats.Add(new GroupStat(f, double.NaN, double.NaN, buffer.Count));
                    continue;
                }

                var mean = 0.0;
                foreach (var v in buffer)
                    mean += v;
                mean /= buffer.Count;

                var ss = 0.0;
                foreach (var v in buffer)
                {
                    var d = v - mean;
                    ss += d * d;
                }

                stats.Add(new GroupStat(f, mean, Math.Sqrt(ss / (buffer.Count - 1)), buffer.Count));
            }
        }

        return stats;
    }
}
