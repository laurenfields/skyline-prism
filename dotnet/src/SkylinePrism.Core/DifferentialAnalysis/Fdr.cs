using System;
using System.Collections.Generic;
using SkylinePrism.Core.Numerics;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// Benjamini-Hochberg false-discovery-rate correction, ported from the PRISM Differential
/// Explorer's <c>_bh</c> (prism_diff_explorer.py). NaN-safe: NaN p-values are excluded from the
/// correction (m counts only finite p-values) and passed through as NaN, bit-identical to that
/// reference. This is the standalone BH used for the "family" significance basis.
///
/// The explorer's genome-wide <c>adj.P.Val</c> instead comes from statsmodels
/// <c>multipletests(method="fdr_bh")</c>, which agrees with this to within rounding on finite
/// p-values (it scales by <c>rank/m</c> rather than <c>m/rank</c>) but, unlike this, returns all
/// NaN if any input is NaN and counts NaN toward m. The moderated-t path feeds only finite
/// p-values (complete rows, positive posterior variance), where the two are equivalent.
/// </summary>
public static class Fdr
{
    /// <summary>
    /// Benjamini-Hochberg adjusted p-values, same order as the input. NaN entries stay NaN and do
    /// not count toward the number of tests. Adjusted values are the monotone (reverse cumulative
    /// minimum) step-up transform, clipped to [0, 1].
    /// </summary>
    public static double[] BenjaminiHochberg(ReadOnlySpan<double> pValues)
    {
        var n = pValues.Length;
        var adj = new double[n];

        // Indices of the finite (testable) p-values; NaN positions are set aside as NaN.
        var okIdx = new List<int>(n);
        for (var i = 0; i < n; i++)
        {
            if (double.IsNaN(pValues[i]))
                adj[i] = double.NaN;
            else
                okIdx.Add(i);
        }

        var m = okIdx.Count;
        if (m == 0)
            return adj;

        var p = new double[m];
        for (var k = 0; k < m; k++)
            p[k] = pValues[okIdx[k]];

        // order[r] = index into p of the r-th smallest p-value (ascending).
        var order = Stats.ArgSort(p);

        // ranked[r] = p_(r) * m / (r+1)  -- the raw step-up values in ascending-p order.
        var ranked = new double[m];
        for (var r = 0; r < m; r++)
            ranked[r] = p[order[r]] * m / (r + 1);

        // Enforce monotonicity via the reverse cumulative minimum (np.minimum.accumulate on the
        // reversed array), so a smaller downstream value pulls earlier ones down.
        var cmin = double.PositiveInfinity;
        for (var r = m - 1; r >= 0; r--)
        {
            if (ranked[r] < cmin)
                cmin = ranked[r];
            ranked[r] = cmin;
        }

        // Clip to [0, 1] and scatter back to the original positions.
        for (var r = 0; r < m; r++)
        {
            var v = ranked[r];
            if (v < 0.0)
                v = 0.0;
            else if (v > 1.0)
                v = 1.0;
            adj[okIdx[order[r]]] = v;
        }

        return adj;
    }

    /// <summary>
    /// Adjusted p-values under the requested <paramref name="method"/>, same order as the input.
    /// </summary>
    /// <remarks>
    /// Every method here keeps <see cref="BenjaminiHochberg"/>'s NaN policy: a NaN passes through as
    /// NaN and does not count toward the number of tests. That is deliberate and differs from
    /// statsmodels, which returns all-NaN if any input is NaN - see the type summary. Whichever
    /// method is chosen, a feature that could not be tested must not silently make every other
    /// feature's correction harsher.
    /// </remarks>
    public static double[] Adjust(ReadOnlySpan<double> pValues, MultipleTesting method) => method switch
    {
        MultipleTesting.BenjaminiHochberg => BenjaminiHochberg(pValues),
        MultipleTesting.BenjaminiYekutieli => BenjaminiYekutieli(pValues),
        MultipleTesting.Bonferroni => Bonferroni(pValues),
        MultipleTesting.Holm => Holm(pValues),
        MultipleTesting.None => Uncorrected(pValues),
        _ => throw new ArgumentOutOfRangeException(nameof(method), method, "unknown correction"),
    };

    /// <summary>
    /// Benjamini-Yekutieli: BH scaled by the harmonic number <c>c(m) = sum(1/i)</c>, which makes it
    /// valid under ARBITRARY dependence between tests rather than only positive dependence.
    /// </summary>
    /// <remarks>
    /// Worth having beside BH specifically for peptide-level work: peptides from one protein are
    /// strongly correlated, which is the case BH's assumptions do not cover and the reason
    /// peptide-level q-values are described elsewhere as anti-conservative.
    /// </remarks>
    public static double[] BenjaminiYekutieli(ReadOnlySpan<double> pValues)
    {
        var finite = 0;
        foreach (var p in pValues)
            if (!double.IsNaN(p))
                finite++;

        var c = 0.0;
        for (var i = 1; i <= finite; i++)
            c += 1.0 / i;

        var adj = BenjaminiHochberg(pValues);
        for (var i = 0; i < adj.Length; i++)
            if (!double.IsNaN(adj[i]))
                adj[i] = Math.Min(adj[i] * c, 1.0);
        return adj;
    }

    /// <summary>Bonferroni: <c>p * m</c>, clipped to 1.</summary>
    public static double[] Bonferroni(ReadOnlySpan<double> pValues)
    {
        var m = CountFinite(pValues);
        var adj = new double[pValues.Length];
        for (var i = 0; i < pValues.Length; i++)
            adj[i] = double.IsNaN(pValues[i]) ? double.NaN : Math.Min(pValues[i] * m, 1.0);
        return adj;
    }

    /// <summary>
    /// Holm step-down: the i-th smallest p is scaled by <c>m - i</c>, then made monotone by a running
    /// maximum. Uniformly more powerful than Bonferroni and controls the same family-wise error rate.
    /// </summary>
    public static double[] Holm(ReadOnlySpan<double> pValues)
    {
        var n = pValues.Length;
        var adj = new double[n];
        var okIdx = new List<int>(n);
        for (var i = 0; i < n; i++)
        {
            adj[i] = double.NaN;
            if (!double.IsNaN(pValues[i]))
                okIdx.Add(i);
        }

        var m = okIdx.Count;
        if (m == 0)
            return adj;

        var ok = new double[m];
        for (var i = 0; i < m; i++)
            ok[i] = pValues[okIdx[i]];

        var order = Stats.ArgSort(ok); // ascending, stable
        var running = 0.0;
        for (var r = 0; r < m; r++)
        {
            var v = ok[order[r]] * (m - r);
            // Step-down monotonicity: an adjusted value can never fall below one already assigned to
            // a smaller raw p.
            running = Math.Max(running, v);
            adj[okIdx[order[r]]] = Math.Min(running, 1.0);
        }

        return adj;
    }

    /// <summary>The raw p-values, copied. "None" is a real choice, not an absent one.</summary>
    public static double[] Uncorrected(ReadOnlySpan<double> pValues)
    {
        var adj = new double[pValues.Length];
        for (var i = 0; i < pValues.Length; i++)
            adj[i] = pValues[i];
        return adj;
    }

    private static int CountFinite(ReadOnlySpan<double> pValues)
    {
        var m = 0;
        foreach (var p in pValues)
            if (!double.IsNaN(p))
                m++;
        return m;
    }
}
