using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.Distributions;
using SkylinePrism.Core.Numerics;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>One feature's differential-abundance result (a volcano point / hit-table row).</summary>
public sealed record DifferentialRow(
    string FeatureId,
    double LogFc,
    double Fc,
    double AveExpr,
    double T,
    double PValue,
    double AdjPValue,
    double MeanA,
    double MeanB);

/// <summary>A covariate to adjust the contrast for (Sex, PMI, batch, ...).</summary>
public abstract class Covariate
{
    protected Covariate(string name) => Name = name;

    /// <summary>Covariate name, used to label design columns.</summary>
    public string Name { get; }

    /// <summary>
    /// Build a covariate from raw metadata strings, inferring the type the way pandas dtype inference
    /// does: numeric if every non-null value parses as an invariant-culture number (so an integer-coded
    /// batch is centered, not dummy-coded), otherwise categorical. Null entries are missing.
    /// </summary>
    public static Covariate FromMetadata(string name, string?[] values)
    {
        var anyNonNull = false;
        foreach (var v in values)
        {
            if (v is null)
                continue;
            anyNonNull = true;
            if (!double.TryParse(v, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
                return new CategoricalCovariate(name, values);
        }

        if (!anyNonNull)
            return new CategoricalCovariate(name, values);

        var nums = new double[values.Length];
        for (var i = 0; i < values.Length; i++)
            nums[i] = values[i] is null
                ? double.NaN
                : double.Parse(values[i]!, NumberStyles.Float, CultureInfo.InvariantCulture);
        return new NumericCovariate(name, nums);
    }
}

/// <summary>
/// A continuous covariate. <see cref="Values"/> is indexed by matrix column (all samples), aligned to
/// the abundance matrix; <see cref="double.NaN"/> marks a missing value. A covariate missing in any
/// selected sample is skipped, as is a constant one.
/// </summary>
public sealed class NumericCovariate : Covariate
{
    public NumericCovariate(string name, double[] values) : base(name) => Values = values;

    /// <summary>Per-column values (aligned to the abundance matrix columns).</summary>
    public double[] Values { get; }
}

/// <summary>
/// A categorical covariate, dummy-encoded with the first level dropped (like limma's factor handling).
/// <see cref="Values"/> is indexed by matrix column; <c>null</c> marks a missing value (skips the
/// whole covariate). A dummy level collinear with the group contrast is dropped.
/// </summary>
public sealed class CategoricalCovariate : Covariate
{
    public CategoricalCovariate(string name, string?[] values) : base(name) => Values = values;

    /// <summary>Per-column category labels (aligned to the abundance matrix columns).</summary>
    public string?[] Values { get; }
}

/// <summary>Result of <see cref="Differential.Run"/>: per-feature rows plus run-level summary.</summary>
public sealed class DifferentialResult
{
    internal DifferentialResult(IReadOnlyList<DifferentialRow> rows, int nA, int nB,
        int nFeaturesTotal, int nFeaturesTested, double dfResidual, double dfPrior,
        string variancePrior, IReadOnlyList<string> covariatesUsed, IReadOnlyList<string> messages,
        IReadOnlyList<string> warnings)
    {
        Rows = rows;
        NA = nA;
        NB = nB;
        NFeaturesTotal = nFeaturesTotal;
        NFeaturesTested = nFeaturesTested;
        DfResidual = dfResidual;
        DfPrior = dfPrior;
        VariancePrior = variancePrior;
        CovariatesUsed = covariatesUsed;
        Messages = messages;
        Warnings = warnings;
    }

    /// <summary>Tested features, ascending by <see cref="DifferentialRow.PValue"/>.</summary>
    public IReadOnlyList<DifferentialRow> Rows { get; }

    /// <summary>Group A / group B sample counts.</summary>
    public int NA { get; }

    /// <summary>Group A / group B sample counts.</summary>
    public int NB { get; }

    /// <summary>Total features supplied.</summary>
    public int NFeaturesTotal { get; }

    /// <summary>Features tested (complete across every selected sample).</summary>
    public int NFeaturesTested { get; }

    /// <summary>Features dropped for having a missing value in a selected sample.</summary>
    public int NFeaturesDropped => NFeaturesTotal - NFeaturesTested;

    /// <summary>Residual degrees of freedom, <c>n_samples - n_coef</c>.</summary>
    public double DfResidual { get; }

    /// <summary>Empirical-Bayes prior degrees of freedom. May be <see cref="double.PositiveInfinity"/>.</summary>
    public double DfPrior { get; }

    /// <summary>Which variance prior was fitted (<c>global</c>, <c>intensity_trend</c>, ...).</summary>
    public string VariancePrior { get; }

    /// <summary>Design columns actually used, excluding the intercept and group term.</summary>
    public IReadOnlyList<string> CovariatesUsed { get; }

    /// <summary>Notes about covariates skipped or levels dropped during design construction.</summary>
    public IReadOnlyList<string> Messages { get; }

    /// <summary>Diagnostic messages from the variance-moderation step.</summary>
    public IReadOnlyList<string> Warnings { get; }
}

/// <summary>
/// Two-group differential abundance via a limma empirical-Bayes moderated t-test, ported from the
/// PRISM Differential Explorer's <c>differential</c> (prism_diff_explorer.py). Fits the shared design
/// [intercept, groupB] per feature (<see cref="LinearModel"/>), moderates the residual variances
/// (<see cref="EmpiricalBayes"/>, under the prior <see cref="DifferentialOptions.Prior"/> selects), and
/// reports the moderated t, its two-sided p-value on <c>df_residual + df_prior</c> degrees of freedom,
/// and the adjusted p-value under <see cref="DifferentialOptions.Correction"/>. B is the treatment arm,
/// so a positive <c>logFC</c> is higher in B.
///
/// <para>The simple tests (<see cref="SimpleTests"/>) and the paired designs route through the same
/// <see cref="Run(double[,], IReadOnlyList{string}, IReadOnlyList{int}, IReadOnlyList{int}, DifferentialOptions)"/>
/// entry point, so every caller gets the same row shape whichever estimator ran.</para>
/// </summary>
public static class Differential
{
    /// <summary>
    /// Run the moderated-t contrast on a LOG2 abundance matrix. <paramref name="groupAColumns"/> and
    /// <paramref name="groupBColumns"/> are disjoint column indices into
    /// <paramref name="exprLog2FeaturesBySamples"/>; only features observed (non-NaN) in every selected
    /// sample are tested. Optional <paramref name="covariates"/> adjust the contrast; the intensity-trend
    /// variance prior is not yet supported.
    /// </summary>
    public static DifferentialResult Run(
        double[,] exprLog2FeaturesBySamples,
        IReadOnlyList<string> featureIds,
        IReadOnlyList<int> groupAColumns,
        IReadOnlyList<int> groupBColumns,
        int minPerGroup = 2,
        IReadOnlyList<Covariate>? covariates = null,
        bool trend = false)
        => Run(exprLog2FeaturesBySamples, featureIds, groupAColumns, groupBColumns,
            new DifferentialOptions
            {
                MinPerGroup = minPerGroup,
                Covariates = covariates,
                // The historical meaning of this flag: limma's trend, not the toolkit's. Callers that
                // want the lab's default prior pass DifferentialOptions instead.
                Prior = trend ? VariancePrior.LimmaTrend : VariancePrior.Global,
            });

    /// <summary>
    /// The moderated-t contrast under an explicit <see cref="DifferentialOptions"/>.
    /// </summary>
    public static DifferentialResult Run(
        double[,] exprLog2FeaturesBySamples,
        IReadOnlyList<string> featureIds,
        IReadOnlyList<int> groupAColumns,
        IReadOnlyList<int> groupBColumns,
        DifferentialOptions options)
    {
        var minPerGroup = options.MinPerGroup;
        var covariates = options.Covariates;
        var pairingMessages = new List<string>();
        IReadOnlyList<SamplePair>? pairs = null;

        if (options.Design == DifferentialDesign.Paired)
        {
            if (options.SubjectLabels is null)
                throw new ArgumentException(
                    "A paired design needs a pairing column (DifferentialOptions.SubjectLabels).");

            var resolved = PairedSamples.Resolve(options.SubjectLabels, groupAColumns, groupBColumns);
            pairs = resolved.Pairs;
            pairingMessages.AddRange(resolved.Messages);
            if (pairs.Count < minPerGroup)
                throw new ArgumentException(
                    $"A paired design needs at least {minPerGroup} matched subjects; {pairs.Count} "
                    + "matched. " + string.Join(" ", resolved.Messages));

            // From here the arms ARE the matched pairs, in a common subject order. Everything
            // downstream - the design's subject block, the complete-case filter, the per-group means
            // - then lines up by position, and an unmatched sample cannot leak into one arm only.
            groupAColumns = pairs.Select(pair => pair.AColumn).ToList();
            groupBColumns = pairs.Select(pair => pair.BColumn).ToList();
        }

        if (options.Test is DifferentialTest.PairedT or DifferentialTest.Wilcoxon)
        {
            if (pairs is null)
                throw new ArgumentException(
                    $"{options.Test} is a paired test and needs the Paired design.");
            return SimpleTests.RunPaired(exprLog2FeaturesBySamples, featureIds, pairs, options.Test,
                options.Correction, minPerGroup, Messages(options, pairingMessages));
        }

        if (options.Test != DifferentialTest.ModeratedT)
            return RunSimple(exprLog2FeaturesBySamples, featureIds, groupAColumns, groupBColumns,
                options, pairingMessages);
        var nFeatures = exprLog2FeaturesBySamples.GetLength(0);
        var nColumns = exprLog2FeaturesBySamples.GetLength(1);
        if (featureIds.Count != nFeatures)
            throw new ArgumentException(
                $"featureIds has {featureIds.Count} entries but the matrix has {nFeatures} features.",
                nameof(featureIds));

        var nA = groupAColumns.Count;
        var nB = groupBColumns.Count;
        if (nA < minPerGroup || nB < minPerGroup)
            throw new ArgumentException(
                $"Each group needs at least {minPerGroup} samples (A has {nA}, B has {nB}).");

        var seen = new HashSet<int>();
        foreach (var c in groupAColumns.Concat(groupBColumns))
        {
            if (c < 0 || c >= nColumns)
                throw new ArgumentException($"Sample column index {c} is out of range.");
            if (!seen.Add(c))
                throw new ArgumentException("Groups A and B must be disjoint (a sample is in both).");
        }

        // Selected columns in order [A..., B...]; design is [intercept, groupB].
        var cols = new int[nA + nB];
        for (var i = 0; i < nA; i++)
            cols[i] = groupAColumns[i];
        for (var i = 0; i < nB; i++)
            cols[nA + i] = groupBColumns[i];
        var nSamples = cols.Length;

        var (design, covariatesUsed, messages) = BuildDesign(nA, nB, cols, covariates);
        messages.InsertRange(0, pairingMessages);
        if (pairs is not null)
        {
            // Sex, age, genotype, diagnosis - the most natural things to tick in a paired design -
            // are all constant within a subject, and therefore exactly collinear with the subject
            // block. Left in, the design is rank-deficient and the run dies on a check that names
            // neither the covariate nor a block the user never asked for. Dropping them here, and
            // saying so, matches how BuildDesign already handles a dummy collinear with the group.
            // Note this is not a limitation of the implementation: a within-subject contrast cannot
            // estimate a between-subject effect, because the subject block has already absorbed it.
            (design, covariatesUsed) =
                DropSubjectCollinear(design, covariatesUsed, pairs.Count, messages);
            design = WithSubjectBlock(design, pairs.Count);
        }
        const int coefIdx = 1; // groupB is the second design column
        var nParams = design.GetLength(1);
        if (nSamples - nParams < 1)
            throw new ArgumentException(
                $"Not enough residual degrees of freedom (n={nSamples}, params={nParams}). " +
                "Use more samples or fewer covariates.");
        if (LinAlg.MatrixRank(design) < nParams)
            throw new ArgumentException(
                "Design matrix is rank-deficient (covariates collinear with each other or with group).");

        // Complete-case filter: a feature is tested only if observed in every selected sample.
        var tested = new List<int>(nFeatures);
        for (var f = 0; f < nFeatures; f++)
        {
            var complete = true;
            for (var s = 0; s < nSamples; s++)
                if (double.IsNaN(exprLog2FeaturesBySamples[f, cols[s]]))
                {
                    complete = false;
                    break;
                }

            if (complete)
                tested.Add(f);
        }

        if (tested.Count == 0)
            throw new InvalidOperationException("No feature is observed across every selected sample.");

        var nTested = tested.Count;
        var mk = new double[nTested, nSamples];
        for (var i = 0; i < nTested; i++)
        for (var s = 0; s < nSamples; s++)
            mk[i, s] = exprLog2FeaturesBySamples[tested[i], cols[s]];

        var fit = LinearModel.Fit(mk, design);
        var variances = new double[nTested];
        for (var i = 0; i < nTested; i++)
            variances[i] = fit.Sigma[i] * fit.Sigma[i];

        var (squeezed, variancePrior) = FitPrior(
            options, exprLog2FeaturesBySamples, cols, nA, nB, variances, fit, tested, messages);

        var dfTotal = fit.DfResidual + squeezed.DfPrior;
        var stdevUnscaled = fit.StdevUnscaled[coefIdx];

        var pValues = new double[nTested];
        var rows = new DifferentialRow[nTested];
        var groupB = new double[nB];
        var groupA = new double[nA];
        for (var i = 0; i < nTested; i++)
        {
            var coef = fit.Coefficients[i, coefIdx];
            var t = coef / (stdevUnscaled * Math.Sqrt(squeezed.VarPost[i]));
            var p = ModeratedPValue(t, dfTotal);
            pValues[i] = p;

            for (var s = 0; s < nA; s++)
                groupA[s] = mk[i, s];
            for (var s = 0; s < nB; s++)
                groupB[s] = mk[i, nA + s];

            rows[i] = new DifferentialRow(
                featureIds[tested[i]], coef, Math.Pow(2.0, coef), fit.Amean[i], t, p,
                double.NaN, NumpyMath.Mean(groupA), NumpyMath.Mean(groupB));
        }

        var adj = Fdr.Adjust(pValues, options.Correction);
        for (var i = 0; i < nTested; i++)
            rows[i] = rows[i] with { AdjPValue = adj[i] };

        var ordered = rows.OrderBy(r => r.PValue).ToArray();

        return new DifferentialResult(ordered, nA, nB, nFeatures, nTested, fit.DfResidual,
            squeezed.DfPrior, variancePrior, covariatesUsed, messages, squeezed.Warnings);
    }

    /// <summary>
    /// Build the design matrix [intercept, groupB, covariate columns...]. Numeric covariates are
    /// mean-centered; categorical covariates are dummy-encoded with the first (sorted) level dropped.
    /// A covariate missing in any selected sample, or constant, is skipped; a dummy level collinear
    /// with the group is dropped. Every skip/drop is recorded in the returned messages.
    /// </summary>
    /// <summary>
    /// Drop covariate columns that the subject block will make redundant, naming them.
    /// </summary>
    /// <remarks>
    /// A covariate constant within every subject carries no information a within-subject contrast
    /// can use - the subject block absorbs it entirely - so its column is exactly a linear
    /// combination of the block. Detected by the definition rather than by a rank test, because the
    /// definition is what can be explained to a reader.
    /// </remarks>
    private static (double[,] Design, List<string> Used) DropSubjectCollinear(
        double[,] design, List<string> covariatesUsed, int nPairs, List<string> messages)
    {
        var nSamples = design.GetLength(0);
        var nParams = design.GetLength(1);
        const int firstCovariate = 2; // [intercept, group, covariates...]
        if (nParams <= firstCovariate || nPairs < 1)
            return (design, covariatesUsed);

        var keep = new List<int>();
        var dropped = new List<string>();
        for (var c = 0; c < nParams; c++)
        {
            if (c < firstCovariate)
            {
                keep.Add(c);
                continue;
            }

            var constantWithinSubject = true;
            for (var j = 0; j < nPairs && constantWithinSubject; j++)
                if (design[j, c] != design[nPairs + j, c])
                    constantWithinSubject = false;

            if (constantWithinSubject)
            {
                var name = c - firstCovariate < covariatesUsed.Count
                    ? covariatesUsed[c - firstCovariate]
                    : $"column {c}";
                dropped.Add(name);
            }
            else
            {
                keep.Add(c);
            }
        }

        if (dropped.Count == 0)
            return (design, covariatesUsed);

        messages.Add($"Covariate(s) {string.Join(", ", dropped)} are constant within each subject, so "
            + "the paired design's subject block already accounts for them - they were dropped. A "
            + "within-subject contrast cannot estimate a between-subject effect.");

        var reduced = new double[nSamples, keep.Count];
        for (var r = 0; r < nSamples; r++)
            for (var k = 0; k < keep.Count; k++)
                reduced[r, k] = design[r, keep[k]];

        var used = covariatesUsed.Where(n => !dropped.Contains(n)).ToList();
        return (reduced, used);
    }

    /// <summary>
    /// Append a fixed-effect subject block to a paired design: one indicator per subject after the
    /// first, each marking that subject's two samples.
    /// </summary>
    /// <remarks>
    /// <para>This is what makes the contrast a WITHIN-subject one. Without it every subject's overall
    /// level is part of the residual, so between-subject variation - which a paired design exists to
    /// remove - inflates the variance and costs the test its power.</para>
    /// <para>A FIXED effect, not a random one, matching the lab's toolkit
    /// (<c>statistical_analysis.py:1321</c>). A random intercept would be a mixed model, which is a
    /// different estimator and is not implemented. The first subject is dropped, as always with
    /// indicator coding, because the intercept already spans it.</para>
    /// <para>The columns are positional: the caller has already reduced the arms to matched pairs in
    /// a common subject order, so subject j is at row j in arm A and row <c>nPairs + j</c> in arm B.</para>
    /// </remarks>
    private static double[,] WithSubjectBlock(double[,] design, int nPairs)
    {
        var nSamples = design.GetLength(0);
        var nParams = design.GetLength(1);
        var extra = nPairs - 1;
        if (extra <= 0)
            return design;

        var widened = new double[nSamples, nParams + extra];
        for (var r = 0; r < nSamples; r++)
        for (var c = 0; c < nParams; c++)
            widened[r, c] = design[r, c];

        for (var j = 1; j < nPairs; j++)
        {
            widened[j, nParams + j - 1] = 1.0;          // arm A half
            widened[nPairs + j, nParams + j - 1] = 1.0; // arm B half
        }

        return widened;
    }

    /// <summary>Pairing notes plus whatever the options themselves make unusable.</summary>
    private static List<string> Messages(DifferentialOptions options, List<string> pairingMessages)
    {
        var messages = new List<string>(pairingMessages);
        if (options.Covariates is { Count: > 0 } cov)
            messages.Add($"{options.Describe()} cannot adjust for covariates "
                + $"({string.Join(", ", cov.Select(c => c.Name))}) - it has no design matrix to put "
                + "them in. Use the moderated t for an adjusted contrast.");
        return messages;
    }

    private static (double[,] Design, List<string> CovariatesUsed, List<string> Messages) BuildDesign(
        int nA, int nB, int[] cols, IReadOnlyList<Covariate>? covariates)
    {
        var nSamples = nA + nB;
        var grp = new double[nSamples];
        for (var s = nA; s < nSamples; s++)
            grp[s] = 1.0;

        var extra = new List<double[]>();
        var names = new List<string>();
        var messages = new List<string>();

        if (covariates != null)
        {
            foreach (var cov in covariates)
            {
                if (cov is NumericCovariate num)
                {
                    var v = new double[nSamples];
                    var missing = false;
                    for (var s = 0; s < nSamples; s++)
                    {
                        v[s] = num.Values[cols[s]];
                        if (double.IsNaN(v[s]))
                            missing = true;
                    }

                    if (missing)
                    {
                        messages.Add($"Covariate '{num.Name}' has missing values in selected samples - skipped.");
                        continue;
                    }

                    if (Distinct(v) < 2)
                    {
                        messages.Add($"Covariate '{num.Name}' is constant - skipped.");
                        continue;
                    }

                    var mean = NumpyMath.Mean(v);
                    var centered = new double[nSamples];
                    for (var s = 0; s < nSamples; s++)
                        centered[s] = v[s] - mean;
                    extra.Add(centered);
                    names.Add(num.Name);
                }
                else if (cov is CategoricalCovariate cat)
                {
                    var missing = false;
                    for (var s = 0; s < nSamples; s++)
                        if (cat.Values[cols[s]] is null)
                        {
                            missing = true;
                            break;
                        }

                    if (missing)
                    {
                        messages.Add($"Covariate '{cat.Name}' has missing values in selected samples - skipped.");
                        continue;
                    }

                    var vals = new string[nSamples];
                    for (var s = 0; s < nSamples; s++)
                        vals[s] = cat.Values[cols[s]]!;

                    // get_dummies(drop_first=True): sorted levels, drop the first, one indicator each.
                    var levels = vals.Distinct().OrderBy(x => x, StringComparer.Ordinal).ToArray();
                    for (var li = 1; li < levels.Length; li++)
                    {
                        var level = levels[li];
                        var col = new double[nSamples];
                        for (var s = 0; s < nSamples; s++)
                            col[s] = vals[s] == level ? 1.0 : 0.0;

                        if (Distinct(col) < 2)
                            continue;
                        if (AllClose(col, grp, complement: false) || AllClose(col, grp, complement: true))
                        {
                            messages.Add($"Covariate level '{cat.Name}_{level}' is confounded with group - dropped.");
                            continue;
                        }

                        extra.Add(col);
                        names.Add($"{cat.Name}_{level}");
                    }
                }
            }
        }

        var nCoef = 2 + extra.Count;
        var design = new double[nSamples, nCoef];
        for (var s = 0; s < nSamples; s++)
        {
            design[s, 0] = 1.0;
            design[s, 1] = grp[s];
            for (var c = 0; c < extra.Count; c++)
                design[s, 2 + c] = extra[c][s];
        }

        return (design, names, messages);
    }

    /// <summary>Number of exactly-distinct values (numpy.unique semantics), enough to spot a constant.</summary>
    private static int Distinct(double[] v)
    {
        var set = new HashSet<double>();
        foreach (var x in v)
            set.Add(x);
        return set.Count;
    }

    /// <summary>
    /// numpy.allclose(a, complement ? 1-grp : grp) with numpy's default tolerances (rtol 1e-5, atol
    /// 1e-8): true when the dummy column exactly tracks (or complements) the group assignment.
    /// </summary>
    private static bool AllClose(double[] a, double[] grp, bool complement)
    {
        for (var i = 0; i < a.Length; i++)
        {
            var b = complement ? 1.0 - grp[i] : grp[i];
            if (Math.Abs(a[i] - b) > 1e-8 + 1e-5 * Math.Abs(b))
                return false;
        }

        return true;
    }

    /// <summary>
    /// Two-sided moderated-t p-value: <c>2 * cdf(-|t|)</c> on <paramref name="dfTotal"/> degrees of
    /// freedom. An infinite prior gives an infinite total df, where the t-distribution is the standard
    /// normal (matching scipy's <c>t.cdf(df=inf)</c>).
    /// </summary>
    /// <summary>
    /// The non-moderated tests, which share this entry point but none of the linear-model machinery.
    /// </summary>
    /// <remarks>
    /// Covariates are reported as ignored rather than dropped in silence: a t-test has no design
    /// matrix to hold one, so a ticked covariate simply cannot be honored, and a reader who ticked
    /// one is entitled to know the contrast in front of them is unadjusted.
    /// </remarks>
    private static DifferentialResult RunSimple(
        double[,] exprLog2FeaturesBySamples,
        IReadOnlyList<string> featureIds,
        IReadOnlyList<int> groupAColumns,
        IReadOnlyList<int> groupBColumns,
        DifferentialOptions options,
        IReadOnlyList<string> pairingMessages)
    {
        var messages = new List<string>(pairingMessages);
        if (options.Covariates is { Count: > 0 } cov)
            messages.Add($"{options.Describe()} cannot adjust for covariates "
                + $"({string.Join(", ", cov.Select(c => c.Name))}) - it has no design matrix to put "
                + "them in. Use the moderated t for an adjusted contrast.");

        return SimpleTests.Run(exprLog2FeaturesBySamples, featureIds, groupAColumns, groupBColumns,
            options.Test, options.Correction, options.MinPerGroup, messages);
    }

    /// <summary>
    /// The variance prior the options asked for, with the name to report it under. Falls back to the
    /// global prior - saying so in <paramref name="messages"/> - whenever the requested one cannot be
    /// fitted, because a silently substituted prior is a silently different p-value.
    /// </summary>
    private static (SqueezeVarResult Squeezed, string Name) FitPrior(
        DifferentialOptions options, double[,] expr, int[] cols, int nA, int nB, double[] variances,
        LinearModelFit fit, List<int> tested, List<string> messages)
    {
        // Peptide counts arrive per FEATURE of the input matrix; the fit is over the tested subset,
        // so they have to be gathered in the same order or the trend would pair each variance with
        // another feature's count.
        double[]? counts = null;
        if (options.PeptideCounts is { } supplied)
        {
            counts = new double[tested.Count];
            for (var i = 0; i < tested.Count; i++)
                counts[i] = tested[i] < supplied.Count ? supplied[tested[i]] : double.NaN;
        }

        var needsCounts = options.Prior is VariancePrior.PeptideCount;
        if (needsCounts && counts is null)
        {
            messages.Add("This prior needs a peptide count per feature, which only the protein-level "
                + "matrix carries - the global prior was used instead.");
            return (EmpiricalBayes.SqueezeVarGlobal(variances, fit.DfResidual), "global");
        }

        switch (options.Prior)
        {
            case VariancePrior.LimmaTrend when fit.Amean.All(double.IsFinite):
                return (EmpiricalBayes.SqueezeVarTrend(variances, fit.DfResidual, fit.Amean),
                    "limma-trend");

            case VariancePrior.LimmaTrend:
                messages.Add("Mean expression has non-finite values - the limma-trend prior is "
                    + "unavailable, so the global prior was used.");
                break;

            case VariancePrior.IntensityTrend:
            {
                var prior = VariancePriors.IntensityTrend(
                    expr, tested, PriorGroups(options, cols, nA, nB));
                if (prior is not null)
                    return (WithGlobalDf(variances, fit.DfResidual, prior), "intensity-trend");

                messages.Add("Too few usable (feature, group) points to fit the intensity trend - "
                    + "the global prior was used instead.");
                break;
            }

            case VariancePrior.PeptideCount:
            {
                var prior = VariancePriors.PeptideCountTrend(variances, counts!);
                if (prior is not null)
                    return (WithGlobalDf(variances, fit.DfResidual, prior), "peptide-count");
                messages.Add("Too few features carry a usable peptide count to fit the trend - the "
                    + "global prior was used instead.");
                break;
            }
        }

        return (EmpiricalBayes.SqueezeVarGlobal(variances, fit.DfResidual), "global");
    }

    /// <summary>
    /// The groups the variance prior is fitted over: the two contrast arms by default, or whichever
    /// replicates the caller nominated.
    /// </summary>
    /// <remarks>
    /// The override exists because a design group's within-group spread contains inter-subject
    /// BIOLOGY, which inflates the prior and over-shrinks real signal. Pointing it at dedicated QC or
    /// reference injections measures instrument-and-workflow variance instead, which is what the
    /// prior is supposed to describe - and those replicates take no part in the contrast, so their
    /// columns exist only in the FULL matrix, which is why these indices are absolute.
    /// </remarks>
    private static IReadOnlyList<IReadOnlyList<int>> PriorGroups(
        DifferentialOptions options, int[] cols, int nA, int nB) =>
        options.PriorGroupColumns is null
            ? new IReadOnlyList<int>[]
            {
                cols.Take(nA).ToArray(),
                cols.Skip(nA).Take(nB).ToArray(),
            }
            : options.PriorGroupColumns;

    /// <summary>
    /// A per-feature prior scale paired with the GLOBAL prior degrees of freedom.
    /// </summary>
    /// <remarks>
    /// Every toolkit-style prior works this way - it fits a scale and leaves the degrees of freedom
    /// alone - and that is exactly what separates them from limma's trend, which re-estimates both.
    /// Kept in one place so a new prior cannot accidentally re-estimate the df and still call itself
    /// one of these.
    /// </remarks>
    private static SqueezeVarResult WithGlobalDf(double[] variances, double dfResidual, double[] prior)
    {
        var global = EmpiricalBayes.SqueezeVarGlobal(variances, dfResidual);
        return EmpiricalBayes.SqueezeVarWithScale(
            variances, dfResidual, prior, global.DfPrior, global.Warnings);
    }

    private static double ModeratedPValue(double t, double dfTotal)
    {
        return Distributions.TwoSidedT(t, dfTotal);
    }
}
