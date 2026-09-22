using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SkylinePrism.Core.Numerics;

/// <summary>How each feature is rescaled before the decomposition.</summary>
public enum PcaScaling
{
    /// <summary>
    /// Center and divide by the population standard deviation, so every feature contributes
    /// equally (sklearn's <c>StandardScaler</c> + <c>PCA</c>). What the QC scatter wants: it is
    /// looking for sample groupings, and without scaling the handful of most abundant proteins
    /// would set the axes on their own.
    /// </summary>
    Standardize,

    /// <summary>
    /// Center only, leaving each feature at its own scale (<c>numpy.linalg.svd</c> of the centered
    /// matrix). What the differential explorer wants: on log2 abundances the scale carries real
    /// information about which features vary, and standardizing would amplify the flat ones.
    /// </summary>
    CenterOnly,
}

/// <summary>What to do about a feature that is not observed in every sample.</summary>
public enum PcaMissingPolicy
{
    /// <summary>
    /// Impute a missing cell to its feature's mean over the observed samples (so 0 after
    /// centering) and keep the feature, as long as at least two samples observed it. Necessary on
    /// the QC path: with hundreds of samples, dropping every feature with a single gap discards
    /// almost the whole matrix.
    /// </summary>
    ImputeFeatureMean,

    /// <summary>
    /// Drop any feature not observed in every selected sample. Stricter, and the right choice when
    /// the sample set is small and chosen deliberately, where an imputed cell would be a
    /// meaningful fraction of the evidence.
    /// </summary>
    CompleteCase,
}

/// <summary>Options for <see cref="Pca.Fit"/>.</summary>
public sealed class PcaOptions
{
    /// <summary>Components to return. Capped by the sample count and the usable feature count.</summary>
    public int Components { get; init; } = 2;

    /// <summary>Per-feature rescaling. See <see cref="PcaScaling"/>.</summary>
    public PcaScaling Scaling { get; init; } = PcaScaling.Standardize;

    /// <summary>Treatment of unobserved cells. See <see cref="PcaMissingPolicy"/>.</summary>
    public PcaMissingPolicy Missing { get; init; } = PcaMissingPolicy.ImputeFeatureMean;

    /// <summary>
    /// Sample column indices to include, or null for all of them. Statistics are computed over the
    /// selected samples only - a feature's mean is its mean across THESE samples, not across the
    /// cohort - so a subset PCA is a PCA of the subset, not a projection of one.
    /// </summary>
    public IReadOnlyList<int>? SampleColumns { get; init; }

    /// <summary>
    /// Sample ids for EVERY column of the matrix (not just the selected ones), or null. When
    /// given, <see cref="PcaResult.SampleIds"/> comes back holding the selected ids in score-row
    /// order, which is what a plot needs to label or hover a point.
    /// </summary>
    public IReadOnlyList<string>? SampleIds { get; init; }

    /// <summary>
    /// Throw when there is too little data to decompose, instead of returning zero scores.
    ///
    /// <para>The two callers genuinely want opposite things here, which is why this is a flag and
    /// not a decision. A QC report is generated unattended and must not fail the run over one
    /// unplottable panel, so it takes the zeros and draws an empty frame. An interactive pane has
    /// someone sitting in front of it who needs to be told which condition was not met, so it
    /// takes the exception and puts the message in the status line.</para>
    /// </summary>
    public bool RequireSufficientData { get; init; }
}

/// <summary>Result of <see cref="Pca.Fit"/>.</summary>
public sealed class PcaResult
{
    internal PcaResult(double[,] scores, double[] varianceRatio, int nFeaturesUsed,
        string[] sampleIds)
    {
        Scores = scores;
        VarianceRatio = varianceRatio;
        NFeaturesUsed = nFeaturesUsed;
        SampleIds = sampleIds;
    }

    /// <summary>
    /// Sample ids in row order of <see cref="Scores"/>, or empty when the caller supplied none.
    /// </summary>
    public string[] SampleIds { get; }

    /// <summary>
    /// Principal-coordinate scores <c>U * S</c>, indexed <c>[sample, component]</c>, with samples
    /// in the order given by <see cref="PcaOptions.SampleColumns"/>. Each component's SIGN is
    /// arbitrary - an eigenvector and its negation are the same component - so compare magnitudes,
    /// never signs.
    /// </summary>
    public double[,] Scores { get; }

    /// <summary>Fraction of total variance carried by each returned component.</summary>
    public double[] VarianceRatio { get; }

    /// <summary>Features that survived the scaling and missing-value rules.</summary>
    public int NFeaturesUsed { get; }
}

/// <summary>
/// Sample-space PCA: the one implementation behind every PCA PRISM draws - the QC scatter, the
/// dual-control separation metric in <see cref="Qc.ValidationStatus"/>, and the differential
/// explorer's sample plot.
///
/// <para>There were two of these until the QC plot and the differential plot were merged. They
/// differed in what a caller could ask for, not in what PCA is, so the differences became
/// <see cref="PcaOptions"/> rather than a second copy: standardize vs center-only, impute vs
/// complete-case, two components vs k with their variance ratios, all samples vs a chosen subset.
/// Every combination produces exactly the numbers the corresponding implementation produced
/// before, which is what <c>PcaTests</c> and the <c>pca.json</c> goldens pin.</para>
///
/// <para><b>Never form the feature-space SVD.</b> Proteomics matrices have far more features than
/// samples, and MathNet's dense SVD would allocate an (nFeatures x nFeatures) V - 40k x 40k is
/// about 13 GB, which hangs or OOMs the machine. This eigendecomposes the small
/// (nSamples x nSamples) Gram matrix <c>X Xᵀ = U S² Uᵀ</c> instead, so cost scales with the sample
/// count and not the feature count.</para>
/// </summary>
public static class Pca
{
    /// <summary>Features per Gram-accumulation block; bounds the working matrix to nSamples x this.</summary>
    private const int FeatureBlock = 4096;

    /// <summary>
    /// Compute 2-D PCA scores from a <b>[nSamples, nFeatures]</b> matrix, standardizing features
    /// and imputing missing cells to the feature mean. The QC scatter's long-standing behavior.
    /// </summary>
    public static double[,] Fit2D(double[,] samplesByFeatures)
        => Fit(samplesByFeatures, transposed: false, new PcaOptions()).Scores;

    /// <summary>
    /// The same PCA on a <b>[nFeatures, nSamples]</b> matrix - the orientation every PRISM matrix
    /// is already in - without materializing its transpose.
    /// <para>
    /// Worth having rather than calling <c>Transpose</c> first: the transpose is a full second copy
    /// of the biggest object in the pipeline (5.7 GB on a 100-document peptide matrix), it lands on
    /// the large object heap where it inflates the working set well past the point it is dropped,
    /// and the loop below reads the matrix one FEATURE at a time anyway - which is the untransposed
    /// layout.
    /// </para>
    /// </summary>
    public static double[,] Fit2DOfFeaturesBySamples(double[,] featuresBySamples)
        => Fit(featuresBySamples, transposed: true, new PcaOptions()).Scores;

    /// <summary>
    /// Full PCA over a <b>[nFeatures, nSamples]</b> matrix - the orientation PRISM stores - with
    /// the behavior chosen by <paramref name="options"/>.
    /// </summary>
    public static PcaResult Fit(double[,] featuresBySamples, PcaOptions options)
        => Fit(featuresBySamples, transposed: true, options);

    /// <summary>
    /// Full PCA over a <b>[nSamples, nFeatures]</b> matrix.
    /// </summary>
    public static PcaResult FitOfSamplesByFeatures(double[,] samplesByFeatures, PcaOptions options)
        => Fit(samplesByFeatures, transposed: false, options);

    private static PcaResult Fit(double[,] featuresBySamples, bool transposed,
        PcaOptions options)
    {
        var matrix = featuresBySamples;
        var totalSamples = transposed ? matrix.GetLength(1) : matrix.GetLength(0);
        var nFeatures = transposed ? matrix.GetLength(0) : matrix.GetLength(1);

        int[] cols;
        if (options.SampleColumns is null)
        {
            cols = new int[totalSamples];
            for (var i = 0; i < totalSamples; i++)
                cols[i] = i;
        }
        else
        {
            cols = new int[options.SampleColumns.Count];
            for (var i = 0; i < cols.Length; i++)
            {
                var c = options.SampleColumns[i];
                if (c < 0 || c >= totalSamples)
                    throw new ArgumentException(
                        $"Sample column index {c} is out of range (matrix has {totalSamples}).",
                        nameof(options));
                cols[i] = c;
            }
        }

        var nSamples = cols.Length;
        var requested = Math.Max(options.Components, 1);

        if (options.RequireSufficientData && nSamples < 2)
            throw new ArgumentException("PCA needs at least 2 samples.", nameof(options));

        var ids = Array.Empty<string>();
        if (options.SampleIds is not null)
        {
            if (options.SampleIds.Count != totalSamples)
                throw new ArgumentException(
                    $"SampleIds has {options.SampleIds.Count} entries but the matrix has "
                    + $"{totalSamples} columns.", nameof(options));
            ids = new string[nSamples];
            for (var i = 0; i < nSamples; i++)
                ids[i] = options.SampleIds[cols[i]];
        }

        // The Gram matrix is accumulated a BLOCK OF FEATURES at a time rather than by materializing
        // the standardized matrix. X is nSamples x nFeatures, and on a 100-document cohort that is
        // 9,600 x 75,000 - 5.7 GB, which used to be built twice over (a list of per-feature arrays,
        // then a DenseMatrix copy of it) purely to produce a Gram matrix of nSamples x nSamples, a
        // few tens of MB. Blocking keeps the optimized multiply while holding only
        // nSamples x FeatureBlock at a time (~150 MB at 9,600 samples).
        var gram = new DenseMatrix(nSamples, nSamples);
        var block = new DenseMatrix(nSamples, FeatureBlock);
        var col = new double[nSamples];
        var kept = 0;
        var inBlock = 0;

        void FlushBlock()
        {
            if (inBlock == 0)
                return;
            // Xb.Xbᵀ summed over blocks IS X.Xᵀ - the product is a sum over features either way.
            var used = inBlock == FeatureBlock
                ? (Matrix<double>)block
                : block.SubMatrix(0, nSamples, 0, inBlock);
            gram.Add(used.TransposeAndMultiply(used), gram);
            inBlock = 0;
        }

        var completeCase = options.Missing == PcaMissingPolicy.CompleteCase;
        var standardize = options.Scaling == PcaScaling.Standardize;

        for (var j = 0; j < nFeatures; j++)
        {
            double sum = 0;
            var cnt = 0;
            var incomplete = false;
            for (var i = 0; i < nSamples; i++)
            {
                var v = transposed ? matrix[j, cols[i]] : matrix[cols[i], j];
                col[i] = v;
                if (double.IsNaN(v))
                    incomplete = true;
                else
                {
                    sum += v;
                    cnt++;
                }
            }

            if (completeCase)
            {
                if (incomplete)
                    continue;
            }
            else if (cnt < 2)
            {
                // Fewer than two observations says nothing about how the samples covary, and the
                // standard deviation below would be meaningless.
                continue;
            }

            var mean = sum / cnt;

            var scale = 1.0;
            if (standardize)
            {
                double ss = 0;
                for (var i = 0; i < nSamples; i++)
                    if (!double.IsNaN(col[i]))
                    {
                        var d = col[i] - mean;
                        ss += d * d;
                    }

                var std = Math.Sqrt(ss / cnt);
                if (std == 0.0)
                    continue; // a constant feature carries no information and cannot be scaled
                scale = 1.0 / std;
            }

            // Under CenterOnly a constant feature is NOT dropped: it contributes exactly zero to
            // the Gram matrix, so it changes no score, and dropping it would silently move
            // NFeaturesUsed away from the complete-case count the caller was told about.
            for (var i = 0; i < nSamples; i++)
                block[i, inBlock] = double.IsNaN(col[i]) ? 0.0 : (col[i] - mean) * scale;

            kept++;
            if (++inBlock == FeatureBlock)
                FlushBlock();
        }

        FlushBlock();

        if (kept < 2 || nSamples < 2)
        {
            if (options.RequireSufficientData)
                throw new ArgumentException(
                    $"PCA needs at least 2 features usable across the samples; {kept} were.",
                    nameof(featuresBySamples));

            // Otherwise the historical Fit2D contract: too little to decompose gives zeros rather
            // than an exception, because a QC plot with nothing in it beats a crashed report.
            return new PcaResult(new double[nSamples, Math.Min(requested, 2)],
                new double[Math.Min(requested, 2)], kept, ids);
        }

        // X = U S Vᵀ  =>  X Xᵀ = U S² Uᵀ, so the symmetric eigendecomposition of the Gram matrix
        // gives U (eigenvectors) and S² (eigenvalues); scores = U * S.
        var evd = gram.Evd(Symmetricity.Symmetric);
        var evals = evd.EigenValues;   // ascending
        var evecs = evd.EigenVectors;  // columns are the eigenvectors (= U)

        double totalVar = 0;
        for (var i = 0; i < nSamples; i++)
            totalVar += evals[i].Real;

        // numpy's svd(full_matrices=False) yields min(nSamples, nFeatures) singular values; matching
        // it avoids padding the tail with zero-rank components (which show as "-0.0%" variance).
        var k = Math.Min(requested, Math.Min(nSamples, kept));
        var scores = new double[nSamples, k];
        var varianceRatio = new double[k];
        for (var c = 0; c < k; c++)
        {
            var idx = nSamples - 1 - c; // descending
            var lambda = evals[idx].Real;
            var s = Math.Sqrt(Math.Max(lambda, 0.0));
            for (var i = 0; i < nSamples; i++)
                scores[i, c] = evecs[i, idx] * s;
            varianceRatio[c] = totalVar > 0 ? lambda / totalVar : 0.0;
        }

        return new PcaResult(scores, varianceRatio, kept, ids);
    }
}
