using System.Collections.Generic;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>How the samples are related to each other - the design the contrast is estimated under.</summary>
public enum DifferentialDesign
{
    /// <summary>Two independent groups.</summary>
    Unpaired,

    /// <summary>
    /// Two measurements of the same subject. Fitted as a FIXED-effect subject block
    /// (<c>[1, treat, subject one-hot]</c>), which is what the lab's toolkit does - not a random
    /// effect. A random-intercept model is a separate thing and is not implemented yet.
    /// </summary>
    Paired,

    /// <summary>A slope against a numeric time or dose, testing whether it differs from zero.</summary>
    LinearTrend,
}

/// <summary>The estimator applied to the design.</summary>
public enum DifferentialTest
{
    /// <summary>limma empirical-Bayes moderated t. The only test that uses a variance prior.</summary>
    ModeratedT,

    /// <summary>Two-sample t with unequal variances (Welch-Satterthwaite df).</summary>
    WelchT,

    /// <summary>Two-sample t with a pooled variance.</summary>
    StudentT,

    /// <summary>One-sample t on the within-subject differences.</summary>
    PairedT,

    /// <summary>Wilcoxon signed-rank on the within-subject differences.</summary>
    Wilcoxon,

    /// <summary>Mann-Whitney U, two-sided.</summary>
    MannWhitney,
}

/// <summary>
/// Where the moderated-t prior variance comes from. Only meaningful for
/// <see cref="DifferentialTest.ModeratedT"/>.
/// </summary>
public enum VariancePrior
{
    /// <summary>One global <c>(s0^2, d0)</c> for every feature - Smyth (2004).</summary>
    Global,

    /// <summary>
    /// The lab's default: a LOWESS of log(within-group variance) on log(within-group MEAN INTENSITY),
    /// fitted on RAW pre-log intensities, one point per (feature, group), converted back to log space
    /// by the delta method. Only the per-feature scale is replaced; the prior degrees of freedom stay
    /// global.
    /// </summary>
    /// <remarks>
    /// This is <c>proteomics-toolkit</c>'s <c>moderation="intensity_trend"</c>, and it is NOT the same
    /// estimator as <see cref="LimmaTrend"/> despite both being called an intensity trend. See
    /// <c>docs/differential-analysis.md</c>: they disagree by a median 3-7% on p-values and up to 178%
    /// on individual features.
    /// </remarks>
    IntensityTrend,

    /// <summary>
    /// limma's own <c>trend=TRUE</c>: a natural cubic spline of log(residual variance) on mean LOG2
    /// expression, re-estimating the prior degrees of freedom as well as the scale. Kept because it is
    /// what limma does and what the committed inmoose goldens pin, but it is not the lab's default.
    /// </summary>
    LimmaTrend,

    /// <summary>
    /// DEqMS (Zhu 2020): a LOWESS of log(residual variance) on log(peptide count). Needs a peptide
    /// count per feature, so it is protein-level only.
    /// </summary>
    PeptideCount,

    /// <summary>
    /// Additive two-stage: <c>log(var) = f1(log mean intensity) + f2(log peptide count)</c>, the second
    /// LOWESS fitted on the first stage's residual. Protein-level only.
    /// </summary>
    IntensityAndPeptideCount,
}

/// <summary>How the raw p-values are corrected for multiple testing.</summary>
public enum MultipleTesting
{
    /// <summary>Benjamini-Hochberg step-up FDR.</summary>
    BenjaminiHochberg,

    /// <summary>Benjamini-Yekutieli, valid under arbitrary dependence.</summary>
    BenjaminiYekutieli,

    /// <summary>Bonferroni family-wise correction.</summary>
    Bonferroni,

    /// <summary>Holm step-down family-wise correction.</summary>
    Holm,

    /// <summary>No correction - the adjusted column repeats the raw p.</summary>
    None,
}

/// <summary>
/// Everything that selects WHICH analysis <see cref="Differential.Run"/> performs, as opposed to the
/// data it performs it on.
/// </summary>
/// <remarks>
/// <para>Three orthogonal axes, which is how the lab's toolkit models it too: a
/// <see cref="Design"/>, a <see cref="Test"/>, and - for the moderated test alone - a
/// <see cref="Prior"/>. Kept as one record rather than a widening parameter list because the
/// combinations are what a caller reasons about, and because a record can gain an axis without
/// breaking every call site.</para>
/// <para>The defaults are the lab's defaults: an unpaired moderated-t with the intensity-trend prior.
/// Note that <see cref="Prior"/> defaults to <see cref="VariancePrior.IntensityTrend"/> rather than
/// <see cref="VariancePrior.Global"/>, which is a deliberate change from the old <c>trend: false</c>
/// default - see the release notes.</para>
/// </remarks>
public sealed record DifferentialOptions
{
    /// <summary>The lab's defaults.</summary>
    public static readonly DifferentialOptions Default = new();

    /// <summary>The unpaired global-prior moderated t, which is what <c>trend: false</c> used to mean.</summary>
    public static readonly DifferentialOptions GlobalPrior = new() { Prior = VariancePrior.Global };

    public DifferentialDesign Design { get; init; } = DifferentialDesign.Unpaired;

    public DifferentialTest Test { get; init; } = DifferentialTest.ModeratedT;

    public VariancePrior Prior { get; init; } = VariancePrior.IntensityTrend;

    public MultipleTesting Correction { get; init; } = MultipleTesting.BenjaminiHochberg;

    /// <summary>
    /// Pairing key per sample column, for <see cref="DifferentialDesign.Paired"/> and for the paired
    /// tests. Indexed by the ORIGINAL matrix column, so it is read through the selected columns.
    /// </summary>
    public IReadOnlyList<string?>? SubjectLabels { get; init; }

    /// <summary>Numeric time or dose per sample column, for <see cref="DifferentialDesign.LinearTrend"/>.</summary>
    public IReadOnlyList<double>? TimeValues { get; init; }

    /// <summary>
    /// Peptides behind each FEATURE (not sample), for the count-based priors. Parallel to the feature
    /// ids. Null where the matrix does not carry <c>n_peptides</c>, which is the peptide-level case.
    /// </summary>
    public IReadOnlyList<double>? PeptideCounts { get; init; }

    /// <summary>
    /// Sample columns to fit the variance prior on, instead of the contrast arms - the toolkit's
    /// <c>variance_prior_group_column</c>. Pointing this at dedicated QC or reference replicates keeps
    /// inter-subject biology out of the prior, which otherwise inflates it and over-shrinks real
    /// signal. Null uses the contrast arms.
    /// </summary>
    public IReadOnlyList<int>? PriorGroupColumns { get; init; }

    /// <summary>Minimum samples per arm before the contrast is refused.</summary>
    public int MinPerGroup { get; init; } = 2;

    /// <summary>Covariates to adjust the contrast for. Ignored by the rank-based tests.</summary>
    public IReadOnlyList<Covariate>? Covariates { get; init; }

    /// <summary>A short human-readable name for what this asks for, for the status line and provenance.</summary>
    public string Describe()
    {
        var test = Test switch
        {
            DifferentialTest.ModeratedT => "moderated t (" + Prior switch
            {
                VariancePrior.Global => "global prior",
                VariancePrior.IntensityTrend => "intensity-trend prior",
                VariancePrior.LimmaTrend => "limma-trend prior",
                VariancePrior.PeptideCount => "peptide-count prior",
                _ => "intensity + peptide-count prior",
            } + ")",
            DifferentialTest.WelchT => "Welch t",
            DifferentialTest.StudentT => "Student t",
            DifferentialTest.PairedT => "paired t",
            DifferentialTest.Wilcoxon => "Wilcoxon signed-rank",
            _ => "Mann-Whitney U",
        };
        var design = Design switch
        {
            DifferentialDesign.Paired => "paired",
            DifferentialDesign.LinearTrend => "linear trend",
            _ => "unpaired",
        };
        return $"{test}, {design}";
    }
}
