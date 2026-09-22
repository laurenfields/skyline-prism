using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.DifferentialAnalysis.Enrichment;
using SkylinePrism.Core.IO;
using SkylinePrism.Core.Qc;
using SkylinePrism.Core.Visualization;
using SkylinePrism.Skyline;

namespace SkylinePrism.App;

/// <summary>
/// The Differential visualization pane: load the corrected matrix from the current output directory
/// and offer three views over a two-group contrast - a limma moderated-t Volcano, a sample PCA, and a
/// peptide Detection-frequency test. All statistics live in SkylinePrism.Core; this is presentation.
/// </summary>
public partial class MainWindow
{
    private static readonly string[] DiffPalette =
    {
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    };

    private static readonly HashSet<string> ReservedMetaColumns =
        new(StringComparer.Ordinal) { "sample", "sample_type", "batch" };

    /// <summary>
    /// Marker sizes for the scatter views (Volcano and the two Detection plots), the significant
    /// series a little larger so it still reads as the emphasized one.
    /// </summary>
    /// <remarks>
    /// Raised from 6/7. At that size a point on a full-screen Volcano was a two-pixel speck: hard to
    /// see at all on a sparse plot, and much smaller than the 18 px
    /// <see cref="QcPlotChrome.HoverRadiusPx"/> that decides what a click lands on - so the target
    /// was far bigger than the thing it was aiming at, which reads as a plot that does not respond.
    /// Named rather than repeated at each call site, because the three views must agree for the
    /// hover radius to mean the same thing on each.
    /// </remarks>
    private const float DiffPointSize = 9;

    /// <summary><see cref="DiffPointSize"/> for the significant series.</summary>
    private const float DiffSigPointSize = 11;

    private DifferentialDataset? _diffDataset;
    private Dictionary<string, string> _diffLabelById = new(StringComparer.Ordinal);

    /// <summary>Feature id -> gene symbol(s), for enrichment only. No entry where none is known.</summary>
    private Dictionary<string, string> _diffGeneById = new(StringComparer.Ordinal);
    private DetectionMatrixData? _detectionData;
    private string? _detectionDir;
    private string? _clinicalCsvPath;
    private List<QcGroupValue> _diffCovariateValues = new();

    // The two contrast arms, as tick lists. An arm is a SET of metadata values whose samples are
    // pooled, not a single value - which is what lets two control classes be contrasted against the
    // rest as one arm.
    private List<QcGroupValue> _diffAValues = new();
    private List<QcGroupValue> _diffBValues = new();
    private HttpJsonPoster? _diffPoster;
    private bool _diffSuppress;
    private int _diffRequest;

    /// <summary>
    /// Generation counter for the VIEW computations, separate from the loader's.
    /// </summary>
    /// <remarks>
    /// The loader has always had one; the views did not, and every selector in this pane now
    /// re-runs on change - design, test, prior, correction, level, view - so two contrasts are
    /// routinely in flight at once. Without a token the SLOWER one wins whenever it happens to
    /// finish last, painting a plot and a hit table that the controls no longer describe. Checked
    /// after every await, because that is where another run can have started.
    /// </remarks>
    private int _diffViewRequest;
    private bool _diffLoaded;
    private string? _diffLoadedDir;
    private FeatureLevel _diffLoadedLevel;

    // Volcano click-to-boxplot state: each plotted point's location + feature id, the per-feature row for
    // the title, and the two groups' sample indices so the boxplot can pull that feature's abundances.
    private List<(ScottPlot.Coordinates Loc, string FeatureId)> _volcanoPoints = new();
    private Dictionary<string, DifferentialRow> _volcanoRowById = new(StringComparer.Ordinal);
    private List<int> _volcanoGroupA = new();
    private List<int> _volcanoGroupB = new();
    private string _volcanoAName = "A";
    private string _volcanoBName = "B";
    private FeatureDetailWindow? _featureDetailWindow;

    // Hover readout and selection ring, both created hidden by RenderVolcano and then only moved.
    // Moving a plottable and refreshing is far cheaper than re-rendering, and it means the highlight
    // does not depend on keeping the DifferentialResult alive to redraw from.
    private ScottPlot.Plottables.Marker? _volcanoHoverMarker;
    private ScottPlot.Plottables.Text? _volcanoHoverText;
    private ScottPlot.Plottables.Marker? _volcanoSelMarker;

    /// <summary>The feature the plot ring and the grid row are both pointing at, or null.</summary>
    private string? _volcanoSelectedId;

    /// <summary>
    /// Set while one of the two selections is being driven from the other. The grid raises
    /// SelectionChanged when its SelectedItem is set in code, so without this a plot click would
    /// select the row, which would re-enter and select the point, and each hop would re-issue the
    /// Skyline selection and reopen the detail window.
    /// </summary>
    private bool _volcanoSyncing;

    private enum DiffView
    {
        Volcano,
        Detection,
        Enrichment,
    }

    private sealed record VolcanoRow(string Feature, double Log2FC, double P, double AdjP, string FeatureId);

    private sealed record DetRow(string Peptide, double RateA, double RateB, double P, double Q);

    private sealed record DetGlmRow(string Peptide, double RateA, double RateB, double LogOR, double P, double Q);

    private sealed record EnrichRow(string Source, string Term, double PValue, double Fold);

    private HttpJsonPoster DiffPoster => _diffPoster ??= new HttpJsonPoster();

    /// <summary>Forget the loaded matrix so the pane reloads on its next show (new dir / new run).</summary>
    private void InvalidateDifferential()
    {
        _diffLoaded = false;
        _diffLoadedDir = null;
        _diffDataset = null;
        _detectionData = null;
        _detectionDir = null;
    }

    private FeatureLevel DiffSelectedLevel() =>
        (DiffLevelCombo.SelectedItem as ComboBoxItem)?.Content as string == "Peptide"
            ? FeatureLevel.Peptide
            : FeatureLevel.Protein;

    /// <summary>
    /// The variance prior the Prior combo is pointing at, defaulting to the lab's choice before the
    /// combo has been populated (the first render happens during window construction).
    /// </summary>
    private VariancePrior DiffSelectedPrior() =>
        ((DiffPriorCombo.SelectedItem as ComboBoxItem)?.Tag as string) switch
        {
            "Global" => VariancePrior.Global,
            "LimmaTrend" => VariancePrior.LimmaTrend,
            "PeptideCount" => VariancePrior.PeptideCount,
            _ => VariancePrior.IntensityTrend,
        };

    /// <summary>The design the Design combo is pointing at.</summary>
    private DifferentialDesign DiffSelectedDesign() =>
        ((DiffDesignCombo.SelectedItem as ComboBoxItem)?.Tag as string) == "Paired"
            ? DifferentialDesign.Paired
            : DifferentialDesign.Unpaired;

    /// <summary>
    /// The pairing key per sample, aligned to the dataset's sample order, or null when no pairing
    /// column is chosen.
    /// </summary>
    private string?[]? DiffSubjectLabels()
    {
        if (_diffDataset is null || DiffPairByCombo.SelectedItem is not string col)
            return null;
        return _diffDataset.MetadataColumns.Contains(col) ? _diffDataset.MetadataValues(col) : null;
    }

    /// <summary>The estimator the Test combo is pointing at.</summary>
    private DifferentialTest DiffSelectedTest() =>
        ((DiffTestCombo.SelectedItem as ComboBoxItem)?.Tag as string) switch
        {
            "WelchT" => DifferentialTest.WelchT,
            "StudentT" => DifferentialTest.StudentT,
            "MannWhitney" => DifferentialTest.MannWhitney,
            "PairedT" => DifferentialTest.PairedT,
            "Wilcoxon" => DifferentialTest.Wilcoxon,
            _ => DifferentialTest.ModeratedT,
        };

    /// <summary>The multiple-testing correction the Correct combo is pointing at.</summary>
    private MultipleTesting DiffSelectedCorrection() =>
        ((DiffCorrectionCombo.SelectedItem as ComboBoxItem)?.Tag as string) switch
        {
            "BenjaminiYekutieli" => MultipleTesting.BenjaminiYekutieli,
            "Holm" => MultipleTesting.Holm,
            "Bonferroni" => MultipleTesting.Bonferroni,
            "None" => MultipleTesting.None,
            _ => MultipleTesting.BenjaminiHochberg,
        };

    /// <summary>
    /// The run's QC and reference sample columns, for fitting the variance prior on controls rather
    /// than on the contrast groups. Null when there are too few to fit anything.
    /// </summary>
    /// <remarks>
    /// Two is the minimum a variance can be computed from at all, so fewer than that is not a thin
    /// prior but no prior - the caller turns the option off rather than falling back silently.
    /// </remarks>
    private IReadOnlyList<IReadOnlyList<int>>? DiffControlColumns()
    {
        if (_diffDataset is null || !_diffDataset.MetadataColumns.Contains("sample_type"))
            return null;

        var types = _diffDataset.MetadataValues("sample_type");

        // One group per control TYPE, not one pooled set: the variance is computed within a group,
        // so pooling QC and reference - different materials, injected at different amounts - would
        // count the systematic gap between them as measurement noise.
        var groups = Enumerable.Range(0, types.Length)
            .Where(i => QcGroupFilter.IsControlValue(types[i]))
            .GroupBy(i => types[i], StringComparer.Ordinal)
            .Select(g => (IReadOnlyList<int>)g.ToList())
            .Where(g => g.Count >= 2) // a group of one has no variance to contribute
            .ToList();

        return groups.Count > 0 ? groups : null;
    }

    /// <summary>What the contrast views should run: the current selections, as Core sees them.</summary>
    private DifferentialOptions DiffOptions(IReadOnlyList<Covariate>? covariates) =>
        new()
        {
            Covariates = covariates,
            Design = DiffSelectedDesign(),
            Test = DiffSelectedTest(),
            Prior = DiffSelectedPrior(),
            Correction = DiffSelectedCorrection(),
            SubjectLabels = DiffSubjectLabels(),
            PeptideCounts = _diffDataset?.PeptideCounts,
            PriorGroupColumns = DiffPriorFromControlsCheck.IsChecked == true
                ? DiffControlColumns()
                : null,
            MinPerGroup = 2,
        };

    /// <summary>
    /// Show only the controls the selected test actually uses.
    /// </summary>
    /// <remarks>
    /// The same two rules the QC pane's <c>UpdateQcControls</c> documents, for the same reasons.
    /// The variance prior is <b>hidden</b> outside the moderated t, because it means nothing there
    /// and this row is already crowded - a disabled control still invites a click. "Adjust for" is
    /// <b>greyed</b> rather than hidden, because it is a real and common setting that simply cannot
    /// be honored by a test with no design matrix; hiding it would make a ticked covariate vanish
    /// with the control, and leaving it live would imply the contrast had been adjusted when it had
    /// not.
    /// </remarks>
    private void UpdateDiffControls()
    {
        var paired = DiffSelectedDesign() == DifferentialDesign.Paired;

        // The pairing column only means something under the paired design.
        var pairVisibility = paired ? Visibility.Visible : Visibility.Collapsed;
        DiffPairByLabel.Visibility = pairVisibility;
        DiffPairByCombo.Visibility = pairVisibility;

        // A test belongs to one design or the other. Collapsed AND disabled, because WPF's
        // arrow-key and type-ahead selection skip only disabled items, so hiding alone would leave
        // an inapplicable test one keypress away.
        ShowTest(DiffTestWelchItem, !paired);
        ShowTest(DiffTestStudentItem, !paired);
        ShowTest(DiffTestMannWhitneyItem, !paired);
        ShowTest(DiffTestPairedTItem, paired);
        ShowTest(DiffTestWilcoxonItem, paired);

        // Changing the design can strand the selection on a test that no longer applies. Fall back
        // to the moderated t, which is valid under both, rather than leaving an invisible selection.
        if (DiffTestCombo.SelectedItem is ComboBoxItem { IsEnabled: false })
        {
            _diffSuppress = true;
            DiffTestCombo.SelectedIndex = 0;
            _diffSuppress = false;
        }

        var moderated = DiffSelectedTest() == DifferentialTest.ModeratedT;
        var priorVisibility = moderated ? Visibility.Visible : Visibility.Collapsed;
        DiffPriorLabel.Visibility = priorVisibility;
        DiffPriorCombo.Visibility = priorVisibility;

        // The count-based priors need n_peptides, which only the protein matrix carries.
        var hasCounts = _diffDataset?.PeptideCounts is not null;
        ShowTest(DiffPriorPeptideCountItem, hasCounts);
        if (DiffPriorCombo.SelectedItem is ComboBoxItem { IsEnabled: false })
        {
            _diffSuppress = true;
            DiffPriorCombo.SelectedIndex = 0; // Intensity trend, valid at either level
            _diffSuppress = false;
        }

        // Fitting the prior on controls only means something for a prior that HAS a per-feature
        // scale, and only where the run actually has control replicates to fit it on.
        var controls = DiffControlColumns();
        var scaled = DiffSelectedPrior() != VariancePrior.Global;
        DiffPriorFromControlsCheck.Visibility = moderated && scaled
            ? Visibility.Visible
            : Visibility.Collapsed;
        DiffPriorFromControlsCheck.IsEnabled = controls is not null;
        if (controls is null && DiffPriorFromControlsCheck.IsChecked == true)
            DiffPriorFromControlsCheck.IsChecked = false;

        DiffCovariatesLabel.IsEnabled = moderated;
        DiffCovariatesCombo.IsEnabled = moderated;
    }

    private static void ShowTest(ComboBoxItem item, bool applies)
    {
        item.Visibility = applies ? Visibility.Visible : Visibility.Collapsed;
        item.IsEnabled = applies;
    }

    private DiffView DiffSelectedView() =>
        ((DiffViewCombo.SelectedItem as ComboBoxItem)?.Content as string) switch
        {
            "Detection" => DiffView.Detection,
            "Enrichment" => DiffView.Enrichment,
            _ => DiffView.Volcano,
        };

    /// <summary>Load the corrected matrix + metadata for the current output directory and level.</summary>
    private async Task LoadDifferentialAsync()
    {
        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir) || !Directory.Exists(dir))
        {
            InvalidateDifferential();
            ClearDiffOutput();
            DiffStatusText.Text = "Set a PRISM output directory above to run a contrast.";
            return;
        }

        _diffSuppress = true;
        try
        {
            if (DiffViewCombo.SelectedItem is null)
                DiffViewCombo.SelectedIndex = 0;
            if (DiffLevelCombo.SelectedItem is null)
                DiffLevelCombo.SelectedIndex = 0;
            // Index 0 is Intensity trend - the lab's default, and deliberately NOT set in the XAML.
            if (DiffPriorCombo.SelectedItem is null)
                DiffPriorCombo.SelectedIndex = 0;
            if (DiffTestCombo.SelectedItem is null)
                DiffTestCombo.SelectedIndex = 0; // Moderated t
            if (DiffDesignCombo.SelectedItem is null)
                DiffDesignCombo.SelectedIndex = 0; // Unpaired
            if (DiffCorrectionCombo.SelectedItem is null)
                DiffCorrectionCombo.SelectedIndex = 0; // Benjamini-Hochberg
        }
        finally
        {
            _diffSuppress = false;
        }

        var level = DiffSelectedLevel();
        if (_diffLoaded && _diffLoadedDir == dir && _diffLoadedLevel == level)
            return; // already current for this directory + level

        ClearDiffOutput();
        DiffStatusText.Text = "Loading...";
        var request = ++_diffRequest;

        DifferentialDataset ds;
        try
        {
            ds = await Task.Run(() => DifferentialDataset.Load(dir, level));
        }
        catch (Exception ex)
        {
            if (request == _diffRequest)
            {
                InvalidateDifferential();
                DiffStatusText.Text = "Load failed: " + ex.Message;
            }

            return;
        }

        if (request != _diffRequest)
            return; // a newer load superseded this one

        _diffDataset = ds;
        _diffLoaded = true;
        _diffLoadedDir = dir;
        _diffLoadedLevel = level;
        _diffLabelById = new Dictionary<string, string>(StringComparer.Ordinal);
        // Genes are kept in their own map: unlike the label there is NO falling back to the feature
        // id, because a protein group id or a peptide sequence is not a gene symbol and enrichment
        // would submit it as one.
        _diffGeneById = new Dictionary<string, string>(StringComparer.Ordinal);
        for (var i = 0; i < ds.FeatureIds.Length; i++)
        {
            _diffLabelById[ds.FeatureIds[i]] =
                string.IsNullOrEmpty(ds.FeatureLabels[i]) ? ds.FeatureIds[i] : ds.FeatureLabels[i];
            if (!string.IsNullOrEmpty(ds.FeatureGenes[i]))
                _diffGeneById[ds.FeatureIds[i]] = ds.FeatureGenes[i];
        }

        // Re-apply a previously attached clinical CSV to the freshly loaded dataset (best-effort).
        if (_clinicalCsvPath is not null && File.Exists(_clinicalCsvPath))
        {
            try
            {
                ds.AttachClinical(_clinicalCsvPath);
            }
            catch
            {
                // a mismatched clinical file just adds nothing; not fatal
            }
        }

        _diffSuppress = true;
        try
        {
            DiffGroupByCombo.ItemsSource = ds.MetadataColumns;
            DiffGroupByCombo.SelectedItem = DefaultContrastColumn(ds);
        }
        finally
        {
            _diffSuppress = false;
        }

        PopulateDiffGroupValues();
        UpdateDiffControls();
        DiffStatusText.Text =
            $"Loaded {ds.FeatureIds.Length} {level.ToString().ToLowerInvariant()} features x " +
            $"{ds.SampleIds.Length} samples. Pick groups and Run.";
    }

    /// <summary>Prefer the first non-reserved metadata column with at least two values; else sample_type.</summary>
    private static string? DefaultContrastColumn(DifferentialDataset ds)
    {
        foreach (var col in ds.MetadataColumns)
        {
            if (ReservedMetaColumns.Contains(col))
                continue;
            var distinct = ds.MetadataValues(col).Where(v => !string.IsNullOrEmpty(v)).Distinct().Count();
            if (distinct >= 2)
                return col;
        }

        return ds.MetadataColumns.FirstOrDefault(c => c == "sample_type")
            ?? ds.MetadataColumns.FirstOrDefault();
    }

    private void PopulateDiffGroupValues()
    {
        if (_diffDataset is null || DiffGroupByCombo.SelectedItem is not string col)
            return;

        var values = _diffDataset.MetadataValues(col)
            .Where(v => !string.IsNullOrEmpty(v))
            .Distinct()
            .OrderBy(v => v, StringComparer.Ordinal)
            .ToList();

        // Default to the first two values, which is what the single pickers did; anything more is
        // an explicit choice by the user.
        _diffAValues = values
            .Select((v, i) => new QcGroupValue
            {
                Name = v!, IsSelected = i == 0, Changed = UpdateDiffArmSummaries,
            })
            .ToList();
        _diffBValues = values
            .Select((v, i) => new QcGroupValue
            {
                Name = v!, IsSelected = i == 1, Changed = UpdateDiffArmSummaries,
            })
            .ToList();

        DiffACombo.ItemsSource = _diffAValues;
        DiffBCombo.ItemsSource = _diffBValues;
        UpdateDiffArmSummaries();

        // Any metadata column can identify a subject except the one being contrasted - pairing by
        // the contrast column itself would put every subject in one arm.
        var pairCandidates = _diffDataset.MetadataColumns
            .Where(m => !string.Equals(m, col, StringComparison.Ordinal))
            .ToList();
        var keepPair = DiffPairByCombo.SelectedItem as string;
        _diffSuppress = true;
        try
        {
            DiffPairByCombo.ItemsSource = pairCandidates;
            DiffPairByCombo.SelectedItem = keepPair is not null && pairCandidates.Contains(keepPair)
                ? keepPair
                : null;
        }
        finally
        {
            _diffSuppress = false;
        }

        PopulateDiffCovariates(col);
    }

    /// <summary>
    /// The closed-state text of the two arm pickers, and of the covariates picker beside them.
    /// </summary>
    /// <remarks>
    /// A tick-list ComboBox has no SelectedItem, so WPF has nothing to display when it is closed and
    /// the text stays at whatever the XAML set. Every other tick list in this window writes its own
    /// summary; the covariates one did not, so it read "(none)" however many covariates were ticked -
    /// fixed here rather than left as the odd one out.
    /// </remarks>
    private void UpdateDiffArmSummaries()
    {
        DiffACombo.Text = SummarizeArm(_diffAValues);
        DiffBCombo.Text = SummarizeArm(_diffBValues);
        var covariates = _diffCovariateValues.Where(v => v.IsSelected).Select(v => v.Name).ToList();
        DiffCovariatesCombo.Text = covariates.Count == 0 ? "(none)" : string.Join(", ", covariates);
    }

    /// <summary>
    /// " + " rather than ", ": the values are POOLED into one arm, and a comma reads like a list of
    /// separate things to compare.
    /// </summary>
    private static string SummarizeArm(IReadOnlyList<QcGroupValue> values)
    {
        var on = values.Where(v => v.IsSelected).Select(v => v.Name).ToList();
        return on.Count == 0 ? "(pick one or more)" : string.Join(" + ", on);
    }

    private void PopulateDiffCovariates(string groupByColumn)
    {
        if (_diffDataset is null)
            return;

        // Any metadata column can be a covariate except the sample id itself and the contrast column.
        _diffCovariateValues = _diffDataset.MetadataColumns
            .Where(c => c != groupByColumn && c != "sample")
            .Select(c => new QcGroupValue { Name = c, Changed = UpdateDiffArmSummaries })
            .ToList();
        DiffCovariatesCombo.ItemsSource = _diffCovariateValues;
    }

    /// <summary>Ticked covariates, with values aligned to <paramref name="targetSampleIds"/>, or null if none.</summary>
    private IReadOnlyList<Covariate>? SelectedCovariatesFor(IReadOnlyList<string> targetSampleIds)
    {
        if (_diffDataset is null)
            return null;
        var chosen = _diffCovariateValues.Where(v => v.IsSelected).Select(v => v.Name).ToList();
        if (chosen.Count == 0)
            return null;

        var indexById = new Dictionary<string, int>(StringComparer.Ordinal);
        for (var i = 0; i < _diffDataset.SampleIds.Length; i++)
            indexById[_diffDataset.SampleIds[i]] = i;

        var result = new List<Covariate>(chosen.Count);
        foreach (var col in chosen)
        {
            var colValues = _diffDataset.MetadataValues(col);
            var aligned = new string?[targetSampleIds.Count];
            for (var k = 0; k < targetSampleIds.Count; k++)
                aligned[k] = indexById.TryGetValue(targetSampleIds[k], out var di) ? colValues[di] : null;
            result.Add(Covariate.FromMetadata(col, aligned));
        }

        return result;
    }

    private async void OnDiffLevelChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_diffSuppress || !IsInitialized)
            return;
        try
        {
            await LoadDifferentialAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffLevelChanged), ex);
        }
    }

    private void OnDiffGroupByChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_diffSuppress)
            return;
        PopulateDiffGroupValues();
    }

    /// <summary>
    /// Choose an external clinical metadata CSV (a top-level input, beside the metadata report). The
    /// path is remembered and applied whenever a differential dataset is loaded; if one is already
    /// loaded, it is joined immediately and the Differential pane's selectors refresh.
    /// </summary>
    private void OnBrowseClinical(object sender, RoutedEventArgs e)
    {
        var dlg = new Microsoft.Win32.OpenFileDialog
        {
            Title = "Choose clinical metadata CSV",
            Filter = "CSV files (*.csv)|*.csv|All files (*.*)|*.*",
            CheckFileExists = true,
        };
        if (dlg.ShowDialog(this) != true)
            return;

        _clinicalCsvPath = dlg.FileName;
        ClinicalCsvBox.Text = dlg.FileName;

        // If a differential dataset is already loaded, apply the join now; otherwise it will be applied
        // the next time one loads (LoadDifferentialAsync re-applies _clinicalCsvPath).
        if (_diffDataset is not null)
            ApplyClinicalToLoadedDataset();
    }

    /// <summary>
    /// Fold the clinical columns just joined to the differential dataset into the QC pane's
    /// grouping sources, so every plot in the window can be grouped by them.
    /// </summary>
    /// <remarks>
    /// Values are taken from the dataset rather than re-read from the CSV on purpose: the join
    /// picked a key column by best match rate, and re-deriving it here could pick a different one
    /// and colour the plot by a slightly different assignment than the one the user is reading the
    /// Volcano against.
    /// </remarks>
    private void PublishClinicalToQcPane(IReadOnlyList<string> addedColumns)
    {
        if (_diffDataset is null || addedColumns.Count == 0)
            return;

        var values = new Dictionary<string, string?[]>(StringComparer.Ordinal);
        foreach (var c in addedColumns)
            values[c] = _diffDataset.MetadataValues(c);

        var byId = new Dictionary<string, int>(StringComparer.Ordinal);
        for (var i = 0; i < _diffDataset.SampleIds.Length; i++)
            byId[_diffDataset.SampleIds[i]] = i;

        var clinical = SampleAnnotationTable.FromValues(
            _diffDataset.SampleIds, addedColumns,
            (id, col) => byId.TryGetValue(id, out var idx) ? values[col][idx] : null);

        // MergedWith keeps what is already there on a name clash, so a clinical column sharing a
        // name with one Skyline exported does not quietly replace it.
        _qcExtraAnnotations = _qcExtraAnnotations.MergedWith(clinical);
        PopulateGroupCombos();
    }

    /// <summary>Join the remembered clinical CSV to the loaded dataset and refresh the selectors.</summary>
    private void ApplyClinicalToLoadedDataset()
    {
        if (_diffDataset is null || _clinicalCsvPath is null || !File.Exists(_clinicalCsvPath))
            return;

        try
        {
            var result = _diffDataset.AttachClinical(_clinicalCsvPath);
            if (result.KeyColumn is null || result.AddedColumns.Count == 0)
            {
                DiffStatusText.Text =
                    $"Clinical CSV: no column matched the samples (best match rate {result.MatchRate:P0}). " +
                    "Nothing was added.";
                return;
            }

            // The QC pane draws the only sample PCA now, so its Group-by list has to learn about
            // these columns too - otherwise attaching a clinical CSV would enrich the Volcano's
            // covariates and leave the PCA unable to colour by any of them, which is the split
            // that put a second PCA in this pane in the first place.
            PublishClinicalToQcPane(result.AddedColumns);

            // The Markers pane keeps its own cache for the OTHER feature level, and that copy was
            // loaded before this join existed. GetMarkersDatasetAsync reuses the Differential
            // dataset when the directory and level both match - so the stale case is specifically
            // the two panes sitting on different levels, where the marker cache is returned as-is
            // and the clinical columns promised here never appear in it.
            InvalidateMarkers();

            // Refresh the group-by choices so the new clinical columns appear; select the first one.
            _diffSuppress = true;
            try
            {
                DiffGroupByCombo.ItemsSource = null;
                DiffGroupByCombo.ItemsSource = _diffDataset.MetadataColumns;
                DiffGroupByCombo.SelectedItem = result.AddedColumns[0];
            }
            finally
            {
                _diffSuppress = false;
            }

            PopulateDiffGroupValues();
            DiffStatusText.Text =
                $"Attached {result.AddedColumns.Count} clinical column(s) via key '{result.KeyColumn}' " +
                $"(matched {result.MatchRate:P0} of samples): {string.Join(", ", result.AddedColumns)}.";
        }
        catch (Exception ex)
        {
            DiffStatusText.Text = $"Could not attach clinical CSV: {ex.Message}";
        }
    }

    private async void OnDiffViewChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_diffSuppress || _diffDataset is null)
            return;
        try
        {
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffViewChanged), ex);
        }
    }

    private async void OnRunDifferential(object sender, RoutedEventArgs e)
    {
        try
        {
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnRunDifferential), ex);
        }
    }

    private async Task RunCurrentViewAsync()
    {
        if (_diffDataset is null)
        {
            DiffStatusText.Text = "Load a run first (point the output directory at a finished PRISM run).";
            return;
        }

        // Anything already in flight is now stale, whatever order it finishes in.
        var request = ++_diffViewRequest;
        switch (DiffSelectedView())
        {
            case DiffView.Detection:
                await RunDetectionAsync(request);
                break;
            case DiffView.Enrichment:
                await RunEnrichmentAsync(request);
                break;
            default:
                await RunVolcanoAsync(request);
                break;
        }
    }

    /// <summary>Whether this run is still the current one. False means drop everything and paint nothing.</summary>
    private bool StillCurrent(int request) => request == _diffViewRequest;

    private bool TryGetGroups(out string col, out List<int> groupA, out List<int> groupB,
        out string aVal, out string bVal)
    {
        col = string.Empty;
        aVal = string.Empty;
        bVal = string.Empty;
        groupA = new List<int>();
        groupB = new List<int>();
        if (_diffDataset is null || DiffGroupByCombo.SelectedItem is not string c)
            return false;

        var aSet = _diffAValues.Where(v => v.IsSelected).Select(v => v.Name).ToList();
        var bSet = _diffBValues.Where(v => v.IsSelected).Select(v => v.Name).ToList();
        if (aSet.Count == 0 || bSet.Count == 0)
            return false;

        // A value in both arms would put the same samples on both sides of the contrast, which is
        // not a contrast. Refused here rather than silently dropped from one side, because which
        // side it was dropped from would change the answer.
        if (aSet.Intersect(bSet, StringComparer.Ordinal).Any())
            return false;

        col = c;
        aVal = string.Join(" + ", aSet);
        bVal = string.Join(" + ", bSet);
        var inA = new HashSet<string>(aSet, StringComparer.Ordinal);
        var inB = new HashSet<string>(bSet, StringComparer.Ordinal);
        var meta = _diffDataset.MetadataValues(c);
        for (var j = 0; j < meta.Length; j++)
        {
            if (meta[j] is not { } v)
                continue;
            if (inA.Contains(v))
                groupA.Add(j);
            else if (inB.Contains(v))
                groupB.Add(j);
        }

        return groupA.Count > 0 && groupB.Count > 0;
    }

    private async Task RunVolcanoAsync(int request)
    {
        if (!TryGetGroups(out _, out var a, out var b, out var aVal, out var bVal))
        {
            DiffStatusText.Text = "Pick a group-by column, then tick at least one value for each "
                + "arm. A value cannot be in both.";
            return;
        }

        var dataset = _diffDataset!;
        var covariates = SelectedCovariatesFor(dataset.SampleIds);
        var options = DiffOptions(covariates);
        DifferentialResult res;
        try
        {
            res = await Task.Run(() =>
                Differential.Run(dataset.ExprLog2, dataset.FeatureIds, a, b, options));
        }
        catch (Exception ex) when (ex is ArgumentException or InvalidOperationException)
        {
            if (StillCurrent(request))
                DiffStatusText.Text = "Cannot run this contrast: " + ex.Message;
            return;
        }

        if (!StillCurrent(request))
            return;

        _volcanoGroupA = a;
        _volcanoGroupB = b;
        _volcanoAName = aVal;
        _volcanoBName = bVal;
        RenderVolcano(res);
        DiffGrid.ItemsSource = res.Rows.Select(r => new VolcanoRow(
            _diffLabelById.GetValueOrDefault(r.FeatureId, r.FeatureId), r.LogFc, r.PValue, r.AdjPValue,
            r.FeatureId))
            .ToList();

        var nSig = res.Rows.Count(r => r.AdjPValue < 0.05 && Math.Abs(r.LogFc) >= 1.0);
        var adj = res.CovariatesUsed.Count > 0 ? $"; adjusted for {string.Join(", ", res.CovariatesUsed)}" : string.Empty;
        // Name the method. With a menu this size the status line is the only record of what
        // produced a hit list, and any message Core raised (an unhonoured covariate, a prior that
        // could not be fitted) belongs beside it rather than nowhere.
        var note = res.Messages.Count > 0 ? " " + string.Join(" ", res.Messages) : string.Empty;
        // res.NA/NB, not a.Count/b.Count: a paired design drops unmatched subjects, so the arms the
        // test actually used can be smaller than the arms that were picked. Reporting the picked
        // sizes would credit the result with samples that took no part in it.
        DiffStatusText.Text =
            $"{options.Describe()}: {aVal} (n={res.NA}) vs {bVal} (n={res.NB}) - "
            + $"{res.NFeaturesTested} tested, {nSig} significant{adj}.{note} "
            + "Click a point (or a row) for its boxplot; it also selects in Skyline. Hover for the gene.";
    }

    private async Task RunDetectionAsync(int request)
    {
        if (!TryGetGroups(out _, out var a, out var b, out var aVal, out var bVal))
        {
            DiffStatusText.Text = "Pick a group-by column, then tick at least one value for each "
                + "arm. A value cannot be in both.";
            return;
        }

        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir))
        {
            DiffStatusText.Text = "Set a PRISM output directory above.";
            return;
        }

        DetectionMatrixData det;
        if (_detectionData is not null && _detectionDir == dir)
        {
            det = _detectionData;
        }
        else
        {
            if (_isRunning)
            {
                DiffStatusText.Text = "A PRISM run is in progress - wait for it to finish before reading detection.";
                return;
            }

            DiffStatusText.Text = "Loading detection matrix from merged_data...";
            try
            {
                det = await Task.Run(() => DetectionMatrix.Load(dir!, 0.01, null));
                if (!StillCurrent(request))
                    return;
            }
            catch (Exception ex)
            {
                DiffStatusText.Text = "Detection load failed (needs merged_data in the output dir): " + ex.Message;
                return;
            }

            _detectionData = det;
            _detectionDir = dir;
        }

        var dataset = _diffDataset!;
        var detIndex = det.SampleIds.Select((s, i) => (s, i)).ToDictionary(x => x.s, x => x.i, StringComparer.Ordinal);
        var aCols = a.Select(j => dataset.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
        var bCols = b.Select(j => dataset.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
        if (aCols.Count == 0 || bCols.Count == 0)
        {
            DiffStatusText.Text = "The selected samples were not found in the detection matrix.";
            return;
        }

        var dropped = a.Count - aCols.Count + (b.Count - bCols.Count);
        var droppedNote = dropped > 0 ? $" ({dropped} samples not in merged_data)" : string.Empty;

        // A ticked covariate switches to the Firth-penalized GLM (adjusted detection); otherwise Fisher.
        var covariates = SelectedCovariatesFor(det.SampleIds);
        if (covariates is not null)
        {
            DetectionGlmResult glm;
            try
            {
                if (!StillCurrent(request))
                    return;
                glm = await Task.Run(() =>
                    DetectionGlm.Run(det.Matrix, det.PeptideIds, aCols, bCols, covariates));
            }
            catch (Exception ex) when (ex is ArgumentException or InvalidOperationException)
            {
                DiffStatusText.Text = "Adjusted detection failed: " + ex.Message;
                return;
            }

            if (!glm.Identifiable)
            {
                ClearDiffOutput();
                DiffStatusText.Text = "Adjusted detection is not identifiable (group confounded with the "
                    + $"covariates, R^2={glm.GroupCollinearityR2:0.00}). Use the unadjusted view.";
                return;
            }

            RenderDetectionGlm(glm.Rows);
            DiffGrid.ItemsSource = glm.Rows.Take(1000)
                .Select(r => new DetGlmRow(r.PeptideId, r.RateA, r.RateB, r.LogOr, r.P, r.Q)).ToList();
            DiffStatusText.Text =
                $"Adjusted detection (Firth GLM): {aVal} (n={aCols.Count}) vs {bVal} (n={bCols.Count}), "
                + $"adjusted for {string.Join(", ", glm.CovariatesUsed)}, {glm.Rows.Count} peptides{droppedNote}.";
            return;
        }

        IReadOnlyList<DetectionRow> rows;
        try
        {
            rows = await Task.Run(() => DetectionTest.Run(det.Matrix, det.PeptideIds, aCols, bCols));
            if (!StillCurrent(request))
                return;
        }
        catch (Exception ex) when (ex is ArgumentException or InvalidOperationException)
        {
            DiffStatusText.Text = "Detection test failed: " + ex.Message;
            return;
        }

        RenderDetection(rows);
        DiffGrid.ItemsSource = rows.Take(1000)
            .Select(r => new DetRow(r.PeptideId, r.RateA, r.RateB, r.P, r.Q)).ToList();
        DiffStatusText.Text =
            $"Detection (peptide-level, DetectionQValue < 0.01): {aVal} (n={aCols.Count}) vs " +
            $"{bVal} (n={bCols.Count}) over {rows.Count} peptides{droppedNote}.";
    }

    private async Task RunEnrichmentAsync(int request)
    {
        if (!TryGetGroups(out _, out var a, out var b, out var aVal, out var bVal))
        {
            DiffStatusText.Text = "Pick a group-by column, then tick at least one value for each "
                + "arm. A value cannot be in both.";
            return;
        }

        var dataset = _diffDataset!;

        // Enrichment is about GENES, so it reads FeatureGenes and never the display label. At
        // peptide level the label is the modified sequence, and a sequence passes CleanSymbols
        // untouched - it only drops empty/"nan" and splits delimiters - so using the label asked
        // g:Profiler about a list of peptide sequences and captioned the empty answer as a gene
        // enrichment. leading_gene_name is stamped onto corrected_peptides for exactly this, so
        // peptide-level enrichment works; it is only impossible on a peptide file written before
        // that column existed, and then we say so rather than guess.
        if (dataset.FeatureGenes.All(string.IsNullOrEmpty))
        {
            ClearDiffOutput();
            DiffStatusText.Text =
                "Enrichment needs gene symbols and this run carries none - its "
                + (dataset.Level == FeatureLevel.Protein ? "protein" : "peptide")
                + " matrix has no leading_gene_name column. Re-run the pipeline to add it.";
            return;
        }

        var covariates = SelectedCovariatesFor(dataset.SampleIds);
        var options = DiffOptions(covariates);
        DifferentialResult res;
        try
        {
            res = await Task.Run(() =>
                Differential.Run(dataset.ExprLog2, dataset.FeatureIds, a, b, options));
        }
        catch (Exception ex) when (ex is ArgumentException or InvalidOperationException)
        {
            if (StillCurrent(request))
                DiffStatusText.Text = "Cannot run this contrast: " + ex.Message;
            return;
        }

        if (!StillCurrent(request))
            return;

        var (sig, background) = Enrichment.SigAndBackgroundGenes(
            res, fid => _diffGeneById.GetValueOrDefault(fid), 0.05, 1.0);
        if (sig.Count == 0)
        {
            ClearDiffOutput();
            DiffStatusText.Text =
                $"No significant genes (adj.P < 0.05, |log2FC| >= 1) for {aVal} vs {bVal} - nothing to enrich.";
            return;
        }

        DiffStatusText.Text = $"Querying g:Profiler for {sig.Count} genes...";
        List<EnrichmentTerm> terms;
        try
        {
            terms = await Task.Run(() => Enrichment.GProfiler(sig, background, DiffPoster));
            if (!StillCurrent(request))
                return;
        }
        catch (Exception ex)
        {
            DiffStatusText.Text = "Enrichment request failed (needs internet access): " + ex.Message;
            return;
        }

        RenderEnrichment(terms);
        DiffGrid.ItemsSource = terms.Take(1000)
            .Select(t => new EnrichRow(t.Source, t.TermName, t.PValue, t.FoldEnrichment)).ToList();
        DiffStatusText.Text = terms.Count == 0
            ? $"No enriched terms for {sig.Count} significant genes (background {background.Count})."
            : $"{terms.Count} enriched terms for {sig.Count} significant genes (background {background.Count}).";
    }

    private void RenderEnrichment(IReadOnlyList<EnrichmentTerm> terms)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;
        if (terms.Count > 0)
        {
            var ys = terms.Take(15).Select(t => -Math.Log10(Math.Max(t.PValue, 1e-300))).ToArray();
            plt.Add.Bars(ys);
            plt.XLabel("top enriched terms (ranked by p)");
            plt.YLabel("-log10 p (g:SCS)");
        }

        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }

    private void RenderDetectionGlm(IReadOnlyList<DetectionGlmRow> rows)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;

        var bgX = new List<double>();
        var bgY = new List<double>();
        var sigX = new List<double>();
        var sigY = new List<double>();
        foreach (var r in rows)
        {
            if (!double.IsFinite(r.LogOr) || !double.IsFinite(r.P))
                continue;
            var y = -Math.Log10(Math.Max(r.P, 1e-300));
            if (r.Q < 0.05)
            {
                sigX.Add(r.LogOr);
                sigY.Add(y);
            }
            else
            {
                bgX.Add(r.LogOr);
                bgY.Add(y);
            }
        }

        AddMarkers(plt, bgX, bgY, "#b8c4d0", DiffPointSize, "q >= 0.05");
        AddMarkers(plt, sigX, sigY, "#2ca02c", DiffSigPointSize, "q < 0.05");
        plt.Add.VerticalLine(0.0);
        plt.ShowLegend();
        plt.XLabel("log odds ratio (B / A)");
        plt.YLabel("-log10 P");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }

    private void OnDiffGridAutoGeneratingColumn(object sender, DataGridAutoGeneratingColumnEventArgs e)
    {
        // FeatureId backs click-to-boxplot from the grid; it is not a column the user needs to see.
        if (e.PropertyName == "FeatureId")
        {
            e.Cancel = true;
            return;
        }

        var (header, format) = e.PropertyName switch
        {
            "Log2FC" => ("log2FC", "0.###"),
            "AdjP" => ("adj.P", "0.##e0"),
            "P" => ("P", "0.##e0"),
            "PValue" => ("p", "0.##e0"),
            "Q" => ("q", "0.##e0"),
            "LogOR" => ("logOR", "0.###"),
            "RateA" => ("rate A", "0.00"),
            "RateB" => ("rate B", "0.00"),
            "VariancePct" => ("variance %", "0.0"),
            "Fold" => ("fold", "0.0"),
            "Source" => ("source", null),
            "Term" => ("term", null),
            _ => (e.PropertyName, (string?)null),
        };

        e.Column.Header = header;
        if (format is not null && e.Column is DataGridTextColumn text && text.Binding is Binding binding)
            binding.StringFormat = format;
    }

    private void ClearDiffOutput()
    {
        DiffGrid.ItemsSource = null;
        DiffPlot.Reset();
        DiffPlot.Refresh();
    }

    /// <summary>
    /// On the Volcano view, a click near a point makes that feature the selection: it rings the
    /// point, selects its row in the hit table, opens the per-feature boxplot, and selects the
    /// protein or peptide in Skyline.
    /// </summary>
    private void OnDiffPlotMouseDown(object sender, System.Windows.Input.MouseButtonEventArgs e)
    {
        try
        {
            if (DiffSelectedView() != DiffView.Volcano || _volcanoPoints.Count == 0 || _diffDataset is null)
                return;

            var plt = DiffPlot.Plot;
            var pos = e.GetPosition(DiffPlot);
            var scale = DiffPlot.DisplayScale;
            var cursor = new ScottPlot.Pixel(pos.X * scale, pos.Y * scale);
            var idx = QcPlotChrome.NearestPoint(
                _volcanoPoints.Select(p => plt.GetPixel(p.Loc)).ToList(), cursor);
            if (idx < 0)
                return;

            SelectVolcanoFeature(_volcanoPoints[idx].FeatureId);
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffPlotMouseDown), ex);
        }
    }

    /// <summary>
    /// Changing the design changes which tests apply and whether a pairing column is needed, so the
    /// controls are updated before anything is re-run.
    /// </summary>
    private async void OnDiffDesignChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            if (!IsInitialized)
                return;
            UpdateDiffControls();
            if (_diffSuppress || _diffDataset is null)
                return;
            // A paired design with no pairing column yet is an incomplete request, not an error:
            // say what is missing and wait rather than throwing from the run.
            if (DiffSelectedDesign() == DifferentialDesign.Paired && DiffSubjectLabels() is null)
            {
                DiffStatusText.Text = "Paired: choose the metadata column that identifies the subject "
                    + "under 'Pair by'.";
                return;
            }

            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffDesignChanged), ex);
        }
    }

    /// <summary>Choosing the pairing column completes a paired request, so run it.</summary>
    private async void OnDiffPairByChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            if (!IsInitialized || _diffSuppress || _diffDataset is null)
                return;
            if (DiffSelectedDesign() != DifferentialDesign.Paired)
                return; // the column is only consulted by the paired design
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffPairByChanged), ex);
        }
    }

    /// <summary>
    /// Changing the estimator re-runs the current view, and shows or hides the controls that only
    /// some estimators use.
    /// </summary>
    private async void OnDiffTestChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            if (!IsInitialized)
                return;
            // Outside the suppress guard: the controls must follow the combo even while a load is
            // populating the pane, or they would be left describing the previous test.
            UpdateDiffControls();
            if (_diffSuppress || _diffDataset is null)
                return;
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffTestChanged), ex);
        }
    }

    /// <summary>
    /// Fitting the prior on controls rather than on the contrast groups changes every moderated
    /// p-value, so it re-runs.
    /// </summary>
    private async void OnDiffPriorSourceChanged(object sender, RoutedEventArgs e)
    {
        try
        {
            if (!IsInitialized || _diffSuppress || _diffDataset is null)
                return;
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffPriorSourceChanged), ex);
        }
    }

    /// <summary>Changing the correction re-runs: it changes every adjusted p in the table.</summary>
    private async void OnDiffCorrectionChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            if (!IsInitialized || _diffSuppress || _diffDataset is null)
                return;
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffCorrectionChanged), ex);
        }
    }

    /// <summary>
    /// Changing the variance prior re-runs the current view, the way changing the View already does -
    /// it is a different estimator, so the plot on screen is no longer the one the controls describe.
    /// </summary>
    private async void OnDiffPriorChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            // IsInitialized as well as the suppress flag: this combo has a SelectionChanged handler,
            // and if it ever gains a XAML-side default it would fire part way through
            // InitializeComponent. See XamlInitializationOrderTests.
            if (!IsInitialized)
                return;
            UpdateDiffControls(); // the control-prior checkbox does not apply to the global prior
            if (_diffSuppress || _diffDataset is null)
                return;
            await RunCurrentViewAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffPriorChanged), ex);
        }
    }

    /// <summary>
    /// Selecting a hit row does everything a volcano click does - the same selection, reached from
    /// the other side - so the ring moves to that feature's point and Skyline follows too.
    /// </summary>
    private void OnDiffGridSelectionChanged(object sender, SelectionChangedEventArgs e)
    {
        try
        {
            if (DiffGrid.SelectedItem is VolcanoRow vr)
                SelectVolcanoFeature(vr.FeatureId);
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffGridSelectionChanged), ex);
        }
    }

    /// <summary>
    /// Select the feature's protein or peptide in the running Skyline document, and say what
    /// happened in the status line.
    /// </summary>
    /// <remarks>
    /// <para>Shares the Dynamic Range pane's locator machinery rather than repeating it: the same
    /// cached document tree, the same "try each of PRISM's groups then fall back to the sequence's
    /// first occurrence" precedence, and the same distinction between "could not read the tree" and
    /// "this one is not in it" - which look identical to a user but need opposite responses.</para>
    /// <para>Nothing happens at all when PRISM is running standalone; that is reported once, in the
    /// status line, rather than being silent.</para>
    /// </remarks>
    private void SelectFeatureInSkyline(string featureId)
    {
        var identity = _diffDataset?.IdentityOf(featureId);
        if (identity is null)
            return;

        var label = identity.Describe();
        if (_session is null)
        {
            // Said rather than skipped: the click asked for something specific and nothing at all
            // happened, which is indistinguishable from a click that missed the point.
            DiffStatusText.Text = $"{label} - not attached to a running Skyline, so nothing to select.";
            return;
        }

        // The DATASET's level, not _diffLoadedLevel: the identity just came out of this dataset, and
        // resolving it against the other level's document tree would look up a peptide among
        // proteins.
        var level = _diffDataset!.Level == FeatureLevel.Protein
            ? AbundanceLevel.Protein
            : AbundanceLevel.Peptide;

        // ResolveLocator speaks AbundanceEntry, so the identity is presented as one. Only the
        // identity fields matter here - the abundance/rank fields are never read on this path.
        var entry = new AbundanceEntry(
            Key: identity.FeatureId,
            Label: identity.Label,
            Accession: identity.Accessions.FirstOrDefault(),
            Gene: identity.Genes.FirstOrDefault(),
            ProteinName: identity.ProteinNames.FirstOrDefault(),
            MeanAbundance: 0, Log10Abundance: 0, Rank: 0, SamplesUsed: 0)
        {
            ProteinGroups = identity.ProteinGroups,
            ProteinNames = identity.ProteinNames,
        };

        var locator = ResolveLocatorLocked(entry, level, out var viaFallback, out var treeUnavailable);
        if (locator is null)
        {
            DiffStatusText.Text = treeUnavailable
                ? $"{label} - could not read the document tree from Skyline (it may be busy; see the "
                  + "Log tab). Click again to retry."
                : $"{label} - no matching element in the Skyline document (PRISM's protein grouping "
                  + "can differ from the document's).";
            return;
        }

        var driver = new SkylineReportDriver(_session, Log);
        if (!driver.SelectElement(locator))
        {
            DiffStatusText.Text = $"Could not select {label} in Skyline.";
            return;
        }

        DiffStatusText.Text = $"Selected {label} in Skyline"
            + (viaFallback
                ? " - under the first protein in the document tree, since PRISM's grouping did not "
                  + "match a protein node."
                : ".");
    }

    private void ShowFeatureDetail(string featureId)
    {
        if (_diffDataset is null)
            return;
        var row = Array.IndexOf(_diffDataset.FeatureIds, featureId);
        if (row < 0)
            return;

        // Values and replicate names are collected together so they stay index-aligned: the
        // non-finite cells are skipped, and a name list built separately would silently shift.
        var aVals = new List<double>();
        var aReps = new List<string>();
        foreach (var s in _volcanoGroupA)
        {
            var v = _diffDataset.ExprLog2[row, s];
            if (!double.IsFinite(v))
                continue;
            aVals.Add(v);
            aReps.Add(_diffDataset.SampleIds[s]);
        }

        var bVals = new List<double>();
        var bReps = new List<string>();
        foreach (var s in _volcanoGroupB)
        {
            var v = _diffDataset.ExprLog2[row, s];
            if (!double.IsFinite(v))
                continue;
            bVals.Add(v);
            bReps.Add(_diffDataset.SampleIds[s]);
        }

        var label = _diffLabelById.GetValueOrDefault(featureId, featureId);
        _volcanoRowById.TryGetValue(featureId, out var dr);

        if (_featureDetailWindow is null)
        {
            _featureDetailWindow = new FeatureDetailWindow { Owner = this };
            _featureDetailWindow.Closed += (_, _) => _featureDetailWindow = null;
        }

        _featureDetailWindow.ShowFeature(label, featureId, dr?.LogFc ?? double.NaN,
            dr?.AdjPValue ?? double.NaN,
            _volcanoAName, aVals, aReps, _volcanoBName, bVals, bReps);
        _featureDetailWindow.Show();
        _featureDetailWindow.Activate();
    }

    private void RenderVolcano(DifferentialResult res)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;

        var bgX = new List<double>();
        var bgY = new List<double>();
        var sigX = new List<double>();
        var sigY = new List<double>();
        _volcanoPoints = new List<(ScottPlot.Coordinates, string)>(res.Rows.Count);
        _volcanoRowById = new Dictionary<string, DifferentialRow>(StringComparer.Ordinal);
        foreach (var r in res.Rows)
        {
            // The ADJUSTED p-value, because that is what decides a hit here (AdjPValue < 0.05) and
            // what the axis says. Differential.Run always applies Benjamini-Hochberg, so AdjPValue is
            // always a q-value - there is no raw-only mode to fall back to. Plotting raw p while
            // colouring by q put the cut-off line at whatever raw p the weakest surviving hit
            // happened to have, which moved with the data and matched no number the reader could see.
            var y = -Math.Log10(Math.Max(r.AdjPValue, 1e-300));
            if (double.IsFinite(r.LogFc) && double.IsFinite(y))
                _volcanoPoints.Add((new ScottPlot.Coordinates(r.LogFc, y), r.FeatureId));
            _volcanoRowById[r.FeatureId] = r;
            if (r.AdjPValue < 0.05 && Math.Abs(r.LogFc) >= 1.0)
            {
                sigX.Add(r.LogFc);
                sigY.Add(y);
            }
            else
            {
                bgX.Add(r.LogFc);
                bgY.Add(y);
            }
        }

        AddMarkers(plt, bgX, bgY, "#b8c4d0", DiffPointSize, "not significant");
        AddMarkers(plt, sigX, sigY, "#d62728", DiffSigPointSize, "significant");
        plt.Add.VerticalLine(1.0);
        plt.Add.VerticalLine(-1.0);
        // Now an exact, readable threshold rather than a data-dependent one: q = 0.05.
        plt.Add.HorizontalLine(-Math.Log10(0.05));

        // Both overlays belong to this Plot instance, so they are recreated with it and must be
        // re-seeded rather than carried over from the previous contrast.
        AddVolcanoOverlays(plt);

        plt.ShowLegend();
        plt.XLabel("log2 fold change (B / A)");
        // Benjamini-Hochberg is unconditional, so name what is actually on the axis. Ties in the
        // adjusted values are expected and show up as horizontal bands - that is a property of BH,
        // not a rendering fault.
        plt.YLabel("-log10(adjusted p-value)");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }

    /// <summary>
    /// The Volcano's hidden hover readout (ring + text) and its selection ring, seeded anywhere -
    /// they are positioned when something is hovered or selected.
    /// </summary>
    private void AddVolcanoOverlays(ScottPlot.Plot plt)
    {
        var hover = plt.Add.Marker(
            0, 0, ScottPlot.MarkerShape.OpenCircle, DiffPointSize + 11, ScottPlot.Colors.Black);
        hover.IsVisible = false;
        _volcanoHoverMarker = hover;

        var text = plt.Add.Text(" ", 0, 0);
        // 20 to match the QC pane's readout: StyleQcPlot draws tick labels at 24, so a smaller size
        // makes the one piece of text a user leans in to read the smallest on the plot.
        PlotRenderer.StyleTextLabel(text, 20, bold: true);
        text.LabelFontColor = ScottPlot.Colors.Black;
        text.LabelBackgroundColor = ScottPlot.Colors.White.WithAlpha(0.85);
        text.LabelAlignment = ScottPlot.Alignment.LowerLeft;
        text.IsVisible = false;
        _volcanoHoverText = text;

        // Distinguishable from the hover ring at a glance: bigger, thicker, and in the significant
        // colour rather than black, because the two can be on screen at the same time.
        var sel = plt.Add.Marker(
            0, 0, ScottPlot.MarkerShape.OpenCircle, DiffSigPointSize + 15,
            ScottPlot.Color.FromHex("#d62728"));
        sel.MarkerLineWidth = 3;
        sel.IsVisible = false;
        _volcanoSelMarker = sel;

        // A re-render is a new contrast, so nothing is selected until the user picks again - and the
        // old id would point into a hit list that no longer contains it.
        _volcanoSelectedId = null;
    }

    /// <summary>
    /// Show what the cursor is over: the feature's label plus the gene and protein it belongs to,
    /// which for a peptide is the only place that context appears on this plot.
    /// </summary>
    private void OnDiffPlotMouseMove(object sender, System.Windows.Input.MouseEventArgs e)
    {
        try
        {
            if (_volcanoHoverMarker is null || _volcanoHoverText is null)
                return;
            if (DiffSelectedView() != DiffView.Volcano || _volcanoPoints.Count == 0 || _diffDataset is null)
                return;

            var plt = DiffPlot.Plot;
            var pos = e.GetPosition(DiffPlot);
            var scale = DiffPlot.DisplayScale;
            var idx = QcPlotChrome.NearestPoint(
                _volcanoPoints.Select(pt => plt.GetPixel(pt.Loc)).ToList(),
                new ScottPlot.Pixel(pos.X * scale, pos.Y * scale));

            if (idx >= 0)
            {
                var point = _volcanoPoints[idx];
                _volcanoHoverMarker.Location = point.Loc;
                _volcanoHoverMarker.IsVisible = true;
                _volcanoHoverText.Location = point.Loc;
                _volcanoHoverText.LabelText = DescribeFeature(point.FeatureId);
                _volcanoHoverText.IsVisible = true;
                DiffPlot.Refresh();
            }
            else if (_volcanoHoverMarker.IsVisible)
            {
                _volcanoHoverMarker.IsVisible = false;
                _volcanoHoverText.IsVisible = false;
                DiffPlot.Refresh();
            }
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnDiffPlotMouseMove), ex);
        }
    }

    /// <summary>
    /// One line for a feature: its label, and the gene / protein / sharing behind it. Falls back to
    /// the id when the matrix carried no identity columns, rather than inventing any.
    /// </summary>
    private string DescribeFeature(string featureId)
    {
        var identity = _diffDataset?.IdentityOf(featureId);
        if (identity is not null)
            return identity.Describe();
        return _diffLabelById.GetValueOrDefault(featureId, featureId);
    }

    /// <summary>
    /// Make <paramref name="featureId"/> the selected feature everywhere: the ring on the plot, the
    /// row in the hit table, the per-feature boxplot, and the selection in Skyline.
    /// </summary>
    /// <remarks>
    /// One entry point for both directions deliberately. The grid and the plot each raise an event
    /// when the other drives it, so routing both through here with <see cref="_volcanoSyncing"/> set
    /// is what stops a click echoing back and forth - each echo would otherwise re-issue the Skyline
    /// RPC and reopen the detail window.
    /// </remarks>
    private void SelectVolcanoFeature(string featureId)
    {
        if (_volcanoSyncing)
            return;

        _volcanoSyncing = true;
        try
        {
            _volcanoSelectedId = featureId;
            HighlightVolcanoPoint(featureId);
            HighlightGridRow(featureId);
            ShowFeatureDetail(featureId);
            SelectFeatureInSkyline(featureId);
        }
        finally
        {
            _volcanoSyncing = false;
        }
    }

    /// <summary>Move the selection ring onto a feature's point, or hide it if it has none plotted.</summary>
    private void HighlightVolcanoPoint(string featureId)
    {
        if (_volcanoSelMarker is null)
            return;

        foreach (var (loc, id) in _volcanoPoints)
        {
            if (!string.Equals(id, featureId, StringComparison.Ordinal))
                continue;
            _volcanoSelMarker.Location = loc;
            _volcanoSelMarker.IsVisible = true;
            BringIntoView(loc);
            DiffPlot.Refresh();
            return;
        }

        // A feature with a non-finite fold change or p-value is in the table but not on the plot.
        _volcanoSelMarker.IsVisible = false;
        DiffPlot.Refresh();
    }

    /// <summary>
    /// Pan the Volcano so a selected point is actually on screen, keeping the current zoom.
    /// </summary>
    /// <remarks>
    /// Selecting from the hit table can name a point that is outside the view - the table is ranked
    /// by significance and is not affected by zooming, so after zooming into one corner most rows
    /// refer to points that are not visible. Ringing one of those rings nothing the user can see,
    /// which looks identical to the selection having failed.
    /// <para>A minimal pan rather than a re-fit or a re-center: the zoom the user chose is
    /// deliberate, and re-centering on every off-screen pick makes the plot jump further than it
    /// needs to. The point is brought just inside the edge it was past, with a margin so it does not
    /// sit under the axis.</para>
    /// </remarks>
    private void BringIntoView(ScottPlot.Coordinates loc)
    {
        var plt = DiffPlot.Plot;
        var limits = plt.Axes.GetLimits();
        var xSpan = limits.Right - limits.Left;
        var ySpan = limits.Top - limits.Bottom;
        if (!(xSpan > 0) || !(ySpan > 0))
            return; // an unrendered plot has no meaningful limits to preserve

        const double marginFraction = 0.08;
        var xMargin = xSpan * marginFraction;
        var yMargin = ySpan * marginFraction;

        var dx = 0.0;
        if (loc.X < limits.Left + xMargin)
            dx = loc.X - (limits.Left + xMargin);
        else if (loc.X > limits.Right - xMargin)
            dx = loc.X - (limits.Right - xMargin);

        var dy = 0.0;
        if (loc.Y < limits.Bottom + yMargin)
            dy = loc.Y - (limits.Bottom + yMargin);
        else if (loc.Y > limits.Top - yMargin)
            dy = loc.Y - (limits.Top - yMargin);

        if (dx == 0 && dy == 0)
            return; // already comfortably in view: leave the axes exactly as the user set them

        plt.Axes.SetLimits(
            limits.Left + dx, limits.Right + dx,
            limits.Bottom + dy, limits.Top + dy);
    }

    /// <summary>Select a feature's row in the hit table and scroll it into view.</summary>
    private void HighlightGridRow(string featureId)
    {
        if (DiffGrid.ItemsSource is not IEnumerable<VolcanoRow> rows)
            return;

        var row = rows.FirstOrDefault(r => string.Equals(r.FeatureId, featureId, StringComparison.Ordinal));
        if (row is null)
            return;

        DiffGrid.SelectedItem = row;
        DiffGrid.ScrollIntoView(row);
    }

    private void RenderDetection(IReadOnlyList<DetectionRow> rows)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;

        var bgX = new List<double>();
        var bgY = new List<double>();
        var sigX = new List<double>();
        var sigY = new List<double>();
        foreach (var r in rows)
        {
            var x = r.RateB - r.RateA;
            var y = -Math.Log10(Math.Max(r.P, 1e-300));
            if (r.Q < 0.05)
            {
                sigX.Add(x);
                sigY.Add(y);
            }
            else
            {
                bgX.Add(x);
                bgY.Add(y);
            }
        }

        AddMarkers(plt, bgX, bgY, "#b8c4d0", DiffPointSize, "q >= 0.05");
        AddMarkers(plt, sigX, sigY, "#2ca02c", DiffSigPointSize, "q < 0.05");
        plt.Add.VerticalLine(0.0);
        plt.ShowLegend();
        plt.XLabel("detection rate difference (B - A)");
        plt.YLabel("-log10 P");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }

    private static void AddMarkers(ScottPlot.Plot plt, List<double> xs, List<double> ys, string hex,
        float size, string? legend)
    {
        if (xs.Count == 0)
            return;
        var m = plt.Add.Markers(xs.ToArray(), ys.ToArray());
        m.Color = ScottPlot.Color.FromHex(hex);
        m.MarkerSize = size;
        if (legend is not null)
            m.LegendText = legend;
    }
}
