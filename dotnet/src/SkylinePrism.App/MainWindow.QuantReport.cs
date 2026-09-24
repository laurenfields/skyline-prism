using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.DifferentialAnalysis.Enrichment;

namespace SkylinePrism.App;

/// <summary>
/// The "Quant report" button: run the Differential pane's current contrast across every view and write a
/// self-contained quant_report.html plus the result tables (CSV) to a quant/ folder, the way the pipeline
/// writes its own outputs. Reuses the pane's contrast, options and rule, the ticked Markers-pane panels,
/// and (opt-in) g:Profiler enrichment.
/// </summary>
public partial class MainWindow
{
    private readonly string[] _quantEnrichmentSources = { "GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG" };

    private async void OnGenerateQuantReport(object sender, RoutedEventArgs e)
    {
        try
        {
            await GenerateQuantReportAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnGenerateQuantReport), ex);
        }
    }

    private async Task GenerateQuantReportAsync()
    {
        if (_diffDataset is null)
        {
            DiffStatusText.Text = "Load a run and set up a contrast before writing a quant report.";
            return;
        }

        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir) || !Directory.Exists(dir))
        {
            DiffStatusText.Text = "Set a PRISM output directory above.";
            return;
        }

        if (_isRunning)
        {
            DiffStatusText.Text = "A PRISM run is in progress - wait for it to finish.";
            return;
        }

        // Everything that reads a WPF control is gathered here, on the UI thread; the compute below runs
        // on a worker.
        var ds = _diffDataset;
        var level = DiffSelectedLevel();
        var rule = DiffRule();
        var corrected = DiffSelectedCorrection() != MultipleTesting.None;
        var effectName = DiffEffectName();
        var isTrend = DiffIsTrend();
        var covariates = WithoutTestedTerm(SelectedCovariatesFor(ds.SampleIds));
        var options = DiffOptions(covariates);
        var designName = DiffSelectedDesign().ToString();
        var testName = options.Test.ToString();
        var correctionName = DiffSelectedCorrection().ToString();

        QuantContrast contrast;
        string contrastLabel;
        double[]? trendX = null;
        List<int> aCols = new(), bCols = new();
        if (isTrend)
        {
            if (DiffTrendColumn() is not { } tcol || DiffTrendValues() is not { } xv)
            {
                DiffStatusText.Text = "Pick a numeric trend column first.";
                return;
            }

            trendX = xv;
            contrast = new QuantContrast(null, null, null, tcol);
            contrastLabel = $"trend over {tcol}";
        }
        else
        {
            if (!TryGetGroups(out var col, out aCols, out bCols, out var aVal, out var bVal))
            {
                DiffStatusText.Text = "Pick a group-by column and two values before writing a report.";
                return;
            }

            contrast = new QuantContrast(col, aVal, bVal, null);
            contrastLabel = $"{bVal} vs {aVal} by {col}";
        }

        var markerColumn = MarkersGroupByCombo.SelectedItem as string ?? contrast.GroupBy;
        var markerPanels = _markersPanelItems.Where(i => i.IsSelected).Select(i => i.List).ToList();
        // Enrichment always runs; it degrades to a skipped-with-note section if there is no network or
        // no significant genes, so there is nothing for the user to opt into.
        const bool wantEnrichment = true;
        var poster = DiffPoster;
        var geneById = _diffGeneById;
        var labelById = _diffLabelById;
        var cachedDetection = _detectionData is not null && _detectionDir == dir ? _detectionData : null;

        DiffStatusText.Text = "Generating quant report...";
        QuantReportButton.IsEnabled = false;
        try
        {
            var (html, detMatrix, note) = await Task.Run(() => Build(
                dir!, ds, level, rule, corrected, effectName, isTrend, options, designName, testName,
                correctionName, contrast, contrastLabel, aCols, bCols, trendX, markerColumn, markerPanels,
                wantEnrichment, poster, geneById, labelById, cachedDetection));

            // Cache the detection matrix so a later Detection-pane run on the same folder reuses it.
            if (detMatrix is not null)
            {
                _detectionData = detMatrix;
                _detectionDir = dir;
            }

            DiffStatusText.Text = $"Quant report written to {html}.{note}";
            try
            {
                Process.Start(new ProcessStartInfo(html) { UseShellExecute = true });
            }
            catch
            {
                // Opening the browser is best-effort; the file is written regardless.
            }
        }
        finally
        {
            QuantReportButton.IsEnabled = true;
        }
    }

    /// <summary>Runs on a worker thread: computes every requested view and writes the report.</summary>
    private (string Html, DetectionMatrixData? Detection, string Note) Build(
        string dir, DifferentialDataset ds, FeatureLevel level, SignificanceRule rule, bool corrected,
        string effectName, bool isTrend, DifferentialOptions options, string designName, string testName,
        string correctionName, QuantContrast contrast, string contrastLabel, List<int> aCols,
        List<int> bCols, double[]? trendX, string? markerColumn, IReadOnlyList<Core.Qc.ProteinList> markerPanels,
        bool wantEnrichment, HttpJsonPoster poster, Dictionary<string, string> geneById,
        Dictionary<string, string> labelById, DetectionMatrixData? cachedDetection)
    {
        var notes = new List<string>();

        var res = isTrend
            ? Differential.RunTrend(ds.ExprLog2, ds.FeatureIds,
                Enumerable.Range(0, ds.SampleIds.Length).ToArray(), trendX!, options)
            : Differential.Run(ds.ExprLog2, ds.FeatureIds, aCols, bCols, options);

        // Detection: two-group designs only, and only where merged_data is present.
        DetectionMatrixData? detMatrix = null;
        IReadOnlyList<DetectionRow>? detection = null;
        if (!isTrend)
        {
            try
            {
                detMatrix = cachedDetection ?? DetectionMatrix.Load(dir, 0.01, null);
                var detIndex = detMatrix.SampleIds
                    .Select((s, i) => (s, i)).ToDictionary(x => x.s, x => x.i, StringComparer.Ordinal);
                var da = aCols.Select(j => ds.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
                var db = bCols.Select(j => ds.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
                if (da.Count >= 2 && db.Count >= 2)
                    detection = DetectionTest.Run(detMatrix.Matrix, detMatrix.PeptideIds, da, db);
                else
                    notes.Add(" Detection skipped: the contrast samples were not both present in merged_data.");
            }
            catch (Exception ex)
            {
                notes.Add(" Detection skipped: " + ex.Message);
            }
        }

        // Enrichment: opt-in, needs network.
        IReadOnlyList<EnrichmentTerm>? enrichment = null;
        if (wantEnrichment)
        {
            try
            {
                var (sig, background) = Enrichment.SigAndBackgroundGenes(
                    res, id => geneById.GetValueOrDefault(id), rule);
                if (sig.Count == 0)
                    notes.Add(" Enrichment skipped: no significant genes to submit.");
                else
                    enrichment = Enrichment.GProfiler(sig, background, poster);
            }
            catch (Exception ex)
            {
                notes.Add(" Enrichment skipped (needs internet): " + ex.Message);
            }
        }

        // Markers: the ticked panels, grouped by the Markers-pane column (or the contrast column).
        var markers = new List<MarkerReportSection>();
        if (markerColumn is not null && markerPanels.Count > 0 && ds.MetadataColumns.Contains(markerColumn))
        {
            var groups = ds.MetadataValues(markerColumn);
            var identities = Enumerable.Range(0, ds.FeatureIds.Length).Select(ds.IdentityOf).ToArray();
            foreach (var panel in markerPanels)
            {
                var r = MarkerPanel.Evaluate(ds.ExprLog2, identities, groups, ds.SampleIds, panel, false);
                markers.Add(new MarkerReportSection(panel.Name, markerColumn, r));
            }
        }
        else if (markerPanels.Count == 0)
        {
            notes.Add(" Markers omitted: no panels ticked in the Markers pane.");
        }

        var quant = new QuantConfig(
            level == FeatureLevel.Peptide ? "peptide" : "protein", contrast, designName, testName,
            res.VariancePrior, correctionName, res.CovariatesUsed, rule.Describe(effectName),
            detection is not null, 0.01, enrichment is not null, _quantEnrichmentSources, "both",
            markers.Select(m => m.PanelName).ToList());

        var inputs = new QuantReportInputs
        {
            Differential = res,
            Rule = rule,
            Corrected = corrected,
            Contrast = contrastLabel,
            EffectName = effectName,
            LabelFor = id => labelById.GetValueOrDefault(id, id),
            Detection = detection,
            Enrichment = enrichment,
            Markers = markers,
            Dataset = ds,
            // Raw per-sample values are written for a two-group contrast, where the columns ARE the two
            // arms. A trend has no arms and its fitted sample subset is not recoverable from the result
            // (some samples can be dropped for a non-finite x or a within-subject singleton), so a
            // per-sample matrix here would list samples that were not in the fit - omitted rather than
            // overclaimed. differential.csv (per feature) is still written for a trend.
            ContrastColumns = isTrend ? null : aCols.Concat(bCols).ToList(),
        };

        var html = QuantReport.Write(dir, quant, inputs);
        return (html, detMatrix, string.Concat(notes));
    }
}
