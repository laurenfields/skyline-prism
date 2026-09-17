using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.Visualization;

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

    private DifferentialDataset? _diffDataset;
    private Dictionary<string, string> _diffLabelById = new(StringComparer.Ordinal);
    private DetectionMatrixData? _detectionData;
    private string? _detectionDir;
    private bool _diffSuppress;

    private enum DiffView
    {
        Volcano,
        Pca,
        Detection,
    }

    private sealed record VolcanoRow(string Feature, string log2FC, string P, string adj_P);

    private sealed record PcaVarRow(string Component, string Variance);

    private sealed record DetRow(string Peptide, string rate_A, string rate_B, string p, string q);

    private FeatureLevel DiffSelectedLevel() =>
        (DiffLevelCombo.SelectedItem as ComboBoxItem)?.Content as string == "Peptide"
            ? FeatureLevel.Peptide
            : FeatureLevel.Protein;

    private DiffView DiffSelectedView() =>
        ((DiffViewCombo.SelectedItem as ComboBoxItem)?.Content as string) switch
        {
            "PCA" => DiffView.Pca,
            "Detection" => DiffView.Detection,
            _ => DiffView.Volcano,
        };

    /// <summary>Load the corrected matrix + metadata for the current output directory and level.</summary>
    private async Task LoadDifferentialAsync()
    {
        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir) || !Directory.Exists(dir))
        {
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

            var level = DiffSelectedLevel();
            DifferentialDataset ds;
            try
            {
                ds = await Task.Run(() => DifferentialDataset.Load(dir, level));
            }
            catch (Exception ex)
            {
                _diffDataset = null;
                DiffStatusText.Text = "Load failed: " + ex.Message;
                return;
            }

            _diffDataset = ds;
            _detectionData = null; // level/dir changed; the detection matrix is reloaded on demand
            _diffLabelById = new Dictionary<string, string>(StringComparer.Ordinal);
            for (var i = 0; i < ds.FeatureIds.Length; i++)
                _diffLabelById[ds.FeatureIds[i]] =
                    string.IsNullOrEmpty(ds.FeatureLabels[i]) ? ds.FeatureIds[i] : ds.FeatureLabels[i];

            DiffGroupByCombo.ItemsSource = ds.MetadataColumns;
            DiffGroupByCombo.SelectedItem = ds.MetadataColumns.FirstOrDefault(c => c == "sample_type")
                ?? ds.MetadataColumns.FirstOrDefault();
            PopulateDiffGroupValues();

            DiffStatusText.Text =
                $"Loaded {ds.FeatureIds.Length} {level.ToString().ToLowerInvariant()} features x " +
                $"{ds.SampleIds.Length} samples. Pick groups and Run.";
        }
        finally
        {
            _diffSuppress = false;
        }
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

        DiffACombo.ItemsSource = values;
        DiffBCombo.ItemsSource = values;
        if (values.Count >= 2)
        {
            DiffACombo.SelectedIndex = 0;
            DiffBCombo.SelectedIndex = 1;
        }
        else
        {
            DiffACombo.SelectedIndex = -1;
            DiffBCombo.SelectedIndex = -1;
        }
    }

    private async void OnDiffLevelChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_diffSuppress || !IsLoaded)
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
            DiffStatusText.Text = "Run failed: " + ex.Message;
        }
    }

    private async Task RunCurrentViewAsync()
    {
        if (_diffDataset is null)
        {
            DiffStatusText.Text = "Load a run first (point the output directory at a finished PRISM run).";
            return;
        }

        switch (DiffSelectedView())
        {
            case DiffView.Pca:
                RunPca();
                break;
            case DiffView.Detection:
                await RunDetectionAsync();
                break;
            default:
                RunVolcano();
                break;
        }
    }

    private bool TryGetGroups(out string col, out List<int> groupA, out List<int> groupB,
        out string aVal, out string bVal)
    {
        col = string.Empty;
        aVal = string.Empty;
        bVal = string.Empty;
        groupA = new List<int>();
        groupB = new List<int>();
        if (_diffDataset is null || DiffGroupByCombo.SelectedItem is not string c
            || DiffACombo.SelectedItem is not string a || DiffBCombo.SelectedItem is not string b || a == b)
            return false;

        col = c;
        aVal = a;
        bVal = b;
        var meta = _diffDataset.MetadataValues(c);
        groupA = Enumerable.Range(0, meta.Length).Where(j => meta[j] == a).ToList();
        groupB = Enumerable.Range(0, meta.Length).Where(j => meta[j] == b).ToList();
        return true;
    }

    private void RunVolcano()
    {
        if (!TryGetGroups(out _, out var a, out var b, out var aVal, out var bVal))
        {
            DiffStatusText.Text = "Pick a group-by column and two different values.";
            return;
        }

        var res = Differential.Run(_diffDataset!.ExprLog2, _diffDataset.FeatureIds, a, b);
        RenderVolcano(res);

        var inv = CultureInfo.InvariantCulture;
        DiffGrid.ItemsSource = res.Rows.Select(r => new VolcanoRow(
            _diffLabelById.GetValueOrDefault(r.FeatureId, r.FeatureId),
            r.LogFc.ToString("0.###", inv),
            r.PValue.ToString("0.##e0", inv),
            r.AdjPValue.ToString("0.##e0", inv))).ToList();

        var nSig = res.Rows.Count(r => r.AdjPValue < 0.05 && Math.Abs(r.LogFc) >= 1.0);
        DiffStatusText.Text =
            $"{aVal} (n={a.Count}) vs {bVal} (n={b.Count}) - {res.NFeaturesTested} tested, {nSig} significant.";
    }

    private void RunPca()
    {
        if (DiffGroupByCombo.SelectedItem is not string col)
        {
            DiffStatusText.Text = "Pick a group-by column to colour the PCA by.";
            return;
        }

        var all = Enumerable.Range(0, _diffDataset!.SampleIds.Length).ToList();
        PcaResult pca;
        try
        {
            pca = DifferentialPca.Compute(_diffDataset.ExprLog2, _diffDataset.SampleIds, all);
        }
        catch (Exception ex)
        {
            DiffStatusText.Text = "PCA failed: " + ex.Message;
            return;
        }

        RenderPca(pca, col);

        var inv = CultureInfo.InvariantCulture;
        DiffGrid.ItemsSource = pca.VarianceRatio
            .Select((v, i) => new PcaVarRow($"PC{i + 1}", (v * 100).ToString("0.0", inv) + "%"))
            .ToList();
        DiffStatusText.Text =
            $"PCA over {all.Count} samples, {pca.NFeaturesUsed} complete features, coloured by {col}.";
    }

    private async Task RunDetectionAsync()
    {
        if (!TryGetGroups(out _, out var a, out var b, out var aVal, out var bVal))
        {
            DiffStatusText.Text = "Pick a group-by column and two different values.";
            return;
        }

        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir))
            return;

        DetectionMatrixData det;
        if (_detectionData is not null && _detectionDir == dir)
        {
            det = _detectionData;
        }
        else
        {
            DiffStatusText.Text = "Loading detection matrix from merged_data...";
            try
            {
                det = await Task.Run(() => DetectionMatrix.Load(dir!, 0.01, null));
            }
            catch (Exception ex)
            {
                DiffStatusText.Text = "Detection load failed (needs merged_data in the output dir): " + ex.Message;
                return;
            }

            _detectionData = det;
            _detectionDir = dir;
        }

        var detIndex = det.SampleIds.Select((s, i) => (s, i)).ToDictionary(x => x.s, x => x.i, StringComparer.Ordinal);
        var aCols = a.Select(j => _diffDataset!.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
        var bCols = b.Select(j => _diffDataset!.SampleIds[j]).Where(detIndex.ContainsKey).Select(s => detIndex[s]).ToList();
        if (aCols.Count == 0 || bCols.Count == 0)
        {
            DiffStatusText.Text = "The selected samples were not found in the detection matrix.";
            return;
        }

        var rows = DetectionTest.Run(det.Matrix, det.PeptideIds, aCols, bCols);
        RenderDetection(rows);

        var inv = CultureInfo.InvariantCulture;
        DiffGrid.ItemsSource = rows.Take(1000).Select(r => new DetRow(
            r.PeptideId,
            r.RateA.ToString("0.00", inv),
            r.RateB.ToString("0.00", inv),
            r.P.ToString("0.##e0", inv),
            r.Q.ToString("0.##e0", inv))).ToList();
        DiffStatusText.Text =
            $"Detection: {aVal} (n={aCols.Count}) vs {bVal} (n={bCols.Count}) over {rows.Count} peptides.";
    }

    private void RenderVolcano(DifferentialResult res)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;

        var bgX = new List<double>();
        var bgY = new List<double>();
        var sigX = new List<double>();
        var sigY = new List<double>();
        var pThresh = double.NaN;
        foreach (var r in res.Rows)
        {
            var y = -Math.Log10(Math.Max(r.PValue, 1e-300));
            if (r.AdjPValue < 0.05 && Math.Abs(r.LogFc) >= 1.0)
            {
                sigX.Add(r.LogFc);
                sigY.Add(y);
                if (double.IsNaN(pThresh) || r.PValue > pThresh)
                    pThresh = r.PValue;
            }
            else
            {
                bgX.Add(r.LogFc);
                bgY.Add(y);
            }
        }

        AddMarkers(plt, bgX, bgY, "#b8c4d0", 6, null);
        AddMarkers(plt, sigX, sigY, "#d62728", 7, "significant");
        plt.Add.VerticalLine(1.0);
        plt.Add.VerticalLine(-1.0);
        if (!double.IsNaN(pThresh))
            plt.Add.HorizontalLine(-Math.Log10(Math.Max(pThresh, 1e-300)));

        plt.XLabel("log2 fold change (B / A)");
        plt.YLabel("-log10 P");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }

    private void RenderPca(PcaResult pca, string colorColumn)
    {
        DiffPlot.Reset();
        var plt = DiffPlot.Plot;

        if (pca.Scores.GetLength(1) < 2)
        {
            DiffPlot.Refresh();
            return;
        }

        var labels = _diffDataset!.MetadataValues(colorColumn);
        var groups = new Dictionary<string, (List<double> X, List<double> Y)>();
        for (var i = 0; i < pca.SampleIds.Length; i++)
        {
            var g = string.IsNullOrEmpty(labels[i]) ? "(none)" : labels[i]!;
            if (!groups.TryGetValue(g, out var lists))
                groups[g] = lists = (new List<double>(), new List<double>());
            lists.X.Add(pca.Scores[i, 0]);
            lists.Y.Add(pca.Scores[i, 1]);
        }

        var ci = 0;
        foreach (var (g, lists) in groups.OrderBy(kv => kv.Key, StringComparer.Ordinal))
            AddMarkers(plt, lists.X, lists.Y, DiffPalette[ci++ % DiffPalette.Length], 11, g);

        plt.ShowLegend();
        var inv = CultureInfo.InvariantCulture;
        plt.XLabel($"PC1 ({(pca.VarianceRatio[0] * 100).ToString("0.0", inv)}%)");
        plt.YLabel($"PC2 ({(pca.VarianceRatio[1] * 100).ToString("0.0", inv)}%)");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
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

        AddMarkers(plt, bgX, bgY, "#b8c4d0", 6, null);
        AddMarkers(plt, sigX, sigY, "#2ca02c", 7, "q < 0.05");
        plt.Add.VerticalLine(0.0);
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
