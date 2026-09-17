using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.Visualization;

namespace SkylinePrism.App;

/// <summary>
/// The Differential visualization pane: load the corrected matrix from the current output directory,
/// pick a two-group contrast off a metadata column, run the limma moderated-t test
/// (<see cref="Differential"/>), and show a volcano plot + ranked hit table. Reads its inputs through
/// <see cref="DifferentialDataset"/>; the statistics live entirely in SkylinePrism.Core.
/// </summary>
public partial class MainWindow
{
    private DifferentialDataset? _diffDataset;
    private Dictionary<string, string> _diffLabelById = new(StringComparer.Ordinal);
    private bool _diffSuppress;

    private sealed record DiffRow(string Label, double LogFc, double P, double AdjP)
    {
        public string LogFcText => LogFc.ToString("0.###", CultureInfo.InvariantCulture);
        public string PText => P.ToString("0.##e0", CultureInfo.InvariantCulture);
        public string AdjPText => AdjP.ToString("0.##e0", CultureInfo.InvariantCulture);
    }

    private FeatureLevel DiffSelectedLevel() =>
        (DiffLevelCombo.SelectedItem as ComboBoxItem)?.Content as string == "Peptide"
            ? FeatureLevel.Peptide
            : FeatureLevel.Protein;

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
                $"{ds.SampleIds.Length} samples. Pick two groups and Run.";
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

    private void OnRunDifferential(object sender, RoutedEventArgs e)
    {
        if (_diffDataset is null)
        {
            DiffStatusText.Text = "Load a run first (point the output directory at a finished PRISM run).";
            return;
        }

        if (DiffGroupByCombo.SelectedItem is not string col
            || DiffACombo.SelectedItem is not string aVal
            || DiffBCombo.SelectedItem is not string bVal)
        {
            DiffStatusText.Text = "Pick a group-by column and two values.";
            return;
        }

        if (aVal == bVal)
        {
            DiffStatusText.Text = "Groups A and B must be different.";
            return;
        }

        var meta = _diffDataset.MetadataValues(col);
        var groupA = Enumerable.Range(0, meta.Length).Where(j => meta[j] == aVal).ToList();
        var groupB = Enumerable.Range(0, meta.Length).Where(j => meta[j] == bVal).ToList();

        try
        {
            var res = Differential.Run(_diffDataset.ExprLog2, _diffDataset.FeatureIds, groupA, groupB);
            RenderVolcano(res);
            DiffGrid.ItemsSource = res.Rows
                .Select(r => new DiffRow(
                    _diffLabelById.GetValueOrDefault(r.FeatureId, r.FeatureId),
                    r.LogFc, r.PValue, r.AdjPValue))
                .ToList();
            DiffStatusText.Text =
                $"{aVal} (n={groupA.Count}) vs {bVal} (n={groupB.Count}) - {res.NFeaturesTested} features tested, " +
                $"{res.Rows.Count(r => r.AdjPValue < 0.05 && Math.Abs(r.LogFc) >= 1.0)} significant.";
        }
        catch (Exception ex)
        {
            DiffStatusText.Text = "Run failed: " + ex.Message;
        }
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

        if (bgX.Count > 0)
        {
            var m = plt.Add.Markers(bgX.ToArray(), bgY.ToArray());
            m.Color = ScottPlot.Color.FromHex("#b8c4d0");
            m.MarkerSize = 6;
        }

        if (sigX.Count > 0)
        {
            var m = plt.Add.Markers(sigX.ToArray(), sigY.ToArray());
            m.Color = ScottPlot.Color.FromHex("#d62728");
            m.MarkerSize = 7;
            m.LegendText = "significant";
        }

        plt.Add.VerticalLine(1.0);
        plt.Add.VerticalLine(-1.0);
        if (!double.IsNaN(pThresh))
            plt.Add.HorizontalLine(-Math.Log10(Math.Max(pThresh, 1e-300)));

        plt.XLabel("log2 fold change (B / A)");
        plt.YLabel("-log10 P");
        PlotRenderer.StyleQcPlot(plt);
        DiffPlot.Refresh();
    }
}
