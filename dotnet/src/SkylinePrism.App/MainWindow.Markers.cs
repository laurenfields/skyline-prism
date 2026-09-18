using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Controls;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.Numerics;
using SkylinePrism.Core.Qc;
using SkylinePrism.Core.Visualization;

namespace SkylinePrism.App;

/// <summary>
/// The Markers pane: evaluate a protein list (marker panel) against the corrected matrix. Shows a row
/// z-scored heatmap of the panel's members across groups (or samples) and a per-group boxplot of the
/// panel score. Reuses the Dynamic Range protein lists as the panel options.
/// </summary>
public partial class MainWindow
{
    private DifferentialDataset? _markersDataset;
    private string? _markersDir;
    private FeatureLevel _markersLevel;
    private bool _markersSuppress;

    private sealed record MarkerPanelChoice(string Display, ProteinList List);

    private async Task LoadMarkersAsync()
    {
        var dir = OutputDirBox.Text?.Trim();
        if (string.IsNullOrEmpty(dir) || !Directory.Exists(dir))
        {
            MarkersStatusText.Text = "Set a PRISM output directory above to evaluate a marker panel.";
            return;
        }

        _markersSuppress = true;
        try
        {
            if (MarkersLevelCombo.SelectedItem is null)
                MarkersLevelCombo.SelectedIndex = 0;
            if (MarkersViewCombo.SelectedItem is null)
                MarkersViewCombo.SelectedIndex = 0;
        }
        finally
        {
            _markersSuppress = false;
        }

        var level = MarkersSelectedLevel();
        DifferentialDataset ds;
        try
        {
            ds = await GetMarkersDatasetAsync(dir, level);
        }
        catch (Exception ex)
        {
            MarkersStatusText.Text = "Load failed: " + ex.Message;
            return;
        }

        _markersSuppress = true;
        try
        {
            PopulateMarkersPanels();
            MarkersGroupByCombo.ItemsSource = ds.MetadataColumns;
            if (MarkersGroupByCombo.SelectedItem is null)
                MarkersGroupByCombo.SelectedItem = DefaultContrastColumn(ds);
        }
        finally
        {
            _markersSuppress = false;
        }

        RenderMarkers();
    }

    private async Task<DifferentialDataset> GetMarkersDatasetAsync(string dir, FeatureLevel level)
    {
        // Reuse the Differential pane's dataset when it is the same directory + level (it already has any
        // clinical CSV attached); otherwise load our own and attach the clinical CSV the same way.
        if (_diffDataset is not null && _diffLoadedDir == dir && _diffLoadedLevel == level)
            return _diffDataset;
        if (_markersDataset is not null && _markersDir == dir && _markersLevel == level)
            return _markersDataset;

        var ds = await Task.Run(() => DifferentialDataset.Load(dir, level));
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

        _markersDataset = ds;
        _markersDir = dir;
        _markersLevel = level;
        return ds;
    }

    private void PopulateMarkersPanels()
    {
        var selected = (MarkersPanelCombo.SelectedItem as MarkerPanelChoice)?.List.Name;
        var choices = ProteinListSet.Load().WithBuiltIns()
            .Where(l => l.Members.Count > 0)
            .OrderBy(l => l.Category, StringComparer.Ordinal)
            .ThenBy(l => l.Name, StringComparer.Ordinal)
            .Select(l => new MarkerPanelChoice(
                string.IsNullOrEmpty(l.Category) ? l.Name : $"{l.Category}: {l.Name}", l))
            .ToList();

        MarkersPanelCombo.DisplayMemberPath = nameof(MarkerPanelChoice.Display);
        MarkersPanelCombo.ItemsSource = choices;
        var restore = choices.FirstOrDefault(c => c.List.Name == selected);
        MarkersPanelCombo.SelectedItem = restore ?? choices.FirstOrDefault();
    }

    private FeatureLevel MarkersSelectedLevel() =>
        (MarkersLevelCombo.SelectedItem as ComboBoxItem)?.Content as string == "Peptide"
            ? FeatureLevel.Peptide
            : FeatureLevel.Protein;

    private async void OnMarkersLevelChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_markersSuppress || !IsInitialized)
            return;
        try
        {
            await LoadMarkersAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnMarkersLevelChanged), ex);
        }
    }

    private void OnMarkersPanelChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_markersSuppress)
            return;
        RenderMarkers();
    }

    private void OnMarkersGroupByChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_markersSuppress)
            return;
        RenderMarkers();
    }

    private void OnMarkersViewChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_markersSuppress)
            return;
        RenderMarkers();
    }

    private async void OnMarkersReload(object sender, System.Windows.RoutedEventArgs e)
    {
        try
        {
            _markersDataset = null;
            _markersDir = null;
            await LoadMarkersAsync();
        }
        catch (Exception ex)
        {
            ReportHandlerFailure(nameof(OnMarkersReload), ex);
        }
    }

    private void RenderMarkers()
    {
        var ds = (_diffDataset is not null && _diffLoadedDir == _markersDir
                  && _diffLoadedLevel == _markersLevel)
            ? _diffDataset
            : _markersDataset;
        if (ds is null)
            return;
        if (MarkersPanelCombo.SelectedItem is not MarkerPanelChoice choice
            || MarkersGroupByCombo.SelectedItem is not string groupCol)
            return;

        string?[] groups;
        try
        {
            groups = ds.MetadataValues(groupCol);
        }
        catch (ArgumentException)
        {
            return;
        }

        var perSample = (MarkersViewCombo.SelectedItem as ComboBoxItem)?.Content as string == "Per sample";
        var result = MarkerPanel.Evaluate(ds.ExprLog2, ds.FeatureIds, ds.FeatureLabels, groups,
            ds.SampleIds, choice.List, perSample);

        if (result.MarkerLabels.Length == 0)
        {
            MarkersHeatPlot.Plot.Clear();
            MarkersHeatPlot.Refresh();
            MarkersBoxPlot.Plot.Clear();
            MarkersBoxPlot.Refresh();
            MarkersStatusText.Text =
                $"None of {choice.List.Name}'s {result.Total} members matched a {_markersLevel.ToString().ToLowerInvariant()} feature.";
            MarkersNoteText.Text = result.NotDetected.Count == 0
                ? "No members detected in this run."
                : "Not detected: " + string.Join(", ", result.NotDetected);
            return;
        }

        var heat = MarkersHeatPlot.Plot;
        heat.Clear();
        PlotRenderer.DrawValueHeatmap(heat, result.Heatmap, result.ColumnLabels, result.MarkerLabels,
            result.SymmetricMax, "row z-score", annotate: !perSample && result.MarkerLabels.Length <= 30);
        heat.Title($"{choice.List.Name} (row z-scored log2) - "
            + (perSample ? "per sample" : $"group means by {groupCol}"));
        MarkersHeatPlot.Refresh();

        DrawMarkerBoxplot(result, groupCol);

        MarkersStatusText.Text =
            $"{choice.List.Name}: found {result.Found}/{result.Total} members across "
            + $"{result.GroupNames.Length} groups.";
        MarkersNoteText.Text = result.NotDetected.Count == 0
            ? "All panel members detected. Boxplot below is each sample's mean marker z-score per group."
            : $"Not detected ({result.NotDetected.Count}): " + string.Join(", ", result.NotDetected);
    }

    private void DrawMarkerBoxplot(MarkerPanelResult result, string groupCol)
    {
        var plt = MarkersBoxPlot.Plot;
        plt.Clear();

        var boxes = new List<ScottPlot.Box>();
        for (var g = 0; g < result.GroupNames.Length; g++)
        {
            var vals = result.PanelScoreByGroup[g];
            var hex = DiffPalette[g % DiffPalette.Length];
            if (vals.Length > 0)
            {
                var q1 = Stats.PercentileLinear(vals, 25);
                var med = Stats.PercentileLinear(vals, 50);
                var q3 = Stats.PercentileLinear(vals, 75);
                var iqr = q3 - q1;
                double dataMin = double.PositiveInfinity, dataMax = double.NegativeInfinity;
                foreach (var v in vals)
                {
                    if (v < dataMin) dataMin = v;
                    if (v > dataMax) dataMax = v;
                }

                boxes.Add(new ScottPlot.Box
                {
                    Position = g,
                    Width = 0.6,
                    BoxMin = q1,
                    BoxMiddle = med,
                    BoxMax = q3,
                    WhiskerMin = Math.Max(dataMin, q1 - 1.5 * iqr),
                    WhiskerMax = Math.Min(dataMax, q3 + 1.5 * iqr),
                    FillColor = ScottPlot.Color.FromHex(hex).WithAlpha((byte)90),
                    LineColor = ScottPlot.Color.FromHex(hex),
                });
            }

            // Jittered per-sample points, seeded for a stable layout.
            var rng = new Random(g * 7919 + result.MarkerLabels.Length);
            var xs = new double[vals.Length];
            for (var i = 0; i < vals.Length; i++)
                xs[i] = g - 0.22 + rng.NextDouble() * 0.44;
            if (vals.Length > 0)
            {
                var m = plt.Add.Markers(xs, vals);
                m.Color = ScottPlot.Color.FromHex(hex);
                m.MarkerSize = 5;
            }
        }

        if (boxes.Count > 0)
            plt.Add.Boxes(boxes);

        var pos = new double[result.GroupNames.Length];
        for (var g = 0; g < pos.Length; g++)
            pos[g] = g;
        plt.Axes.Bottom.TickGenerator = new ScottPlot.TickGenerators.NumericManual(pos, result.GroupNames);
        plt.Axes.SetLimitsX(-0.6, result.GroupNames.Length - 0.4);
        plt.XLabel(groupCol);
        plt.YLabel("mean marker z-score");
        PlotRenderer.StyleQcPlot(plt);
        MarkersBoxPlot.Refresh();
    }
}
