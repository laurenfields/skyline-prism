using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows;
using ScottPlot;
using SkylinePrism.Core.Numerics;
using SkylinePrism.Core.Visualization;

namespace SkylinePrism.App;

/// <summary>
/// A per-feature detail popup: a boxplot of one feature's log2 abundance split by the two contrast
/// groups, with jittered per-sample points overlaid. Opened by clicking a point on the Volcano view,
/// mirroring the explorer's "Feature detail" panel. Reused across clicks (one window, redrawn).
/// </summary>
public partial class FeatureDetailWindow : Window
{
    private const string AColor = "#1f77b4"; // group A: blue
    private const string BColor = "#d62728"; // group B: red

    /// <summary>
    /// Each drawn point's location and the replicate it came from, for the hover readout. A boxplot
    /// hides which sample is which by construction, so an outlier in the swarm is not identifiable
    /// without this - and identifying one is usually the reason for opening this window.
    /// </summary>
    private readonly List<(Coordinates Loc, string Replicate)> _points = new();

    private ScottPlot.Plottables.Marker? _hoverMarker;
    private ScottPlot.Plottables.Text? _hoverText;

    public FeatureDetailWindow()
    {
        InitializeComponent();
        DetailPlot.MouseMove += OnDetailPlotMouseMove;
    }

    /// <summary>Draw (or redraw) the boxplot for one feature.</summary>
    /// <param name="aReplicates">
    /// Replicate names for <paramref name="aValues"/>, same order and length - the caller filters
    /// non-finite values out of both together, so the two stay aligned.
    /// </param>
    /// <param name="bReplicates">Replicate names for <paramref name="bValues"/>.</param>
    public void ShowFeature(string label, string featureId, double logFc, double q,
        string aName, IReadOnlyList<double> aValues, IReadOnlyList<string> aReplicates,
        string bName, IReadOnlyList<double> bValues, IReadOnlyList<string> bReplicates)
    {
        _points.Clear();
        var inv = CultureInfo.InvariantCulture;
        var fc = double.IsFinite(logFc) ? logFc.ToString("0.###", inv) : "n/a";
        var qs = double.IsFinite(q) ? q.ToString("0.##e0", inv) : "n/a";
        HeaderText.Text = label == featureId ? label : $"{label} ({featureId})";
        SubText.Text =
            // "non-missing", never "detected": a Skyline export integrates imputed peak boundaries
            // for every replicate, so the corrected matrix is dense and a finite value here may
            // carry no real signal at all. Calling these counts detections would contradict the
            // pane's own caveat and overstate what the boxplot shows.
            $"log2FC (B/A) = {fc}   adj.P = {qs}   -   {aName}: n={aValues.Count}, "
            + $"{bName}: n={bValues.Count} (non-missing values)";

        var plt = DetailPlot.Plot;
        plt.Clear();

        var boxes = new List<ScottPlot.Box>();
        AddBox(boxes, 0, aValues, AColor);
        AddBox(boxes, 1, bValues, BColor);
        if (boxes.Count > 0)
            plt.Add.Boxes(boxes);

        // Jittered per-sample points, seeded so the layout is stable across redraws of the same feature.
        var rng = new Random(HashCode.Combine(featureId, aValues.Count, bValues.Count) & 0x7fffffff);
        AddJitter(plt, 0, aValues, aReplicates, AColor, rng);
        AddJitter(plt, 1, bValues, bReplicates, BColor, rng);
        AddHoverOverlay(plt);

        plt.Axes.Bottom.TickGenerator = new ScottPlot.TickGenerators.NumericManual(
            new double[] { 0, 1 }, new[] { aName, bName });
        plt.Axes.SetLimitsX(-0.6, 1.6);
        plt.XLabel("group");
        plt.YLabel("log2 abundance");
        PlotRenderer.StyleQcPlot(plt);
        DetailPlot.Refresh();
    }

    private static void AddBox(List<ScottPlot.Box> boxes, double position, IReadOnlyList<double> values,
        string colorHex)
    {
        if (values.Count == 0)
            return;
        var arr = new double[values.Count];
        for (var i = 0; i < values.Count; i++)
            arr[i] = values[i];

        var q1 = Stats.PercentileLinear(arr, 25);
        var med = Stats.PercentileLinear(arr, 50);
        var q3 = Stats.PercentileLinear(arr, 75);
        var iqr = q3 - q1;
        double dataMin = double.PositiveInfinity, dataMax = double.NegativeInfinity;
        foreach (var v in arr)
        {
            if (v < dataMin) dataMin = v;
            if (v > dataMax) dataMax = v;
        }

        boxes.Add(new ScottPlot.Box
        {
            Position = position,
            Width = 0.6,
            BoxMin = q1,
            BoxMiddle = med,
            BoxMax = q3,
            WhiskerMin = Math.Max(dataMin, q1 - 1.5 * iqr),
            WhiskerMax = Math.Min(dataMax, q3 + 1.5 * iqr),
            FillColor = ScottPlot.Color.FromHex(colorHex).WithAlpha((byte)90),
            LineColor = ScottPlot.Color.FromHex(colorHex),
        });
    }

    private void AddJitter(ScottPlot.Plot plt, double center, IReadOnlyList<double> values,
        IReadOnlyList<string> replicates, string colorHex, Random rng)
    {
        if (values.Count == 0)
            return;
        var xs = new double[values.Count];
        var ys = new double[values.Count];
        for (var i = 0; i < values.Count; i++)
        {
            xs[i] = center - 0.28 + rng.NextDouble() * 0.56;
            ys[i] = values[i];
            // The jittered x is what was DRAWN, so it is what the cursor has to be matched against.
            _points.Add((new Coordinates(xs[i], ys[i]),
                i < replicates.Count ? replicates[i] : string.Empty));
        }

        var m = plt.Add.Markers(xs, ys);
        m.Color = ScottPlot.Color.FromHex(colorHex);
        m.MarkerSize = 8;
    }

    /// <summary>The hidden hover readout: a ring and a label, positioned when a point is hovered.</summary>
    private void AddHoverOverlay(ScottPlot.Plot plt)
    {
        var marker = plt.Add.Marker(0, 0, MarkerShape.OpenCircle, 19, Colors.Black);
        marker.IsVisible = false;
        _hoverMarker = marker;

        var text = plt.Add.Text(" ", 0, 0);
        PlotRenderer.StyleTextLabel(text, 18, bold: true);
        text.LabelFontColor = Colors.Black;
        text.LabelBackgroundColor = Colors.White.WithAlpha(0.85);
        text.LabelAlignment = Alignment.LowerLeft;
        text.IsVisible = false;
        _hoverText = text;
    }

    /// <summary>
    /// Name the replicate under the cursor. Uses the same radius as the other panes' readouts so
    /// hovering feels the same everywhere.
    /// </summary>
    private void OnDetailPlotMouseMove(object sender, System.Windows.Input.MouseEventArgs e)
    {
        if (_hoverMarker is null || _hoverText is null || _points.Count == 0)
            return;

        var plt = DetailPlot.Plot;
        var pos = e.GetPosition(DetailPlot);
        var scale = DetailPlot.DisplayScale;
        var cursor = new Pixel(pos.X * scale, pos.Y * scale);

        var best = -1;
        var bestDistance = double.MaxValue;
        for (var i = 0; i < _points.Count; i++)
        {
            var px = plt.GetPixel(_points[i].Loc);
            double dx = px.X - cursor.X, dy = px.Y - cursor.Y;
            var d2 = dx * dx + dy * dy;
            if (d2 < bestDistance)
            {
                bestDistance = d2;
                best = i;
            }
        }

        if (best >= 0 && bestDistance <= QcPlotChrome.HoverRadiusPx * QcPlotChrome.HoverRadiusPx
                      && !string.IsNullOrEmpty(_points[best].Replicate))
        {
            var point = _points[best];
            _hoverMarker.Location = point.Loc;
            _hoverMarker.IsVisible = true;
            _hoverText.Location = point.Loc;
            _hoverText.LabelText = point.Replicate;
            _hoverText.IsVisible = true;
            DetailPlot.Refresh();
        }
        else if (_hoverMarker.IsVisible)
        {
            _hoverMarker.IsVisible = false;
            _hoverText.IsVisible = false;
            DetailPlot.Refresh();
        }
    }
}
