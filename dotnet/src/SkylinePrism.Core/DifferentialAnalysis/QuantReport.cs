using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using SkylinePrism.Core.Config;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.DifferentialAnalysis.Enrichment;
using SkylinePrism.Core.Pipeline;
using SkylinePrism.Core.Qc;
using SkylinePrism.Core.Visualization;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>One evaluated marker panel, for the report's Markers section.</summary>
public sealed record MarkerReportSection(string PanelName, string GroupColumn, MarkerPanelResult Result);

/// <summary>Everything one quantification report renders. Sections left null are omitted.</summary>
public sealed class QuantReportInputs
{
    public required DifferentialResult Differential { get; init; }
    public required SignificanceRule Rule { get; init; }
    public required bool Corrected { get; init; }
    public required string Contrast { get; init; }
    public required string EffectName { get; init; }
    public required Func<string, string> LabelFor { get; init; }

    public IReadOnlyList<DetectionRow>? Detection { get; init; }
    public IReadOnlyList<EnrichmentTerm>? Enrichment { get; init; }
    public IReadOnlyList<MarkerReportSection>? Markers { get; init; }

    /// <summary>
    /// The loaded dataset, so the export can emit the RAW per-sample abundances behind the summaries
    /// (linear, matching PRISM's corrected parquet). Null omits the raw-value files.
    /// </summary>
    public DifferentialDataset? Dataset { get; init; }

    /// <summary>The sample columns of the contrast (group A then B), for the raw differential values.</summary>
    public IReadOnlyList<int>? ContrastColumns { get; init; }
}

/// <summary>
/// Builds the self-contained <c>quant_report.html</c> - the quantification counterpart to
/// <c>qc_report.html</c>. It shares the QC report's Analysis Information header (via
/// <see cref="QcReport.AppendRunInfo"/>, read from the run's <c>parameters.json</c>), records the quant
/// analysis parameters as a table and re-runnable YAML, and renders a section per requested view
/// (Differential, Detection, Enrichment, Markers) with embedded plots. Alongside it, the same result
/// tables are written as CSVs under a <c>quant/</c> folder, the way the pipeline writes its own outputs.
/// </summary>
public static class QuantReport
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>
    /// Render the report and its companion files. <paramref name="outputDir"/> is the finished PRISM run
    /// (its <c>parameters.json</c> supplies the shared header); the report and CSVs are written under
    /// <c>&lt;outputDir&gt;/quant/</c>. Returns the HTML path.
    /// </summary>
    public static string Write(string outputDir, QuantConfig quant, QuantReportInputs inputs,
        int maxHitRows = 25)
    {
        var quantDir = Path.Combine(outputDir, "quant");
        Directory.CreateDirectory(quantDir);

        File.WriteAllText(Path.Combine(quantDir, "quant_parameters.yaml"), quant.ToYaml());
        File.WriteAllText(Path.Combine(quantDir, "quant_parameters.json"), quant.ToJson());
        WriteCompanionCsvs(quantDir, inputs);

        var runInfo = ReadRunInfoSafe(outputDir);
        var pipelineConfig = ReadPipelineConfigSafe(outputDir);
        var html = BuildHtml(quant, inputs, runInfo, pipelineConfig, maxHitRows);

        var htmlPath = Path.Combine(quantDir, "quant_report.html");
        File.WriteAllText(htmlPath, html);
        return htmlPath;
    }

    private static string BuildHtml(QuantConfig quant, QuantReportInputs inputs,
        Provenance.RunInfo? runInfo, PrismConfig pipelineConfig, int maxHitRows)
    {
        var generatedAt = DateTime.Now.ToString("yyyy-MM-dd HH:mm", Inv);
        var rule = inputs.Rule;
        var res = inputs.Differential;

        var sb = new StringBuilder();
        sb.Append("<!DOCTYPE html>\n<html><head><meta charset=\"utf-8\">"
            + "<title>PRISM Quantification Report</title>\n");
        AppendStyle(sb);
        sb.Append("</head><body><div class=\"container\">");
        sb.Append("<h1>PRISM Quantification Report</h1>");

        // The SAME Analysis Information block the QC report shows: how the data was produced.
        QcReport.AppendRunInfo(sb, runInfo, pipelineConfig, generatedAt);

        AppendQuantParameters(sb, quant, res, rule, inputs.EffectName);

        // --- Differential ---
        sb.Append("<div class=\"section-header\">Differential abundance</div>");
        var volcano = PlotRenderer.DifferentialVolcanoPng(res, rule, inputs.Corrected, inputs.EffectName,
            inputs.LabelFor, inputs.Contrast);
        AppendImage(sb, volcano, "Volcano: " + inputs.Contrast);
        var hits = res.Rows.Where(rule.IsSignificant).ToList();
        sb.Append($"<h3>Significant hits ({hits.Count} of {res.NFeaturesTested} tested)</h3>");
        if (hits.Count == 0)
            sb.Append($"<p class=\"note\">No feature met {HtmlEncode(rule.Describe(inputs.EffectName))}.</p>");
        else
        {
            AppendHitTable(sb, hits.Take(maxHitRows), inputs.EffectName, inputs.LabelFor);
            AppendPreviewNote(sb, hits.Count, maxHitRows, "differential.csv (every tested feature)");
        }

        AppendMessages(sb, "Notes", res.Messages);
        AppendMessages(sb, "Warnings", res.Warnings);

        // --- Detection ---
        if (inputs.Detection is { Count: > 0 } det)
        {
            sb.Append("<div class=\"section-header\">Detection frequency</div>");
            sb.Append("<p class=\"note\">On/off detection from the transition-level data, not the dense "
                + "abundance matrix. A cell counts as detected where DetectionQValue is below the "
                + "threshold; the test asks whether the detection rate differs between the two groups.</p>");
            AppendImage(sb, PlotRenderer.DetectionVolcanoPng(det, rule, inputs.Corrected, inputs.Contrast),
                "Detection volcano");
            var detHits = det.Where(r => (rule.UseAdjusted ? r.Q : r.P) < rule.PThreshold)
                .OrderBy(r => rule.UseAdjusted ? r.Q : r.P).ToList();
            sb.Append($"<h3>Detection differences ({detHits.Count})</h3>");
            AppendDetectionTable(sb, detHits.Take(maxHitRows));
            AppendPreviewNote(sb, detHits.Count, maxHitRows, "detection.csv (every peptide)");
        }

        // --- Enrichment ---
        if (inputs.Enrichment is { Count: > 0 } terms)
        {
            sb.Append("<div class=\"section-header\">Functional enrichment (g:Profiler)</div>");
            AppendImage(sb, PlotRenderer.EnrichmentBarsPng(terms), "Top enriched terms");
            sb.Append($"<h3>Enriched terms ({terms.Count})</h3>");
            AppendEnrichmentTable(sb, terms.Take(maxHitRows));
            AppendPreviewNote(sb, terms.Count, maxHitRows, "enrichment_terms.csv (every term, full gene lists)");
        }

        // --- Markers ---
        if (inputs.Markers is { Count: > 0 } markers)
        {
            sb.Append("<div class=\"section-header\">Marker panels</div>");
            foreach (var m in markers)
            {
                sb.Append($"<h3>{HtmlEncode(m.PanelName)}</h3>");
                sb.Append($"<p class=\"note\">Found {m.Result.Found} of {m.Result.Total} members across "
                    + $"{m.Result.GroupNames.Length} groups"
                    + (m.Result.NotDetected.Count > 0
                        ? $"; not detected: {HtmlEncode(string.Join(", ", m.Result.NotDetected.Take(30)))}"
                        : string.Empty) + ".</p>");
                if (m.Result.GroupNames.Length > 20)
                    sb.Append($"<p class=\"note\">Grouped by <code>{HtmlEncode(m.GroupColumn)}</code>, which "
                        + $"has {m.Result.GroupNames.Length} values - too many to label legibly. A coarser "
                        + "grouping column (e.g. condition rather than a per-vial id) reads far better.</p>");
                if (m.Result.MarkerLabels.Length > 0)
                {
                    AppendImage(sb, PlotRenderer.MarkerHeatmapPng(m.Result,
                        $"{m.PanelName} (row z-scored) by {m.GroupColumn}"), $"{m.PanelName} heatmap");
                    AppendImage(sb, PlotRenderer.MarkerBoxplotPng(m.Result, m.GroupColumn),
                        $"{m.PanelName} panel score");
                }
            }
        }

        sb.Append($"<p class=\"footer\">Generated by Skyline-PRISM (C#) at {generatedAt}. "
            + "Result tables are in the <code>quant/</code> folder beside this report.</p>");
        sb.Append("</div></body></html>");
        return sb.ToString();
    }

    private static void AppendQuantParameters(
        StringBuilder sb, QuantConfig quant, DifferentialResult res, SignificanceRule rule, string effectName)
    {
        sb.Append("<div class=\"box\"><h2>Quantification Parameters</h2><table class=\"kv\">");
        Kv(sb, "Contrast", quant.Contrast.Describe());
        Kv(sb, "Level", quant.Level);
        Kv(sb, "Method", quant.Design + ", " + quant.Test + ", " + quant.Prior + " prior");
        Kv(sb, "Multiple testing", quant.Correction);
        Kv(sb, "Hit rule", rule.Describe(effectName));
        Kv(sb, "Groups", res.IsTrend
            ? $"n = {res.NA}" + (res.NSubjects > 0 ? $" ({res.NSubjects} subjects)" : string.Empty)
            : $"A n = {res.NA}, B n = {res.NB}");
        Kv(sb, "Features tested", $"{res.NFeaturesTested} of {res.NFeaturesTotal}");
        if (res.CovariatesUsed.Count > 0)
            Kv(sb, "Adjusted for", string.Join(", ", res.CovariatesUsed));
        sb.Append("</table>");
        sb.Append("<details><summary>Quantification parameters (YAML)</summary><pre>")
          .Append(HtmlEncode(quant.ToYaml())).Append("</pre></details>");
        sb.Append("</div>");
    }

    // -------- companion CSVs --------

    private static void WriteCompanionCsvs(string quantDir, QuantReportInputs inputs)
    {
        using (var w = new StreamWriter(Path.Combine(quantDir, "differential.csv")))
        {
            w.WriteLine("feature_id,label," + Csv(inputs.EffectName) + ",fc,ave_expr,statistic,p_value,adj_p_value,mean_a,mean_b");
            foreach (var r in inputs.Differential.Rows)
                w.WriteLine(string.Join(",", Csv(r.FeatureId), Csv(inputs.LabelFor(r.FeatureId)),
                    N(r.LogFc), N(r.Fc), N(r.AveExpr), N(r.T), N(r.PValue), N(r.AdjPValue), N(r.MeanA), N(r.MeanB)));
        }

        if (inputs.Detection is { Count: > 0 } det)
            using (var w = new StreamWriter(Path.Combine(quantDir, "detection.csv")))
            {
                w.WriteLine("peptide,detected_a,n_a,detected_b,n_b,rate_a,rate_b,p_value,adj_p_value");
                foreach (var r in det)
                    w.WriteLine(string.Join(",", Csv(r.PeptideId), r.DetA, r.NA, r.DetB, r.NB,
                        N(r.RateA), N(r.RateB), N(r.P), N(r.Q)));
            }

        if (inputs.Enrichment is { Count: > 0 } terms)
            using (var w = new StreamWriter(Path.Combine(quantDir, "enrichment_terms.csv")))
            {
                w.WriteLine("source,term_id,term_name,p_value,fold_enrichment,intersection_size,genes");
                foreach (var t in terms)
                    w.WriteLine(string.Join(",", Csv(t.Source), Csv(t.TermId), Csv(t.TermName),
                        N(t.PValue), N(t.FoldEnrichment), t.IntersectionSize,
                        Csv(string.Join(";", t.IntersectingGenes))));
            }

        if (inputs.Markers is { Count: > 0 } markers)
            foreach (var m in markers)
            {
                var safe = string.Concat(m.PanelName.Select(c => char.IsLetterOrDigit(c) ? c : '_'));
                using var w = new StreamWriter(Path.Combine(quantDir, $"markers_{safe}_zscores.csv"));
                w.WriteLine("marker," + string.Join(",", m.Result.ColumnLabels.Select(Csv)));
                for (var i = 0; i < m.Result.MarkerLabels.Length; i++)
                {
                    var cells = new List<string> { Csv(m.Result.MarkerLabels[i]) };
                    for (var j = 0; j < m.Result.ColumnLabels.Length; j++)
                        cells.Add(N(m.Result.Heatmap[i, j]));
                    w.WriteLine(string.Join(",", cells));
                }
            }

        // Raw per-sample abundances (LINEAR, matching corrected_*.parquet), so the export stands on its
        // own for reanalysis rather than carrying only fold changes and z-scores.
        if (inputs.Dataset is { } ds)
        {
            var rowOf = new Dictionary<string, int>(StringComparer.Ordinal);
            for (var i = 0; i < ds.FeatureIds.Length; i++)
                rowOf[ds.FeatureIds[i]] = i;

            if (inputs.ContrastColumns is { Count: > 0 } cols)
                WriteValuesMatrix(Path.Combine(quantDir, "differential_values.csv"),
                    ds, rowOf, inputs.Differential.Rows.Select(r => (r.FeatureId, inputs.LabelFor(r.FeatureId))), cols);

            if (inputs.Markers is { Count: > 0 } ms)
                foreach (var m in ms)
                {
                    var safe = string.Concat(m.PanelName.Select(c => char.IsLetterOrDigit(c) ? c : '_'));
                    var groups = ds.MetadataValues(m.GroupColumn);
                    var included = Enumerable.Range(0, groups.Length)
                        .Where(s => !string.IsNullOrEmpty(groups[s])).ToList();
                    var rows = m.Result.MarkerFeatureIds
                        .Select((id, i) => (id, m.Result.MarkerLabels[i]));
                    WriteValuesMatrix(Path.Combine(quantDir, $"markers_{safe}_values.csv"),
                        ds, rowOf, rows, included);
                }
        }
    }

    /// <summary>
    /// Write a feature x sample matrix of LINEAR abundances (2^log2), one row per (featureId, label)
    /// over <paramref name="columns"/> of the matrix. Features absent from the matrix are skipped.
    /// </summary>
    private static void WriteValuesMatrix(string path, DifferentialDataset ds,
        IReadOnlyDictionary<string, int> rowOf, IEnumerable<(string Id, string Label)> features,
        IReadOnlyList<int> columns)
    {
        using var w = new StreamWriter(path);
        w.WriteLine("feature_id,label," + string.Join(",", columns.Select(c => Csv(ds.SampleIds[c]))));
        foreach (var (id, label) in features)
        {
            if (!rowOf.TryGetValue(id, out var row))
                continue;
            var cells = new List<string> { Csv(id), Csv(label) };
            foreach (var c in columns)
            {
                var v = ds.ExprLog2[row, c];
                cells.Add(double.IsFinite(v) ? Math.Pow(2.0, v).ToString("R", Inv) : "");
            }

            w.WriteLine(string.Join(",", cells));
        }
    }

    // -------- html helpers --------

    private static void AppendHitTable(
        StringBuilder sb, IEnumerable<DifferentialRow> rows, string effectName, Func<string, string> labelFor)
    {
        sb.Append("<table><tr><th>feature</th><th>").Append(HtmlEncode(effectName))
          .Append("</th><th>fold change</th><th>mean A</th><th>mean B</th><th>P</th><th>adj.P</th></tr>");
        foreach (var r in rows)
            sb.Append("<tr><td>").Append(HtmlEncode(labelFor(r.FeatureId)))
              .Append("</td><td>").Append(Num(r.LogFc, "0.###"))
              .Append("</td><td>").Append(Num(r.Fc, "0.###"))
              .Append("</td><td>").Append(Num(r.MeanA, "0.##"))
              .Append("</td><td>").Append(Num(r.MeanB, "0.##"))
              .Append("</td><td>").Append(Num(r.PValue, "0.##e0"))
              .Append("</td><td>").Append(Num(r.AdjPValue, "0.##e0"))
              .Append("</td></tr>");
        sb.Append("</table>");
    }

    private static void AppendDetectionTable(StringBuilder sb, IEnumerable<DetectionRow> rows)
    {
        sb.Append("<table><tr><th>peptide</th><th>rate A</th><th>rate B</th><th>P</th><th>adj.P</th></tr>");
        foreach (var r in rows)
            sb.Append("<tr><td>").Append(HtmlEncode(r.PeptideId))
              .Append("</td><td>").Append(Num(r.RateA, "0.00"))
              .Append("</td><td>").Append(Num(r.RateB, "0.00"))
              .Append("</td><td>").Append(Num(r.P, "0.##e0"))
              .Append("</td><td>").Append(Num(r.Q, "0.##e0"))
              .Append("</td></tr>");
        sb.Append("</table>");
    }

    private static void AppendEnrichmentTable(StringBuilder sb, IEnumerable<EnrichmentTerm> terms)
    {
        sb.Append("<table><tr><th>source</th><th>term</th><th>P</th><th>fold</th><th>genes</th></tr>");
        foreach (var t in terms)
            sb.Append("<tr><td>").Append(HtmlEncode(t.Source))
              .Append("</td><td>").Append(HtmlEncode(t.TermName))
              .Append("</td><td>").Append(Num(t.PValue, "0.##e0"))
              .Append("</td><td>").Append(Num(t.FoldEnrichment, "0.0"))
              .Append("</td><td style=\"text-align:left\">")
              .Append(HtmlEncode(GeneList(t.IntersectingGenes)))
              .Append("</td></tr>");
        sb.Append("</table>");
        sb.Append("<p class=\"note\">Gene lists are truncated here; the full members of each term are in "
            + "<code>enrichment_terms.csv</code>.</p>");
    }

    /// <summary>A term's genes for the HTML table: the first few, then a count, so an immunoglobulin
    /// term of 130 genes does not blow out the row. The full list is in the companion CSV.</summary>
    private static string GeneList(IReadOnlyList<string> genes, int show = 20)
    {
        if (genes.Count <= show)
            return string.Join(", ", genes);
        return string.Join(", ", genes.Take(show)) + $", ... (+{genes.Count - show} more)";
    }

    /// <summary>Note that a table is a preview and where the full rows live, shown only when it was capped.</summary>
    private static void AppendPreviewNote(StringBuilder sb, int total, int shown, string csv)
    {
        if (total > shown)
            sb.Append($"<p class=\"note\">Showing the top {shown} of {total}. The full table is in "
                + $"<code>quant/{HtmlEncode(csv)}</code>.</p>");
    }

    private static void AppendImage(StringBuilder sb, byte[] png, string alt) =>
        sb.Append("<div class=\"plot\"><img src=\"data:image/png;base64,")
          .Append(Convert.ToBase64String(png)).Append("\" alt=\"").Append(HtmlEncode(alt)).Append("\" /></div>");

    private static void AppendMessages(StringBuilder sb, string heading, IReadOnlyList<string> messages)
    {
        if (messages.Count == 0)
            return;
        sb.Append($"<div class=\"notes\"><strong>{heading}</strong><ul>");
        foreach (var m in messages)
            sb.Append("<li>").Append(HtmlEncode(m)).Append("</li>");
        sb.Append("</ul></div>");
    }

    private static void Kv(StringBuilder sb, string key, string value) =>
        sb.Append("<tr><td>").Append(HtmlEncode(key)).Append("</td><td>")
          .Append(HtmlEncode(value)).Append("</td></tr>");

    private static Provenance.RunInfo? ReadRunInfoSafe(string outputDir)
    {
        try { return Provenance.ReadRunInfo(Path.Combine(outputDir, "parameters.json")); }
        catch { return null; }
    }

    private static PrismConfig ReadPipelineConfigSafe(string outputDir)
    {
        var path = Path.Combine(outputDir, "parameters.json");
        try { return File.Exists(path) ? Provenance.LoadConfig(path) : new PrismConfig(); }
        catch { return new PrismConfig(); }
    }

    private static string Num(double v, string format) => double.IsNaN(v) ? "n/a" : v.ToString(format, Inv);

    private static string N(double v) => double.IsNaN(v) ? "" : v.ToString("R", Inv);

    private static string Csv(string s) =>
        s.IndexOfAny(new[] { ',', '"', '\n' }) >= 0 ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;

    private static string HtmlEncode(string s) => s
        .Replace("&", "&amp;").Replace("<", "&lt;").Replace(">", "&gt;").Replace("\"", "&quot;");

    private static void AppendStyle(StringBuilder sb)
    {
        sb.Append("<style>\nbody { font-family: ").Append(PlotRenderer.HtmlFontStack)
          .Append("; color: #222; margin: 0; padding: 24px; }\n");
        sb.Append("""
.container { max-width: 1400px; margin: 0 auto; }
h1 { color: #1a3c6e; }
h2 { color: #1a3c6e; border-bottom: 2px solid #dfe6ef; padding-bottom: 4px; margin-top: 32px; }
h3 { color: #1a3c6e; margin-top: 20px; }
.box { background: #f6f8fb; border: 1px solid #dfe6ef; border-radius: 6px; padding: 12px 16px; margin: 12px 0; }
table { border-collapse: collapse; margin: 8px 0; }
th, td { border: 1px solid #cfd8e3; padding: 6px 12px; text-align: right; }
th { background: #eaf0f7; }
td:first-child, th:first-child { text-align: left; }
table.kv td:nth-child(2) { text-align: left; }
.section-header { background: linear-gradient(90deg,#1a3c6e,#3a6ea5); color:#fff; padding:8px 14px; border-radius:6px; margin-top:24px; }
.plot img { max-width: 100%; border: 1px solid #dfe6ef; border-radius: 4px; margin: 8px 0; }
.note { color: #555; font-size: 13px; max-width: 900px; margin: 6px 0 12px; line-height: 1.45; }
.notes { background: #f2f6fc; border: 1px solid #cfd8e3; border-radius: 6px; padding: 8px 14px; margin: 8px 0; }
.notes li { color: #40536e; }
.footer { color: #888; font-size: 12px; margin-top: 32px; }
details { margin: 10px 0; }
summary { cursor: pointer; color: #1a3c6e; font-weight: 600; }
pre { background: #f6f8fb; border: 1px solid #dfe6ef; border-radius: 6px; padding: 10px 14px; overflow-x: auto; font-size: 12.5px; line-height: 1.4; }
</style>
""");
    }
}
