using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>The two arms of a contrast, or a trend column - what the analysis compares.</summary>
public sealed record QuantContrast(string? GroupBy, string? GroupA, string? GroupB, string? TrendOver)
{
    /// <summary>A one-line description for the report and the status line.</summary>
    public string Describe() => TrendOver is not null
        ? $"trend over {TrendOver}"
        : $"{GroupB} vs {GroupA} by {GroupBy}";
}

/// <summary>
/// The parameters of one quantification analysis - what the quant report was produced from, and enough
/// to reproduce it. Serializes to YAML (shown in the report and written beside its outputs) and JSON
/// (machine-readable provenance), mirroring how the pipeline records its own config in
/// <c>parameters.json</c>. This is a DESCRIPTION built by the caller from its settings, not the
/// execution engine - the engine is <see cref="DifferentialOptions"/>.
/// </summary>
public sealed record QuantConfig(
    string Level,
    QuantContrast Contrast,
    string Design,
    string Test,
    string Prior,
    string Correction,
    IReadOnlyList<string> Covariates,
    string HitRule,
    bool DetectionEnabled,
    double DetectionQ,
    bool EnrichmentEnabled,
    IReadOnlyList<string> EnrichmentSources,
    string EnrichmentDirection,
    IReadOnlyList<string> MarkerPanels)
{
    private static readonly JsonSerializerOptions Json = new()
    {
        WriteIndented = true,
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
    };

    public string ToJson() => JsonSerializer.Serialize(this, Json);

    /// <summary>The parameters as a compact, human-readable YAML block - the "YAML like this" the
    /// report shows and writes to <c>quant_parameters.yaml</c>.</summary>
    public string ToYaml()
    {
        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.Append("level: ").Append(Level).Append('\n');
        sb.Append("contrast:\n");
        if (Contrast.TrendOver is not null)
            sb.Append("  trend_over: ").Append(Contrast.TrendOver).Append('\n');
        else
        {
            sb.Append("  group_by: ").Append(Yaml(Contrast.GroupBy)).Append('\n');
            sb.Append("  group_a: ").Append(Yaml(Contrast.GroupA)).Append('\n');
            sb.Append("  group_b: ").Append(Yaml(Contrast.GroupB)).Append('\n');
        }

        sb.Append("design: ").Append(Design).Append('\n');
        sb.Append("test: ").Append(Test).Append('\n');
        sb.Append("prior: ").Append(Prior).Append('\n');
        sb.Append("correction: ").Append(Correction).Append('\n');
        sb.Append("covariates: ").Append(YamlList(Covariates)).Append('\n');
        sb.Append("hit_rule: ").Append(Yaml(HitRule)).Append('\n');
        sb.Append("detection:\n");
        sb.Append("  enabled: ").Append(DetectionEnabled ? "true" : "false").Append('\n');
        sb.Append("  q_threshold: ").Append(DetectionQ.ToString("0.####", inv)).Append('\n');
        sb.Append("enrichment:\n");
        sb.Append("  enabled: ").Append(EnrichmentEnabled ? "true" : "false").Append('\n');
        sb.Append("  sources: ").Append(YamlList(EnrichmentSources)).Append('\n');
        sb.Append("  direction: ").Append(EnrichmentDirection).Append('\n');
        sb.Append("markers:\n");
        sb.Append("  panels: ").Append(YamlList(MarkerPanels)).Append('\n');
        return sb.ToString();
    }

    // Characters that, appearing anywhere, force quoting so the value cannot be read as structure.
    private static readonly char[] YamlSpecials = { ':', '#', '\'', '"', '\n', ',', '[', ']', '{', '}' };

    // Characters that only force quoting when they LEAD the value (a YAML indicator in first position
    // starts a flow collection, an alias/anchor, a tag, a block scalar, etc.).
    private const string YamlLeadIndicators = "-?:,[]{}#&*!|>'\"%@`";

    private static string Yaml(string? v)
    {
        if (string.IsNullOrEmpty(v))
            return "\"\"";
        var needsQuote = v.IndexOfAny(YamlSpecials) >= 0
            || YamlLeadIndicators.IndexOf(v[0]) >= 0
            || char.IsWhiteSpace(v[0]) || char.IsWhiteSpace(v[^1]);
        return needsQuote ? "\"" + v.Replace("\"", "\\\"") + "\"" : v;
    }

    private static string YamlList(IReadOnlyList<string> items) =>
        items.Count == 0 ? "[]" : "[" + string.Join(", ", items.Select(Yaml)) + "]";
}
