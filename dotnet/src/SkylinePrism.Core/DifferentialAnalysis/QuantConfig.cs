using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>The two arms of a contrast, or a trend column - what the analysis compares.</summary>
/// <param name="GroupBy">The grouping column (two-arm contrasts).</param>
/// <param name="GroupA">Arm A's levels - a list, as <c>-a</c> takes one, so a union arm records as one.</param>
/// <param name="GroupB">Arm B's levels.</param>
/// <param name="TrendOver">The numeric column a trend is fitted against (trend designs).</param>
public sealed record QuantContrast(
    string? GroupBy, IReadOnlyList<string>? GroupA, IReadOnlyList<string>? GroupB, string? TrendOver)
{
    /// <summary>A one-line description for the report and the status line.</summary>
    public string Describe() => TrendOver is not null
        ? $"trend over {TrendOver}"
        : $"{ContrastArms.Describe(GroupB ?? Array.Empty<string>())} vs "
          + $"{ContrastArms.Describe(GroupA ?? Array.Empty<string>())} by {GroupBy}";
}

/// <summary>
/// The parameters of one quantification analysis - what the quant report was produced from, and enough
/// to reproduce it. Serializes to YAML (shown in the report and written beside its outputs) and JSON
/// (machine-readable provenance), mirroring how the pipeline records its own config in
/// <c>parameters.json</c>. This is a DESCRIPTION built by the caller from its settings, not the
/// execution engine - the engine is <see cref="DifferentialOptions"/>.
/// </summary>
/// <param name="Prior">The <c>--prior</c> value requested - what a re-run passes.</param>
/// <param name="PriorUsed">The prior that actually ran, with its source (<see cref="DifferentialResult.VariancePrior"/>,
/// e.g. "intensity-trend from controls"); it differs from the request when a requested prior had to fall back.</param>
public sealed record QuantConfig(
    string Level,
    QuantContrast Contrast,
    string Design,
    string Test,
    string Prior,
    string PriorUsed,
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
            sb.Append("  group_a: ").Append(YamlList(Contrast.GroupA ?? Array.Empty<string>())).Append('\n');
            sb.Append("  group_b: ").Append(YamlList(Contrast.GroupB ?? Array.Empty<string>())).Append('\n');
        }

        sb.Append("design: ").Append(Design).Append('\n');
        sb.Append("test: ").Append(Test).Append('\n');
        sb.Append("prior: ").Append(Prior).Append('\n');
        sb.Append("prior_used: ").Append(Yaml(PriorUsed)).Append('\n');
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
