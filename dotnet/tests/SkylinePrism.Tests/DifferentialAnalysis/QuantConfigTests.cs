using System.Collections.Generic;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Tests for <see cref="QuantConfig"/> serialization - the YAML shown in the quant report and written to
/// <c>quant_parameters.yaml</c>, and the JSON written to <c>quant_parameters.json</c>. These are a
/// contract with anyone reproducing an analysis, so the shape and the escaping are pinned.
/// </summary>
public class QuantConfigTests
{
    private static QuantConfig TwoGroup() => new(
        Level: "protein",
        Contrast: new QuantContrast("condition", "control", "disease", null),
        Design: "TwoGroup",
        Test: "ModeratedT",
        Prior: "IntensityTrend",
        Correction: "BH",
        Covariates: new[] { "age", "sex" },
        HitRule: "adj.P < 0.05, |log2FC| >= 1",
        DetectionEnabled: true,
        DetectionQ: 0.01,
        EnrichmentEnabled: true,
        EnrichmentSources: new[] { "GO:BP", "GO:MF" },
        EnrichmentDirection: "both",
        MarkerPanels: new[] { "EV markers" });

    [Fact]
    public void ToYaml_TwoGroup_HasContrastArmsAndSections()
    {
        var yaml = TwoGroup().ToYaml();

        Assert.Contains("level: protein", yaml);
        Assert.Contains("group_by: condition", yaml);
        Assert.Contains("group_a: control", yaml);
        Assert.Contains("group_b: disease", yaml);
        Assert.Contains("correction: BH", yaml);
        Assert.Contains("covariates: [age, sex]", yaml);
        Assert.Contains("q_threshold: 0.01", yaml);
        Assert.Contains("panels: [EV markers]", yaml);
        // A trend line must not appear for a two-arm contrast.
        Assert.DoesNotContain("trend_over", yaml);
    }

    [Fact]
    public void ToYaml_QuotesValuesThatCarryYamlSpecialCharacters()
    {
        // A g:Profiler source such as "GO:BP" contains a colon and must be quoted to stay valid YAML.
        var yaml = TwoGroup().ToYaml();
        Assert.Contains("\"GO:BP\"", yaml);
        Assert.Contains("\"GO:MF\"", yaml);
    }

    [Fact]
    public void ToYaml_Trend_UsesTrendOverAndDropsArms()
    {
        var cfg = TwoGroup() with { Contrast = new QuantContrast(null, null, null, "week") };
        var yaml = cfg.ToYaml();

        Assert.Contains("trend_over: week", yaml);
        Assert.DoesNotContain("group_by:", yaml);
    }

    [Fact]
    public void ToYaml_QuotesValuesThatLeadWithAYamlIndicator()
    {
        // A value beginning with a flow/indicator character would otherwise reparse as a list or an
        // alias rather than a string, so it must be quoted.
        var cfg = TwoGroup() with { MarkerPanels = new[] { "[draft] panel", "@internal" } };
        var yaml = cfg.ToYaml();
        Assert.Contains("\"[draft] panel\"", yaml);
        Assert.Contains("\"@internal\"", yaml);
    }

    [Fact]
    public void ToYaml_EmptyLists_RenderAsEmptyBrackets()
    {
        var cfg = TwoGroup() with { Covariates = new List<string>() };
        Assert.Contains("covariates: []", cfg.ToYaml());
    }

    [Fact]
    public void ToJson_IsSnakeCaseAndOmitsNullTrend()
    {
        var json = TwoGroup().ToJson();

        Assert.Contains("\"level\": \"protein\"", json);
        Assert.Contains("\"group_by\": \"condition\"", json);
        Assert.Contains("\"detection_enabled\": true", json);
        Assert.Contains("\"detection_q\": 0.01", json);
        // TrendOver is null on a two-arm contrast and is omitted rather than serialized as null.
        Assert.DoesNotContain("trend_over", json);
    }
}
