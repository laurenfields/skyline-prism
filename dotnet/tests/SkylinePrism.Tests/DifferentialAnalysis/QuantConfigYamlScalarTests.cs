using System.Collections.Generic;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// <c>quant_parameters.yaml</c> is offered as a re-runnable record of the analysis, so every scalar in
/// it has to come back as the string that went in. Two ways a plain scalar fails to: it can RE-TYPE
/// (a case/control column spelled 0/1 reads back as integers, and no/yes/on/null as booleans and
/// nulls), and a quoted one can be mis-escaped (inside a YAML double-quoted scalar the backslash is an
/// escape introducer, so a value holding "\run" became a carriage return).
/// </summary>
public class QuantConfigYamlScalarTests
{
    private static QuantConfig WithArms(string a, string b, string groupBy = "condition") => new(
        Level: "protein",
        Contrast: new QuantContrast(groupBy, new[] { a }, new[] { b }, null),
        Design: "TwoGroup",
        Test: "ModeratedT",
        Prior: "Controls",
        PriorUsed: "Controls",
        Correction: "BH",
        Covariates: new string[0],
        HitRule: "adj.P < 0.05",
        DetectionEnabled: false,
        DetectionQ: 0.05,
        EnrichmentEnabled: false,
        EnrichmentSources: new string[0],
        EnrichmentDirection: "both",
        MarkerPanels: new string[0]);

    [Theory]
    // Numbers: a case/control column really is often spelled this way.
    [InlineData("0")]
    [InlineData("1")]
    [InlineData("1.0")]
    [InlineData("1e5")]
    [InlineData("-3")]
    [InlineData("0755")]
    [InlineData("0x1f")]
    // YAML 1.1 booleans and nulls.
    [InlineData("no")]
    [InlineData("NO")]
    [InlineData("yes")]
    [InlineData("on")]
    [InlineData("off")]
    [InlineData("y")]
    [InlineData("n")]
    [InlineData("null")]
    [InlineData("~")]
    public void AScalarThatWouldRetype_IsQuoted(string value)
    {
        var yaml = WithArms(value, "other").ToYaml();
        Assert.Contains($"group_a: [\"{value}\"]", yaml);
    }

    [Theory]
    [InlineData("control")]
    [InlineData("disease")]
    [InlineData("Group A1")]
    [InlineData("T1D")]
    public void AnOrdinaryScalar_IsLeftUnquoted(string value)
    {
        var yaml = WithArms(value, "other").ToYaml();
        Assert.Contains($"group_a: [{value}]\n", yaml);
    }

    [Fact]
    public void ABackslash_IsEscapedBeforeTheQuote()
    {
        // Escaping only the quote left \r to be read as a carriage return.
        var yaml = WithArms(@"C:\data\run", "other").ToYaml();
        Assert.Contains(@"group_a: [""C:\\data\\run""]", yaml);
    }

    [Fact]
    public void ControlCharacters_AreEscapedRatherThanEmittedRaw()
    {
        // A raw newline inside a double-quoted scalar ends the line and breaks the document.
        var yaml = WithArms("two\nlines\there", "other").ToYaml();
        Assert.Contains(@"group_a: [""two\nlines\there""]", yaml);
        // The value must not have broken the one-key-per-line shape.
        Assert.Contains("group_b: [other]\n", yaml);
    }

    [Fact]
    public void AQuote_IsStillEscaped()
    {
        var yaml = WithArms("say \"hi\"", "other").ToYaml();
        Assert.Contains(@"group_a: [""say \""hi\""""]", yaml);
    }

    [Fact]
    public void TheTrendColumn_FollowsTheSameRules()
    {
        // trend_over is the one scalar in the contrast block that used to bypass Yaml() entirely.
        var c = WithArms("a", "b") with { Contrast = new QuantContrast(null, null, null, "1.0") };
        Assert.Contains("trend_over: \"1.0\"", c.ToYaml());

        var ok = WithArms("a", "b") with { Contrast = new QuantContrast(null, null, null, "Age") };
        Assert.Contains("trend_over: Age\n", ok.ToYaml());
    }

    [Fact]
    public void ListItems_FollowTheSameRules()
    {
        var c = WithArms("control", "disease") with { Covariates = new[] { "0", "age", @"a\b" } };
        var yaml = c.ToYaml();
        Assert.Contains(@"covariates: [""0"", age, ""a\\b""]", yaml);
    }
}
