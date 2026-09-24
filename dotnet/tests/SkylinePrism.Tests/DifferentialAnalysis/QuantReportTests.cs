using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.DifferentialAnalysis.Enrichment;
using SkylinePrism.Core.Qc;
using SkylinePrism.Tests.TestSupport;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Tests for <see cref="QuantReport"/> against the committed mini output fixture: it must write the
/// self-contained HTML plus its companion CSVs, cap long tables to a preview, and export the raw
/// per-sample abundances in LINEAR scale (2^log2), the same scale as PRISM's corrected parquet.
/// </summary>
public class QuantReportTests
{
    private static string MiniOutput => Fixtures.Path2("mini", "e2e-sum", "output");

    private static (DifferentialDataset Ds, DifferentialResult Res, List<int> Cols) RunContrast()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);
        var types = ds.MetadataValues("sample_type");
        var experimental = Enumerable.Range(0, ds.SampleIds.Length)
            .Where(j => types[j] == "experimental").ToList();
        var half = experimental.Count / 2;
        var groupA = experimental.Take(half).ToList();
        var groupB = experimental.Skip(half).ToList();
        var res = Differential.Run(ds.ExprLog2, ds.FeatureIds, groupA, groupB);
        return (ds, res, groupA.Concat(groupB).ToList());
    }

    private static QuantConfig ConfigFor(DifferentialResult res) => new(
        Level: "protein",
        Contrast: new QuantContrast("sample_type", "experimental (A)", "experimental (B)", null),
        Design: "TwoGroup",
        Test: "ModeratedT",
        Prior: res.VariancePrior,
        Correction: "BH",
        Covariates: res.CovariatesUsed,
        HitRule: "adj.P < 0.05, |log2FC| >= 1",
        DetectionEnabled: false,
        DetectionQ: 0.01,
        EnrichmentEnabled: false,
        EnrichmentSources: Array.Empty<string>(),
        EnrichmentDirection: "both",
        MarkerPanels: Array.Empty<string>());

    private static QuantReportInputs InputsFor(DifferentialDataset ds, DifferentialResult res,
        List<int> cols, SignificanceRule rule) => new()
    {
        Differential = res,
        Rule = rule,
        Corrected = true,
        Contrast = "experimental split A vs B",
        EffectName = "log2FC",
        LabelFor = id => id,
        Dataset = ds,
        ContrastColumns = cols,
    };

    [Fact]
    public void Write_ProducesReportAndCompanionFiles()
    {
        var (ds, res, cols) = RunContrast();
        var dir = Path.Combine(Path.GetTempPath(), $"prism-quant-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        try
        {
            // No parameters.json in this temp dir: the header falls back gracefully rather than throwing.
            var htmlPath = QuantReport.Write(dir, ConfigFor(res), InputsFor(ds, res, cols, SignificanceRule.Default));

            Assert.EndsWith("quant_report.html", htmlPath);
            var quant = Path.Combine(dir, "quant");
            Assert.True(File.Exists(htmlPath));
            Assert.True(File.Exists(Path.Combine(quant, "quant_parameters.yaml")));
            Assert.True(File.Exists(Path.Combine(quant, "quant_parameters.json")));
            Assert.True(File.Exists(Path.Combine(quant, "differential.csv")));
            Assert.True(File.Exists(Path.Combine(quant, "differential_values.csv")));

            // differential.csv: one header + one row per tested feature.
            var diffLines = File.ReadAllLines(Path.Combine(quant, "differential.csv"));
            Assert.Equal(res.Rows.Count + 1, diffLines.Length);
            Assert.StartsWith("feature_id,label,", diffLines[0]);

            // The HTML is self-contained: title, contrast, parameters block, and an embedded volcano PNG.
            var html = File.ReadAllText(htmlPath);
            Assert.Contains("PRISM Quantification Report", html);
            Assert.Contains("Quantification Parameters", html);
            Assert.Contains("experimental split A vs B", html);
            Assert.Contains("data:image/png;base64,", html);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Write_RawValueMatrix_IsLinearNotLog2()
    {
        var (ds, res, cols) = RunContrast();
        var dir = Path.Combine(Path.GetTempPath(), $"prism-quant-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        try
        {
            QuantReport.Write(dir, ConfigFor(res), InputsFor(ds, res, cols, SignificanceRule.Default));

            var lines = File.ReadAllLines(Path.Combine(dir, "quant", "differential_values.csv"));
            var header = lines[0].Split(',');
            // header: feature_id,label,<sample id of cols[0]>,...
            Assert.Equal("feature_id", header[0]);
            Assert.Equal(ds.SampleIds[cols[0]], header[2]);
            Assert.Equal(cols.Count + 2, header.Length);
            Assert.Equal(res.Rows.Count + 1, lines.Length);

            // Spot-check the first data row against the log2 matrix: the exported cell must be 2^log2,
            // i.e. the linear abundance, not the log2 value on disk in peptides_log2_internal.
            var rowOf = new Dictionary<string, int>(StringComparer.Ordinal);
            for (var i = 0; i < ds.FeatureIds.Length; i++)
                rowOf[ds.FeatureIds[i]] = i;

            var firstId = res.Rows[0].FeatureId;
            var dataLine = lines.Skip(1).First(l => l.StartsWith(firstId + ",", StringComparison.Ordinal));
            var actual = double.Parse(dataLine.Split(',')[2], CultureInfo.InvariantCulture);
            var expected = Math.Pow(2.0, ds.ExprLog2[rowOf[firstId], cols[0]]);
            Assert.Equal(expected, actual, 9);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Write_WithAllSections_RendersEachAndWritesItsCsv()
    {
        var (ds, res, cols) = RunContrast();

        // Detection and enrichment are the already-computed rows the report renders, so they can be
        // built directly here without a merged_data read or a g:Profiler call.
        var detection = new List<DetectionRow>
        {
            new("PEPTIDEA", 4, 5, 1, 5, 0.8, 0.2, 0.04, 0.09),
            new("PEPTIDEB", 2, 5, 3, 5, 0.4, 0.6, 0.5, 0.5),
        };
        var enrichment = new List<EnrichmentTerm>
        {
            new("GO:BP", "GO:0006915", "apoptotic process", 1e-4, 120, 50, 3, 20000, 11.5,
                new[] { "ALB", "CD9", "APOE" }),
        };

        // A real marker section from the fixture: pick a gene the dataset actually carries so the panel
        // finds at least one member and MarkerFeatureIds are real ids (which the values CSV needs).
        var gene = ds.FeatureGenes.First(g => !string.IsNullOrEmpty(g));
        var panel = new ProteinList { Name = "Test panel" };
        panel.Members.Add(gene);
        var identities = Enumerable.Range(0, ds.FeatureIds.Length).Select(ds.IdentityOf).ToArray();
        var groups = ds.MetadataValues("sample_type");
        var markerResult = MarkerPanel.Evaluate(ds.ExprLog2, identities, groups, ds.SampleIds, panel, false);
        Assert.True(markerResult.Found >= 1);
        var markers = new List<MarkerReportSection>
        {
            new("Test panel", "sample_type", markerResult),
        };

        var inputs = new QuantReportInputs
        {
            Differential = res,
            Rule = SignificanceRule.Default,
            Corrected = true,
            Contrast = "experimental split A vs B",
            EffectName = "log2FC",
            LabelFor = id => id,
            Detection = detection,
            Enrichment = enrichment,
            Markers = markers,
            Dataset = ds,
            ContrastColumns = cols,
        };

        var dir = Path.Combine(Path.GetTempPath(), $"prism-quant-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        try
        {
            var htmlPath = QuantReport.Write(dir, ConfigFor(res), inputs);
            var quant = Path.Combine(dir, "quant");

            // Detection: header + one row per detection row, and its HTML section rendered.
            var detCsv = File.ReadAllLines(Path.Combine(quant, "detection.csv"));
            Assert.Equal("peptide,detected_a,n_a,detected_b,n_b,rate_a,rate_b,p_value,adj_p_value", detCsv[0]);
            Assert.Equal(detection.Count + 1, detCsv.Length);

            // Enrichment: the full gene list reaches the CSV (the HTML truncates it).
            var enrCsv = File.ReadAllText(Path.Combine(quant, "enrichment_terms.csv"));
            Assert.Contains("apoptotic process", enrCsv);
            Assert.Contains("ALB;CD9;APOE", enrCsv);

            // Markers: both the z-score and the raw LINEAR value matrices, named after the panel.
            Assert.True(File.Exists(Path.Combine(quant, "markers_Test_panel_zscores.csv")));
            var mvals = File.ReadAllLines(Path.Combine(quant, "markers_Test_panel_values.csv"));
            Assert.True(mvals.Length >= 2); // header + at least the one matched member
            Assert.StartsWith("feature_id,label,", mvals[0]);

            var html = File.ReadAllText(htmlPath);
            Assert.Contains("Detection frequency", html);
            Assert.Contains("Functional enrichment", html);
            Assert.Contains("Marker panels", html);
            Assert.Contains("Test panel", html);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Write_CapsLongTablesToAPreview()
    {
        var (ds, res, cols) = RunContrast();
        // A permissive rule makes every tested feature a hit, so a cap of 1 must trip the preview note.
        var rule = new SignificanceRule { PThreshold = 1.1, Log2FcThreshold = 0.0, UseAdjusted = true };
        var dir = Path.Combine(Path.GetTempPath(), $"prism-quant-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        try
        {
            QuantReport.Write(dir, ConfigFor(res), InputsFor(ds, res, cols, rule), maxHitRows: 1);
            var html = File.ReadAllText(Path.Combine(dir, "quant", "quant_report.html"));
            Assert.Contains($"Showing the top 1 of {res.Rows.Count}", html);
            Assert.Contains("differential.csv", html);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }
}
