using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.Json;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.DifferentialAnalysis.Enrichment;
using SkylinePrism.Core.Qc;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Tests for <see cref="QuantAnalysis"/>, the orchestration behind both the pane's Quant report button
/// and <c>prism differential --report</c>: that each design gets the detection test it calls for, that
/// every view that cannot run is named rather than failing the report, and that the recorded
/// parameters are the words the CLI accepts.
/// </summary>
public class QuantAnalysisTests
{
    /// <summary>A poster with no network behind it.</summary>
    private sealed class OfflinePoster : IJsonPoster
    {
        public JsonElement Post(string url, IReadOnlyDictionary<string, object?> payload) =>
            throw new HttpRequestException("No such host is known.");
    }

    private static string TempDir()
    {
        var dir = Path.Combine(Path.GetTempPath(), $"prism-quantrun-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        return dir;
    }

    private static QuantRequest TwoArm(string dir, DifferentialOptions options, IJsonPoster? poster = null,
        IReadOnlyList<ProteinList>? panels = null, string? markerGroupBy = null, bool cache = true)
    {
        var (ds, det, a, b) = DetectionAnalysisTests.Setup();
        return new QuantRequest
        {
            OutputDir = dir,
            Dataset = ds,
            Options = options,
            Rule = SignificanceRule.Default,
            GroupBy = "sample_type",
            GroupA = a,
            GroupB = b,
            ALevels = new[] { "even" },
            BLevels = new[] { "odd" },
            CachedDetection = cache ? det : null,
            EnrichmentPoster = poster,
            MarkerPanels = panels ?? Array.Empty<ProteinList>(),
            MarkerGroupBy = markerGroupBy,
        };
    }

    [Fact]
    public void TwoArm_WritesTheReport_NamesSkippedViews_AndRecordsCliWords()
    {
        var dir = TempDir();
        try
        {
            var request = TwoArm(dir, new DifferentialOptions { Prior = VariancePrior.Global });
            var r = QuantAnalysis.Run(request);

            Assert.True(File.Exists(r.HtmlPath));
            Assert.Equal(DetectionMethod.FisherExact, r.Detection!.Method);
            // The cached matrix is reused, not re-read from merged_data.
            Assert.Same(request.CachedDetection, r.DetectionMatrix);
            Assert.Contains("Enrichment skipped: not requested.", r.Notes);
            Assert.Contains("Markers omitted: no marker panels were selected.", r.Notes);

            // The recorded parameters are prism differential's own flag values.
            var yaml = File.ReadAllText(Path.Combine(dir, "quant", "quant_parameters.yaml"));
            Assert.Contains("level: protein", yaml);
            Assert.Contains("design: unpaired", yaml);
            Assert.Contains("test: moderated", yaml);
            Assert.Contains("correction: bh", yaml);
            Assert.Contains("group_a: [even]", yaml);
            // Requested prior as a --prior value; what ran recorded beside it.
            Assert.Contains("prior: global\n", yaml);
            Assert.Contains("prior_used: global", yaml);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void AComputedResult_IsReused_NotRecomputed()
    {
        var dir = TempDir();
        try
        {
            var options = new DifferentialOptions { Prior = VariancePrior.Global };
            var request = TwoArm(dir, options);
            var ds = request.Dataset;
            var computed = Differential.Run(ds.ExprLog2, ds.FeatureIds, request.GroupA, request.GroupB, options);

            var r = QuantAnalysis.Run(new QuantRequest
            {
                OutputDir = dir, Dataset = ds, Options = options, Rule = request.Rule,
                GroupBy = request.GroupBy, GroupA = request.GroupA, GroupB = request.GroupB,
                ALevels = request.ALevels, BLevels = request.BLevels, CachedDetection = request.CachedDetection,
                Differential = computed,
            });

            Assert.Same(computed, r.Differential);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    /// <summary>
    /// Arm labels are derived from the levels, never defaulted: an invented "A"/"B" would reach
    /// differential.csv's provenance header as a level that does not exist.
    /// </summary>
    [Fact]
    public void TwoArm_WithoutLevels_IsRefused()
    {
        var dir = TempDir();
        try
        {
            var request = TwoArm(dir, new DifferentialOptions { Prior = VariancePrior.Global });
            Assert.Throws<ArgumentException>(() => QuantAnalysis.Run(new QuantRequest
            {
                OutputDir = dir, Dataset = request.Dataset, Options = request.Options, Rule = request.Rule,
                GroupBy = request.GroupBy, GroupA = request.GroupA, GroupB = request.GroupB,
            }));
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Paired_RunsMcNemar_AndExportsOnlyTheMatchedSubjects()
    {
        var dir = TempDir();
        try
        {
            var (ds, _, a, b) = DetectionAnalysisTests.Setup();
            var options = new DifferentialOptions
            {
                Design = DifferentialDesign.Paired,
                SubjectLabels = DetectionAnalysisTests.Subjects(ds, a, b),
                Prior = VariancePrior.Global,
            };
            var r = QuantAnalysis.Run(TwoArm(dir, options));

            Assert.Equal(DetectionMethod.McNemarPaired, r.Detection!.Method);
            var header = File.ReadLines(Path.Combine(dir, "quant", "differential_values.csv")).First().Split(',');
            Assert.Equal(2 + 20 + 20, header.Length); // feature_id, label, then the 20 matched pairs
            var yaml = File.ReadAllText(Path.Combine(dir, "quant", "quant_parameters.yaml"));
            Assert.Contains("design: paired", yaml);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Covariates_RunTheFirthGlm()
    {
        var dir = TempDir();
        try
        {
            var (ds, _, _, _) = DetectionAnalysisTests.Setup();
            var options = new DifferentialOptions
            {
                Covariates = new[] { Covariate.FromMetadata("batch", ds.MetadataValues("batch")) },
                Prior = VariancePrior.Global,
            };
            var r = QuantAnalysis.Run(TwoArm(dir, options));

            Assert.Equal(DetectionMethod.FirthGlm, r.Detection!.Method);
            Assert.Contains("batch", DetectionAnalysis.Describe(r.Detection));
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Trend_SkipsDetectionAndRawValues_AndSaysSo()
    {
        var dir = TempDir();
        var clinical = Path.Combine(dir, "clinical.csv");
        try
        {
            var (ds, _, _, _) = DetectionAnalysisTests.Setup();
            var names = ds.MetadataValues("sample");
            using (var w = new StreamWriter(clinical))
            {
                w.WriteLine("PatientName,week");
                for (var i = 0; i < names.Length; i++)
                    w.WriteLine($"{names[i]},{i % 5}");
            }

            ds.AttachClinical(clinical);
            var r = QuantAnalysis.Run(new QuantRequest
            {
                OutputDir = dir,
                Dataset = ds,
                Options = new DifferentialOptions
                {
                    Design = DifferentialDesign.LinearTrend,
                    TrendColumn = "week",
                    Prior = VariancePrior.Global,
                },
                Rule = SignificanceRule.Default,
            });

            Assert.Null(r.Detection);
            Assert.Contains(r.Notes, n => n.StartsWith("Detection skipped", StringComparison.Ordinal));
            Assert.Contains(r.Notes, n => n.StartsWith("No differential_values.csv", StringComparison.Ordinal));
            Assert.False(File.Exists(Path.Combine(dir, "quant", "differential_values.csv")));
            Assert.StartsWith("# trend: week from 0 to 4",
                File.ReadLines(Path.Combine(dir, "quant", "differential.csv")).First());
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Offline_And_NoMergedData_BecomeNotes_NotFailures()
    {
        var dir = TempDir();
        try
        {
            // A rule every feature passes, so there ARE genes to submit and the poster is actually used.
            var request = TwoArm(dir, new DifferentialOptions { Prior = VariancePrior.Global },
                poster: new OfflinePoster(), cache: false);
            var r = QuantAnalysis.Run(new QuantRequest
            {
                OutputDir = request.OutputDir, Dataset = request.Dataset, Options = request.Options,
                GroupBy = request.GroupBy, GroupA = request.GroupA, GroupB = request.GroupB,
                ALevels = request.ALevels, BLevels = request.BLevels, EnrichmentPoster = request.EnrichmentPoster,
                Rule = new SignificanceRule { PThreshold = 1.1, Log2FcThreshold = 0 },
            });

            Assert.True(File.Exists(r.HtmlPath));
            Assert.Contains(r.Notes, n => n.StartsWith("Enrichment skipped (needs internet", StringComparison.Ordinal));
            // The temp dir has no merged_data.
            Assert.Contains(r.Notes, n => n.StartsWith("Detection skipped (needs merged_data", StringComparison.Ordinal));
            Assert.Null(r.Detection);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }

    [Fact]
    public void Markers_GroupByTheContrastColumnUnlessToldOtherwise()
    {
        var dir = TempDir();
        try
        {
            var (ds, _, _, _) = DetectionAnalysisTests.Setup();
            var panel = new ProteinList { Name = "Panel one" };
            panel.Members.Add(ds.FeatureGenes.First(g => !string.IsNullOrEmpty(g)));

            var r = QuantAnalysis.Run(TwoArm(dir, new DifferentialOptions { Prior = VariancePrior.Global },
                panels: new[] { panel }));
            Assert.DoesNotContain(r.Notes, n => n.StartsWith("Markers omitted", StringComparison.Ordinal));
            Assert.True(File.Exists(Path.Combine(dir, "quant", "markers_Panel_one_zscores.csv")));

            var bad = QuantAnalysis.Run(TwoArm(dir, new DifferentialOptions { Prior = VariancePrior.Global },
                panels: new[] { panel }, markerGroupBy: "no_such_column"));
            Assert.Contains("Markers omitted: no metadata column 'no_such_column' to group the panels by.", bad.Notes);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }
}
