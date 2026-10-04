using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.Qc;
using SkylinePrism.Tests.TestSupport;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// A marker panel's CSV filename keeps only letters and digits from its name, so two panels whose names
/// differ in punctuation alone once resolved to the same path. StreamWriter truncates rather than
/// failing, so one panel's export silently became a copy of the other's while the report linked each
/// section to it. Panel names are free text, so this needs no exotic input - "EV markers" ships and
/// "EV-markers" is the obvious thing for a user to call their own version.
/// </summary>
public class QuantReportPanelFileNameTests
{
    private static string MiniOutput => Fixtures.Path2("mini", "e2e-sum", "output");

    private static readonly DifferentialOptions Options = new() { Prior = VariancePrior.Global };

    [Fact]
    public void PanelsWhoseNamesDifferOnlyInPunctuation_GetSeparateFiles()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);
        var types = ds.MetadataValues("sample_type");
        var experimental = Enumerable.Range(0, ds.SampleIds.Length)
            .Where(j => types[j] == "experimental").ToList();
        var half = experimental.Count / 2;
        var groupA = experimental.Take(half).ToList();
        var groupB = experimental.Skip(half).ToList();
        var res = Differential.Run(ds.ExprLog2, ds.FeatureIds, groupA, groupB, Options);

        var gene = ds.FeatureGenes.First(g => !string.IsNullOrEmpty(g));
        var identities = Enumerable.Range(0, ds.FeatureIds.Length).Select(ds.IdentityOf).ToArray();
        var groups = ds.MetadataValues("sample_type");

        MarkerPanelResult Evaluate(string name)
        {
            var panel = new ProteinList { Name = name };
            panel.Members.Add(gene);
            return MarkerPanel.Evaluate(ds.ExprLog2, identities, groups, ds.SampleIds, panel, false);
        }

        // Three names that all reduce to the same stem, plus one that collides with the SUFFIXED form
        // of the second - so the disambiguation has to register what it hands out, not just the base.
        var markers = new List<MarkerReportSection>
        {
            new("EV markers", "sample_type", Evaluate("EV markers")),
            new("EV-markers", "sample_type", Evaluate("EV-markers")),
            new("EV.markers", "sample_type", Evaluate("EV.markers")),
            new("EV markers 2", "sample_type", Evaluate("EV markers 2")),
        };

        var inputs = new QuantReportInputs
        {
            Differential = res,
            Rule = SignificanceRule.Default,
            Corrected = true,
            Contrast = "experimental split A vs B",
            EffectName = "log2FC",
            LabelFor = id => id,
            Options = Options,
            GroupBy = "sample_type",
            ALabel = "first half",
            BLabel = "second half",
            Dataset = ds,
            ContrastColumns = groupA.Concat(groupB).ToList(),
            Markers = markers,
        };

        var config = new QuantConfig(
            Level: "protein",
            Contrast: new QuantContrast("sample_type", new[] { "A" }, new[] { "B" }, null),
            Design: "TwoGroup", Test: "ModeratedT", Prior: res.VariancePrior, PriorUsed: res.VariancePrior, Correction: "BH",
            Covariates: res.CovariatesUsed, HitRule: "adj.P < 0.05",
            DetectionEnabled: false, DetectionQ: 0.01,
            EnrichmentEnabled: false, EnrichmentSources: Array.Empty<string>(),
            EnrichmentDirection: "both",
            MarkerPanels: markers.Select(m => m.PanelName).ToArray());

        var dir = Path.Combine(Path.GetTempPath(), $"prism-quant-{Guid.NewGuid():N}");
        Directory.CreateDirectory(dir);
        try
        {
            QuantReport.Write(dir, config, inputs);
            var quant = Path.Combine(dir, "quant");

            // One z-score file and one values file per panel - four of each, not one of each overwritten
            // three times.
            var zscores = Directory.GetFiles(quant, "markers_*_zscores.csv");
            var values = Directory.GetFiles(quant, "markers_*_values.csv");
            Assert.Equal(markers.Count, zscores.Length);
            Assert.Equal(markers.Count, values.Length);

            // A panel's two exports share a stem, so a reader can pair them.
            var zStems = zscores.Select(p => Path.GetFileName(p)!.Replace("_zscores.csv", "")).OrderBy(s => s);
            var vStems = values.Select(p => Path.GetFileName(p)!.Replace("_values.csv", "")).OrderBy(s => s);
            Assert.Equal(zStems, vStems);
        }
        finally
        {
            Directory.Delete(dir, recursive: true);
        }
    }
}
