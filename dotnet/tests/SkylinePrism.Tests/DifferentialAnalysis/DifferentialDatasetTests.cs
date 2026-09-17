using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Tests.TestSupport;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Tests for <see cref="DifferentialDataset"/> against the committed mini output fixture, including an
/// end-to-end check (load -> contrast -> moderated-t) matching the Python load_prism + differential.
/// </summary>
public class DifferentialDatasetTests
{
    private static string MiniOutput => Fixtures.Path2("mini", "e2e-sum", "output");

    private static void AssertRel(double expected, double actual, double rtol)
    {
        Assert.True(Math.Abs(actual - expected) / Math.Abs(expected) <= rtol,
            $"expected {expected:R}, actual {actual:R}");
    }

    [Fact]
    public void Load_Protein_MatchesLoadPrism()
    {
        var d = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);

        Assert.Equal(3, d.FeatureIds.Length);
        Assert.Equal(166, d.SampleIds.Length);
        Assert.Equal("protein_group", d.IdColumn);
        Assert.Equal("leading_gene_name", d.LabelColumn);
        Assert.Equal("PG0001", d.FeatureIds[0]);
        Assert.Equal("IRType-Plasma-201_078__@__mini_plate2", d.SampleIds[0]);
        Assert.Equal(5.821064533147193, d.ExprLog2[0, 0], 9); // linear -> log2

        Assert.Contains("sample_type", d.MetadataColumns);
        Assert.Contains("batch", d.MetadataColumns);
        Assert.Equal(138, d.MetadataValues("sample_type").Count(v => v == "experimental"));
    }

    [Fact]
    public void LoadThenDifferential_MatchesExplorerEndToEnd()
    {
        var d = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);

        // Experimental sample columns in matrix order, split first half vs second half.
        var types = d.MetadataValues("sample_type");
        var experimental = Enumerable.Range(0, d.SampleIds.Length)
            .Where(j => types[j] == "experimental").ToList();
        var half = experimental.Count / 2;
        var groupA = experimental.Take(half).ToList();
        var groupB = experimental.Skip(half).ToList();

        var res = Differential.Run(d.ExprLog2, d.FeatureIds, groupA, groupB);

        Assert.Equal(3, res.NFeaturesTested);
        Assert.Equal(0.2924877138903692, res.DfPrior, 9);

        Assert.Equal("PG0002", res.Rows[0].FeatureId); // smallest p
        var byId = res.Rows.ToDictionary(r => r.FeatureId);
        Assert.Equal(-0.0002532187140976061, byId["PG0002"].LogFc, 9);
        Assert.Equal(-0.9933318101885763, byId["PG0002"].T, 9);
        AssertRel(0.3223085567035089, byId["PG0002"].PValue, 1e-9);
        Assert.Equal(-0.06552109959229187, byId["PG0001"].LogFc, 9);
        AssertRel(0.6155448330123279, byId["PG0001"].PValue, 1e-9);
        AssertRel(0.78184737677963, byId["PG0003"].AdjPValue, 1e-9);
    }
}
