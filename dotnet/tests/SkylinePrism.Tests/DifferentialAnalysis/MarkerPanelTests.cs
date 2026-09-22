using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.Qc;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Tests for <see cref="MarkerPanel"/>: gene-symbol matching, row z-scoring, group-mean vs per-sample
/// heatmap columns, the per-group panel-score for the boxplot, and the found/not-detected report.
/// </summary>
public class MarkerPanelTests
{
    // Two matched markers (ALB, CD9) and one unmatched feature (XXX); the panel also names ZZZ, absent.
    private static readonly string[] Ids = { "PG1", "PG2", "PG3" };
    private static readonly string[] Labels = { "ALB", "CD9", "XXX" };
    private static readonly string?[] Groups = { "A", "A", "B", "B" };
    private static readonly string[] Samples = { "s1", "s2", "s3", "s4" };

    private static double[,] Expr() => new double[,]
    {
        { 1, 1, 3, 3 }, // ALB: mean 2, sd 1 -> z -1,-1,1,1
        { 1, 1, 3, 3 }, // CD9: same
        { 5, 6, 7, 8 }, // XXX: not in panel
    };

    private static ProteinList Panel() =>
        new() { Name = "Test panel", Members = { "ALB", "CD9", "ZZZ" } };

    /// <summary>
    /// The identity a protein row carries when the matrix named only a group id and a gene - the
    /// shape the existing cases were written against.
    /// </summary>
    private static FeatureIdentity[] Identities() => Ids
        .Select((id, i) => new FeatureIdentity(
            id, Labels[i], new[] { id }, System.Array.Empty<string>(),
            System.Array.Empty<string>(), new[] { Labels[i] }))
        .ToArray();

    [Fact]
    public void Evaluate_GroupMeans_MatchesExpected()
    {
        var r = MarkerPanel.Evaluate(Expr(), Identities(), Groups, Samples, Panel(), perSample: false);

        Assert.Equal(new[] { "ALB", "CD9" }, r.MarkerLabels);
        Assert.Equal(new[] { "A", "B" }, r.ColumnLabels);
        Assert.Equal(new[] { "A", "B" }, r.GroupNames);

        // Row z-scored group means: A = -1, B = +1 for both markers.
        Assert.Equal(-1.0, r.Heatmap[0, 0], 9);
        Assert.Equal(1.0, r.Heatmap[0, 1], 9);
        Assert.Equal(-1.0, r.Heatmap[1, 0], 9);
        Assert.Equal(1.0, r.Heatmap[1, 1], 9);

        Assert.Equal(2, r.Found);
        Assert.Equal(3, r.Total);
        Assert.Equal(new[] { "ZZZ" }, r.NotDetected.ToArray());
        Assert.Equal(1.0, r.SymmetricMax, 9);

        // Panel score per group: each sample's mean marker z-score.
        Assert.Equal(new[] { -1.0, -1.0 }, r.PanelScoreByGroup[0]);
        Assert.Equal(new[] { 1.0, 1.0 }, r.PanelScoreByGroup[1]);
    }

    [Fact]
    public void Evaluate_PerSample_ColumnsAreSamplesOrderedByGroup()
    {
        var r = MarkerPanel.Evaluate(Expr(), Identities(), Groups, Samples, Panel(), perSample: true);

        Assert.Equal(new[] { "s1", "s2", "s3", "s4" }, r.ColumnLabels);
        Assert.Equal(-1.0, r.Heatmap[0, 0], 9); // ALB, s1
        Assert.Equal(1.0, r.Heatmap[0, 3], 9);  // ALB, s4
    }

    [Fact]
    public void Evaluate_ZeroVarianceRow_IsNaN()
    {
        var expr = new double[,]
        {
            { 2, 2, 2, 2 }, // ALB: no variance -> NaN row
            { 1, 1, 3, 3 }, // CD9
            { 5, 6, 7, 8 },
        };

        var r = MarkerPanel.Evaluate(expr, Identities(), Groups, Samples, Panel(), perSample: false);
        Assert.True(double.IsNaN(r.Heatmap[0, 0]));
        Assert.True(double.IsNaN(r.Heatmap[0, 1]));
        Assert.Equal(-1.0, r.Heatmap[1, 0], 9);
    }

    /// <summary>
    /// A panel written in ACCESSIONS matches, which it did not when the matcher was handed only the
    /// feature id and its display label.
    /// </summary>
    [Fact]
    public void Evaluate_MatchesAProteinPanelWrittenInAccessions()
    {
        var identities = new[]
        {
            new FeatureIdentity("PG1", "ALB", new[] { "PG1" }, new[] { "P02768" },
                new[] { "ALBU_HUMAN" }, new[] { "ALB" }),
            new FeatureIdentity("PG2", "CD9", new[] { "PG2" }, new[] { "P21926" },
                new[] { "CD9_HUMAN" }, new[] { "CD9" }),
            new FeatureIdentity("PG3", "XXX", new[] { "PG3" }, new[] { "Q99999" },
                new[] { "XXX_HUMAN" }, new[] { "XXX" }),
        };
        var panel = new ProteinList { Name = "By accession", Members = { "P02768", "P21926" } };

        var r = MarkerPanel.Evaluate(Expr(), identities, Groups, Samples, panel, perSample: false);

        Assert.Equal(new[] { "ALB", "CD9" }, r.MarkerLabels);
        Assert.Equal(2, r.Found);
        Assert.Empty(r.NotDetected);
    }

    /// <summary>...and one written in protein NAMES.</summary>
    [Fact]
    public void Evaluate_MatchesAProteinPanelWrittenInNames()
    {
        var identities = new[]
        {
            new FeatureIdentity("PG1", "ALB", new[] { "PG1" }, new[] { "P02768" },
                new[] { "ALBU_HUMAN" }, new[] { "ALB" }),
            new FeatureIdentity("PG2", "CD9", new[] { "PG2" }, new[] { "P21926" },
                new[] { "CD9_HUMAN" }, new[] { "CD9" }),
            new FeatureIdentity("PG3", "XXX", new[] { "PG3" }, new[] { "Q99999" },
                new[] { "XXX_HUMAN" }, new[] { "XXX" }),
        };
        var panel = new ProteinList { Name = "By name", Members = { "ALBU_HUMAN" } };

        var r = MarkerPanel.Evaluate(Expr(), identities, Groups, Samples, panel, perSample: false);

        Assert.Equal(new[] { "ALB" }, r.MarkerLabels);
        Assert.Equal(1, r.Found);
    }

    /// <summary>
    /// A PEPTIDE row is labelled with its modified sequence, which is neither an accession nor a
    /// gene. It has to match on the protein columns the peptide matrix carries beside it.
    /// </summary>
    [Fact]
    public void Evaluate_MatchesPeptideRowsThroughTheirProteinColumns()
    {
        var identities = new[]
        {
            new FeatureIdentity("PEP1", "LVNEVTEFAK", new[] { "PG1" }, new[] { "P02768" },
                new[] { "ALBU_HUMAN" }, new[] { "ALB" }),
            new FeatureIdentity("PEP2", "SC(unimod:4)GLVGK", new[] { "PG2" }, new[] { "P21926" },
                new[] { "CD9_HUMAN" }, new[] { "CD9" }),
            new FeatureIdentity("PEP3", "QQTHPNGGEK", new[] { "PG3" }, new[] { "Q99999" },
                new[] { "XXX_HUMAN" }, new[] { "XXX" }),
        };
        var panel = new ProteinList { Name = "Genes", Members = { "ALB", "CD9", "ZZZ" } };

        var r = MarkerPanel.Evaluate(Expr(), identities, Groups, Samples, panel, perSample: false);

        // The two peptides of listed proteins match, and keep their sequences as labels.
        Assert.Equal(new[] { "LVNEVTEFAK", "SC(unimod:4)GLVGK" }, r.MarkerLabels);
        Assert.Equal(2, r.Found);
        Assert.Equal(3, r.Total);
        Assert.Equal(new[] { "ZZZ" }, r.NotDetected);
    }

    /// <summary>
    /// A SHARED peptide names every group it belongs to, and being in the list under any one of them
    /// is enough - the signal genuinely came from a protein the panel names.
    /// </summary>
    [Fact]
    public void Evaluate_MatchesASharedPeptideOnItsSecondGroup()
    {
        var identities = new[]
        {
            new FeatureIdentity("PEP1", "SHAREDPEPK", new[] { "PG9", "PG1" },
                new[] { "Q00000", "P02768" }, new[] { "OTHER_HUMAN", "ALBU_HUMAN" },
                new[] { "OTHER", "ALB" }),
            new FeatureIdentity("PEP2", "UNRELATEDK", new[] { "PG3" }, new[] { "Q99999" },
                new[] { "XXX_HUMAN" }, new[] { "XXX" }),
            new FeatureIdentity("PEP3", "ALSOUNRELATEDK", new[] { "PG4" }, new[] { "Q88888" },
                new[] { "YYY_HUMAN" }, new[] { "YYY" }),
        };
        var panel = new ProteinList { Name = "ALB only", Members = { "ALB" } };

        var r = MarkerPanel.Evaluate(Expr(), identities, Groups, Samples, panel, perSample: false);

        Assert.Equal(new[] { "SHAREDPEPK" }, r.MarkerLabels);
        Assert.Equal(1, r.Found);
    }

    /// <summary>
    /// A feature the matrix carried no protein annotation for is still findable by its own id, so a
    /// panel of protein-group ids keeps working on a matrix written before those columns existed.
    /// </summary>
    [Fact]
    public void Evaluate_StillMatchesOnTheFeatureIdWhenNothingElseIsKnown()
    {
        var identities = Ids
            .Select(id => new FeatureIdentity(
                id, "", System.Array.Empty<string>(), System.Array.Empty<string>(),
                System.Array.Empty<string>(), System.Array.Empty<string>()))
            .ToArray();
        var panel = new ProteinList { Name = "By group id", Members = { "PG1" } };

        var r = MarkerPanel.Evaluate(Expr(), identities, Groups, Samples, panel, perSample: false);

        Assert.Equal(new[] { "PG1" }, r.MarkerLabels);
        Assert.Equal(1, r.Found);
    }
}
