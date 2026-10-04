using System.Collections.Generic;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Sample alignment used to choose between its two strategies for the WHOLE file - return the exact
/// matches if there were any, otherwise try bare replicate names - so a PARTIAL exact match discarded
/// the rest. The leftovers did not merely drop out: they fell through to the annotation-column list and
/// were treated as feature metadata, so a 96-sample contrast could run on 90 with nothing reported.
///
/// <para>Each column is now resolved on its own, which removes the failure instead of reporting it,
/// and whatever is still left over is refused by the caller.</para>
/// </summary>
public class SampleAlignmentPartialMatchTests
{
    private static readonly string[] MetaCols = { "sample", "sample_type" };

    [Fact]
    public void PartialExactMatch_StillResolvesTheStragglers()
    {
        // Two columns carry the metadata's own stem, one carries an older run's. Before, the two exact
        // matches were returned and the third was silently dropped from the analysis.
        var meta = new Dictionary<string, string?[]>
        {
            ["R1__@__merged_data"] = new string?[] { "R1", "experimental" },
            ["R2__@__merged_data"] = new string?[] { "R2", "qc" },
            ["R3__@__merged_data"] = new string?[] { "R3", "reference" },
        };
        var cols = new[]
        {
            "protein_group", "R1__@__merged_data", "R2__@__merged_data", "R3__@__PRISM",
        };

        var aligned = DifferentialDataset.AlignSampleColumns(cols, meta, MetaCols);

        Assert.Equal(
            new[] { "R1__@__merged_data", "R2__@__merged_data", "R3__@__PRISM" },
            aligned.Select(a => a.Col));
        Assert.Equal("reference", aligned.Single(a => a.Col == "R3__@__PRISM").Meta[1]);
        Assert.Empty(DifferentialDataset.UnmatchedMetadata(meta, aligned));
    }

    [Fact]
    public void ARowTakenByAnExactMatch_IsNotHandedToASecondColumn()
    {
        // The hazard the per-column fallback introduces: the matrix has two columns whose bare name is
        // "R1", and the metadata has only one R1 row. The exact match takes it; the other column must
        // NOT also take it, or two different samples share one sample's type, batch and subject.
        var meta = new Dictionary<string, string?[]>
        {
            ["R1__@__plateA"] = new string?[] { "R1", "reference" },
        };
        var cols = new[] { "R1__@__plateA", "R1__@__plateB" };

        var aligned = DifferentialDataset.AlignSampleColumns(cols, meta, MetaCols);

        var pair = Assert.Single(aligned);
        Assert.Equal("R1__@__plateA", pair.Col);
    }

    [Fact]
    public void UnmatchedMetadata_NamesTheRowsNoColumnResolved()
    {
        var meta = new Dictionary<string, string?[]>
        {
            ["R1__@__b"] = new string?[] { "R1", "experimental" },
            ["R2__@__b"] = new string?[] { "R2", "qc" },
            ["R3__@__b"] = new string?[] { "R3", "reference" },
        };
        // R3 has no column at all - neither its sample_id nor its bare name appears.
        var cols = new[] { "protein_group", "R1__@__b", "R2__@__b" };

        var aligned = DifferentialDataset.AlignSampleColumns(cols, meta, MetaCols);
        var missing = DifferentialDataset.UnmatchedMetadata(meta, aligned);

        Assert.Equal(new[] { "R3__@__b" }, missing);
    }

    [Fact]
    public void AnnotationColumnsAreNotCountedAsUnmatched()
    {
        // Leftover COLUMNS are normal - the matrix carries feature annotations that match no sample by
        // construction. Only leftover metadata ROWS are the error.
        var meta = new Dictionary<string, string?[]>
        {
            ["R1__@__b"] = new string?[] { "R1", "experimental" },
        };
        var cols = new[]
        {
            "protein_group", "leading_protein", "leading_gene_name", "n_peptides", "R1__@__b",
        };

        var aligned = DifferentialDataset.AlignSampleColumns(cols, meta, MetaCols);

        Assert.Single(aligned);
        Assert.Empty(DifferentialDataset.UnmatchedMetadata(meta, aligned));
    }

    [Fact]
    public void TwoIdenticalMetadataRows_AreTrackedAsTwoRows()
    {
        // Rows equal field for field are still two rows, so claiming is by reference. If they were
        // compared by value, the second column would look already-claimed and be dropped.
        var meta = new Dictionary<string, string?[]>
        {
            ["R1__@__b"] = new string?[] { "R1", "experimental" },
            ["R2__@__b"] = new string?[] { "R1", "experimental" },
        };
        var cols = new[] { "R1__@__b", "R2__@__b" };

        var aligned = DifferentialDataset.AlignSampleColumns(cols, meta, MetaCols);

        Assert.Equal(2, aligned.Count);
        Assert.Empty(DifferentialDataset.UnmatchedMetadata(meta, aligned));
    }
}
