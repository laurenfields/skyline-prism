using System;
using System.IO;
using SkylinePrism.Core.IO;
using Xunit;

namespace SkylinePrism.Tests.IO;

/// <summary>
/// A Sample Type COLUMN is not a Sample Type ANNOTATION. Skyline's Replicates grid always offers the
/// column and defaults every replicate to "Unknown", so a document whose sample types were never set
/// exports a full column of unset cells. Recording a type for those rows makes
/// <c>sample_annotations.reference_pattern</c> / <c>qc_pattern</c> unreachable - the pipeline resolves
/// <c>metadata.TypeFor(...) ?? ClassifySampleType(...)</c>, and a non-null "experimental" never falls
/// through - so every Ref/QC injection is reported as experimental, ComBat loses its anchors and the QC
/// report has no controls to validate against.
/// </summary>
public class SampleTypeFallbackTests
{
    private static string WriteReport(string body)
    {
        var path = Path.Combine(Path.GetTempPath(), "rep_" + Guid.NewGuid().ToString("N") + ".csv");
        File.WriteAllText(path, body);
        return path;
    }

    [Theory]
    [InlineData("Unknown")]   // Skyline's default for an unset sample type
    [InlineData("")]          // an empty cell
    [InlineData("#N/A")]      // what Skyline exports for a value it has none of
    public void TryLoad_RecordsNoType_WhenTheCellIsUnset(string cell)
    {
        var path = WriteReport(
            "Replicate,Sample Type,Plate\n" +
            $"Ref_01,{cell},B1\n" +
            $"QC_01,{cell},B1\n" +
            $"Study_01,{cell},B1\n");
        try
        {
            var md = ReplicateMetadata.TryLoad(path, batchColumn: "Plate");
            Assert.NotNull(md);

            // No annotation means no entry, so the caller's name-pattern fallback still runs.
            Assert.False(md!.TypeByReplicate.ContainsKey("Ref_01"));
            Assert.Null(md.TypeFor("Ref_01", "Ref_01"));
            Assert.Null(md.TypeFor("QC_01", "QC_01"));
            Assert.Null(md.TypeFor("Study_01", "Study_01"));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Fact]
    public void TryLoad_KeepsRealTypes_AndStillSkipsTheUnsetRows()
    {
        // The mixed case is the one that must not regress into "all or nothing": a document where only
        // some replicates were annotated keeps those, and lets the rest fall through to the patterns.
        var path = WriteReport(
            "Replicate,Sample Type\n" +
            "Pool_A,Standard\n" +
            "Carl_A,Quality Control\n" +
            "Ref_02,Unknown\n" +
            "Study_01,\n");
        try
        {
            var md = ReplicateMetadata.TryLoad(path);
            Assert.NotNull(md);

            Assert.Equal("reference", md!.TypeFor("Pool_A", "Pool_A"));
            Assert.Equal("qc", md.TypeFor("Carl_A", "Carl_A"));
            Assert.Null(md.TypeFor("Ref_02", "Ref_02"));
            Assert.Null(md.TypeFor("Study_01", "Study_01"));

            // HasTypes stays true - the report really does carry annotations, just not for every row.
            Assert.True(md.HasTypes);
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Fact]
    public void TryLoad_StillClassifiesAnExplicitBlank()
    {
        // "Solvent"/"Blank"/"Double Blank" are REAL annotations that happen to mean "not a sample".
        // They must keep mapping to blank rather than being treated as unset and name-matched.
        var path = WriteReport(
            "Replicate,Sample Type\n" +
            "Blank_01,Blank\n" +
            "Solvent_01,Solvent\n");
        try
        {
            var md = ReplicateMetadata.TryLoad(path);
            Assert.NotNull(md);
            Assert.Equal("blank", md!.TypeFor("Blank_01", "Blank_01"));
            Assert.Equal("blank", md.TypeFor("Solvent_01", "Solvent_01"));
        }
        finally
        {
            File.Delete(path);
        }
    }
}
