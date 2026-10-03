using System;
using System.Collections.Generic;
using System.IO;
using SkylinePrism.Core.IO;
using Xunit;

namespace SkylinePrism.Tests.IO;

/// <summary>
/// A replicate annotation may be called "Sample Type" - in clinical work it very often is, meaning
/// serum vs plasma - so the export carries TWO columns that both match the detector: Skyline's
/// built-in <c>SampleType</c> and the user's annotation. Binding to the wrong one is silent: a column
/// is found, every value maps, and every replicate comes out "experimental".
///
/// <para>Shaped on the run that found this: 96 replicates, built-in SampleType = Unknown x84 /
/// Standard x6 / Quality Control x6, annotation "Sample Type" = Serum x94 / Plasma x2. All 96 were
/// reported experimental, so ComBat lost its anchors and the QC report had no controls.</para>
/// </summary>
public class SampleTypeColumnCollisionTests
{
    private static string Write(string body)
    {
        var path = Path.Combine(Path.GetTempPath(), "rep_" + Guid.NewGuid().ToString("N") + ".csv");
        File.WriteAllText(path, body);
        return path;
    }

    /// <summary>The real column order: built-in first, the clinical annotation much later.</summary>
    private const string BuiltInFirst =
        "Replicate,SampleType,AnalyteConcentration,BRI Subject ID,Sample Type,Volume_mL\n" +
        "IER-CRLs-001-016,Standard,,S1,Serum,1.0\n" +
        "WBC-IBR-PS-001-006,Quality Control,,S2,Serum,1.0\n" +
        "TAM-001-002,Unknown,,S3,Plasma,1.0\n";

    [Fact]
    public void TheBuiltInColumnWins_NotTheClinicalAnnotation()
    {
        var path = Write(BuiltInFirst);
        try
        {
            var md = ReplicateMetadata.TryLoad(path);
            Assert.NotNull(md);

            Assert.Equal("reference", md!.TypeFor("IER-CRLs-001-016", "IER-CRLs-001-016"));
            Assert.Equal("qc", md.TypeFor("WBC-IBR-PS-001-006", "WBC-IBR-PS-001-006"));

            // The annotation is still carried through verbatim - it is real metadata, just not this.
            var values = md.ValuesFor("IER-CRLs-001-016", "IER-CRLs-001-016");
            Assert.Equal("Serum", values!["Sample Type"]);
            Assert.Equal("Standard", values["SampleType"]);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void TheBuiltInColumnWins_WhateverOrderTheColumnsAreIn()
    {
        // Order must not decide it: the same document with the annotation first has the same answer.
        var path = Write(
            "Replicate,Sample Type,SampleType\n" +
            "IER-CRLs-001-016,Serum,Standard\n" +
            "WBC-IBR-PS-001-006,Serum,Quality Control\n" +
            "TAM-001-002,Plasma,Unknown\n");
        try
        {
            var md = ReplicateMetadata.TryLoad(path);
            Assert.Equal("reference", md!.TypeFor("IER-CRLs-001-016", "IER-CRLs-001-016"));
            Assert.Equal("qc", md.TypeFor("WBC-IBR-PS-001-006", "WBC-IBR-PS-001-006"));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void TheCollisionIsReported()
    {
        var path = Write(BuiltInFirst);
        var log = new List<string>();
        try
        {
            ReplicateMetadata.TryLoad(path, log.Add);
            var line = Assert.Single(log, l => l.Contains("could be the sample type", StringComparison.Ordinal));
            Assert.Contains("using 'SampleType'", line, StringComparison.Ordinal);
            Assert.Contains("Ignoring 'Sample Type'", line, StringComparison.Ordinal);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void AnExplicitColumnStillWins_EvenAgainstTheBuiltIn()
    {
        // The override exists so a user can settle it themselves; detection must not second-guess it.
        var path = Write(BuiltInFirst);
        try
        {
            var md = ReplicateMetadata.TryLoad(path, sampleTypeColumn: "Sample Type");
            Assert.NotNull(md);
            // "Serum" is not a Skyline sample type, so it maps to experimental - which is what the
            // user asked for by naming that column.
            Assert.Equal("experimental", md!.TypeFor("IER-CRLs-001-016", "IER-CRLs-001-016"));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void OneColumnIsUnchanged_AndSaysNothing()
    {
        var path = Write(
            "Replicate,SampleType\n" +
            "Pool_A,Standard\n" +
            "Study_01,Unknown\n");
        var log = new List<string>();
        try
        {
            var md = ReplicateMetadata.TryLoad(path, log.Add);
            Assert.Equal("reference", md!.TypeFor("Pool_A", "Pool_A"));
            Assert.DoesNotContain(log, l => l.Contains("could be the sample type", StringComparison.Ordinal));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void WhenNeitherColumnSpeaksSkyline_ItSaysSoRatherThanPickingQuietly()
    {
        var path = Write(
            "Replicate,SampleType,Sample Type\n" +
            "A,Serum,Fasting\n" +
            "B,Plasma,Fed\n");
        var log = new List<string>();
        try
        {
            ReplicateMetadata.TryLoad(path, log.Add);
            var line = Assert.Single(log, l => l.Contains("could be the sample type", StringComparison.Ordinal));
            Assert.Contains("this is a guess", line, StringComparison.Ordinal);
            Assert.Contains("metadata.sample_type_column", line, StringComparison.Ordinal);
        }
        finally { File.Delete(path); }
    }
}
