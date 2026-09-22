using System;
using System.IO;
using System.Linq;
using SkylinePrism.Core.IO;
using Xunit;

namespace SkylinePrism.Tests.IO;

/// <summary>
/// The grouping-column source behind the QC plots' "Group by" dropdown. Its job is to widen what a
/// plot can be coloured by beyond the Skyline Replicates report - which a CLI-produced output
/// directory does not have - so the tests care most about degrading quietly and about not losing a
/// column to a comma inside a quoted value.
/// </summary>
public class SampleAnnotationTableTests : IDisposable
{
    private readonly string _dir =
        Path.Combine(Path.GetTempPath(), "prism_annot_" + Guid.NewGuid().ToString("N"));

    public SampleAnnotationTableTests() => Directory.CreateDirectory(_dir);

    public void Dispose()
    {
        try { Directory.Delete(_dir, recursive: true); } catch (IOException) { }
    }

    private string Write(string name, string content)
    {
        var path = Path.Combine(_dir, name);
        File.WriteAllText(path, content);
        return path;
    }

    [Fact]
    public void Read_SampleMetadata_OffersBatchButNotTheKeyOrExcludedColumns()
    {
        var path = Write("sample_metadata.csv",
            "sample_id,sample,sample_type,batch\n"
            + "R1__@__plate1,R1,experimental,plate1\n"
            + "R2__@__plate2,R2,qc,plate2\n");

        var t = SampleAnnotationTable.Read(path, "sample_id", "sample", "sample_type");

        // sample_type is excluded because the QC pane already offers it under "Sample Type" from
        // its own map; offering it twice under two spellings is the confusion this avoids.
        Assert.Equal(new[] { "batch" }, t.Columns);
        Assert.Equal("plate1", t.ValueOf("R1__@__plate1", "batch"));
        Assert.Equal("plate2", t.ValueOf("R2__@__plate2", "batch"));

        // Excluded from the offered columns, but still readable by a caller that knows the name.
        Assert.Equal("qc", t.ValueOf("R2__@__plate2", "sample_type"));
    }

    [Fact]
    public void ValueOf_UnknownSampleOrColumn_IsEmptyNotAnException()
    {
        var path = Write("m.csv", "sample_id,batch\nR1,plate1\n");
        var t = SampleAnnotationTable.Read(path, "sample_id");

        Assert.Equal("", t.ValueOf("nobody", "batch"));
        Assert.Equal("", t.ValueOf("R1", "no_such_column"));
    }

    [Theory]
    [InlineData("missing.csv", "file that is not there")]
    public void Read_MissingFile_IsEmpty(string name, string _)
    {
        Assert.True(SampleAnnotationTable.Read(Path.Combine(_dir, name), "sample_id").IsEmpty);
    }

    [Fact]
    public void Read_HeaderOnlyOrNoKeyColumn_IsEmpty()
    {
        Assert.True(SampleAnnotationTable.Read(
            Write("a.csv", "sample_id,batch\n"), "sample_id").IsEmpty);
        Assert.True(SampleAnnotationTable.Read(
            Write("b.csv", "something,else\n1,2\n"), "sample_id").IsEmpty);
    }

    [Fact]
    public void Read_KeyColumnMatchesIgnoringCaseSpacesAndUnderscores()
    {
        var path = Write("c.csv", "Sample ID,Braak\nR1,III\n");
        var t = SampleAnnotationTable.Read(path, "sample_id");

        Assert.Equal("III", t.ValueOf("R1", "Braak"));
    }

    [Fact]
    public void Read_QuotedFieldWithComma_StaysOneValue()
    {
        // A clinical table is hand-made often enough that this is not hypothetical, and a naive
        // Split(',') would shift every column after it by one - silently, and into a dropdown.
        var path = Write("d.csv",
            "sample_id,diagnosis,age\n"
            + "R1,\"Alzheimer's disease, sporadic\",81\n");

        var t = SampleAnnotationTable.Read(path, "sample_id");

        Assert.Equal("Alzheimer's disease, sporadic", t.ValueOf("R1", "diagnosis"));
        Assert.Equal("81", t.ValueOf("R1", "age"));
    }

    [Fact]
    public void Read_DoubledQuoteInsideQuotedField_IsOneQuote()
    {
        var path = Write("e.csv", "sample_id,note\nR1,\"said \"\"hello\"\"\"\n");
        Assert.Equal("said \"hello\"", SampleAnnotationTable.Read(path, "sample_id").ValueOf("R1", "note"));
    }

    [Fact]
    public void Read_ShortRow_KeepsWhatIsThereAndSkipsTheRest()
    {
        var path = Write("f.csv", "sample_id,batch,braak\nR1,plate1\n");
        var t = SampleAnnotationTable.Read(path, "sample_id");

        Assert.Equal("plate1", t.ValueOf("R1", "batch"));
        Assert.Equal("", t.ValueOf("R1", "braak"));
    }

    [Fact]
    public void MergedWith_UnionsColumnsAndThisTableWinsOnAClash()
    {
        var run = SampleAnnotationTable.Read(
            Write("g.csv", "sample_id,batch\nR1,from_run\n"), "sample_id");
        var clinical = SampleAnnotationTable.Read(
            Write("h.csv", "sample_id,batch,braak\nR1,from_clinical,IV\n"), "sample_id");

        var merged = run.MergedWith(clinical);

        Assert.Equal(new[] { "batch", "braak" }, merged.Columns);
        // The run's own metadata is authoritative about its own batch; the clinical file adds.
        Assert.Equal("from_run", merged.ValueOf("R1", "batch"));
        Assert.Equal("IV", merged.ValueOf("R1", "braak"));
    }

    [Fact]
    public void MergedWith_KeepsSamplesPresentInOnlyOneSide()
    {
        var a = SampleAnnotationTable.Read(Write("i.csv", "sample_id,x\nR1,1\n"), "sample_id");
        var b = SampleAnnotationTable.Read(Write("j.csv", "sample_id,y\nR2,2\n"), "sample_id");

        var merged = a.MergedWith(b);

        Assert.Equal("1", merged.ValueOf("R1", "x"));
        Assert.Equal("2", merged.ValueOf("R2", "y"));
        Assert.Equal("", merged.ValueOf("R1", "y"));
    }

    [Fact]
    public void MergedWith_Empty_IsIdentityInBothDirections()
    {
        var a = SampleAnnotationTable.Read(Write("k.csv", "sample_id,x\nR1,1\n"), "sample_id");

        Assert.Equal("1", a.MergedWith(SampleAnnotationTable.Empty).ValueOf("R1", "x"));
        Assert.Equal("1", SampleAnnotationTable.Empty.MergedWith(a).ValueOf("R1", "x"));
    }

    [Fact]
    public void FromValues_BuildsATableWithoutTouchingDisk()
    {
        var t = SampleAnnotationTable.FromValues(
            new[] { "R1", "R2" },
            new[] { "braak", "sex" },
            (id, col) => id == "R1" ? (col == "braak" ? "II" : "F") : null);

        Assert.Equal(new[] { "braak", "sex" }, t.Columns);
        Assert.Equal("II", t.ValueOf("R1", "braak"));
        Assert.Equal("F", t.ValueOf("R1", "sex"));
        Assert.Equal("", t.ValueOf("R2", "braak"));
    }

    [Fact]
    public void Read_RealMiniFixture_ExposesBatch()
    {
        var path = Path.Combine(
            AppContext.BaseDirectory, "fixtures", "mini", "e2e-sum", "output", "sample_metadata.csv");
        if (!File.Exists(path))
            return; // the fixture is only present for the end-to-end runs

        var t = SampleAnnotationTable.Read(path, "sample_id", "sample", "sample_type");

        Assert.Contains("batch", t.Columns);
        Assert.All(t.Samples, s => Assert.NotEqual("", t.ValueOf(s, "batch")));
    }
}
