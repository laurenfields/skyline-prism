using System;
using System.IO;
using System.Linq;
using SkylinePrism.App;
using SkylinePrism.Core.IO;
using Xunit;

namespace SkylinePrism.Tests.Windows;

/// <summary>
/// What the QC pane's Group-by dropdown offers. Every way of getting this wrong is silent: a column
/// that fails to arrive is simply a grouping the user cannot pick, with no error anywhere - so the
/// two sources that were added for it are pinned here rather than left to a look at the window.
/// </summary>
public class QcGroupColumnsTests
{
    private static string FixtureOutput =>
        Path.GetFullPath(Path.Combine(
            AppContext.BaseDirectory, "..", "..", "..", "..", "..",
            "tests", "fixtures", "mini", "e2e-sum", "output"));

    /// <summary>
    /// The case the sample_metadata.csv source exists for: a run from the CLI writes no
    /// <c>skyline-reports/</c>, so the Replicates report contributes nothing at all. Grouping by
    /// batch is the first thing anyone reaches for in a batch-correction tool, and before this it
    /// was unavailable on any run PRISM had not exported the reports for itself.
    /// </summary>
    [Fact]
    public void ACliRunWithNoReplicatesReport_StillOffersBatch()
    {
        var extra = SampleAnnotationTable.FromValues(
            new[] { "S1__@__plate1", "S2__@__plate2" },
            new[] { "batch" },
            (id, _) => id.EndsWith("plate1", StringComparison.Ordinal) ? "plate1" : "plate2");

        var columns = QcGroupColumns.Offer(Array.Empty<string>(), extra.Columns);

        Assert.Contains("batch", columns);
        Assert.Contains(QcGroupColumns.SampleType, columns);
    }

    /// <summary>
    /// The real thing, read off a committed golden: a CLI output directory with a sample_metadata.csv
    /// and no skyline-reports/ beside it. A hand-built table can agree with a wrong assumption about
    /// the file's shape; this cannot.
    /// </summary>
    [Fact]
    public void TheCommittedCliGolden_HasNoReplicatesReportAndYieldsBatch()
    {
        var dir = FixtureOutput;
        Assert.True(Directory.Exists(dir), $"fixture output missing: {dir}");
        // The premise of the test: this really is a run with no Replicates report.
        Assert.False(Directory.Exists(Path.Combine(dir, "skyline-reports")));

        var extra = SampleAnnotationTable.Read(
            Path.Combine(dir, "sample_metadata.csv"), "sample_id", "sample", "sample_type");

        var columns = QcGroupColumns.Offer(Array.Empty<string>(), extra.Columns);

        Assert.Contains("batch", columns);
        // sample and sample_type are excluded by the read: the pane already offers the latter as the
        // synthetic "Sample Type", and the same grouping under two spellings is only confusing.
        Assert.DoesNotContain("sample_type", columns);
        Assert.DoesNotContain("sample_id", columns);
    }

    /// <summary>
    /// Attaching a clinical CSV in the Differential pane must reach this list too. Otherwise it
    /// enriches the Volcano's covariates and leaves the PCA - which is now the only sample PCA in
    /// the window - unable to color by any of them.
    /// </summary>
    [Fact]
    public void ClinicalColumnsJoinedInTheDifferentialPane_BecomeGroupings()
    {
        var fromRun = SampleAnnotationTable.FromValues(
            new[] { "S1", "S2" }, new[] { "batch" }, (_, _) => "plate1");
        var clinical = SampleAnnotationTable.FromValues(
            new[] { "S1", "S2" }, new[] { "Diagnosis", "Age" }, (_, col) => col == "Age" ? "70" : "AD");

        var columns = QcGroupColumns.Offer(
            Array.Empty<string>(), fromRun.MergedWith(clinical).Columns);

        Assert.Contains("batch", columns);
        Assert.Contains("Diagnosis", columns);
        Assert.Contains("Age", columns);
    }

    /// <summary>
    /// Sample Type is the documented default, so it is offered even when nothing else supplies it,
    /// and it is what the pane starts on.
    /// </summary>
    [Fact]
    public void SampleTypeIsAlwaysOfferedAndIsTheDefault()
    {
        var columns = QcGroupColumns.Offer(Array.Empty<string>(), Array.Empty<string>());

        Assert.Equal(new[] { QcGroupColumns.SampleType }, columns);
        Assert.Equal(0, QcGroupColumns.DefaultIndex(columns));
    }

    /// <summary>
    /// A Skyline-exported "Sample Type" must not end up listed twice - once from the report and once
    /// synthetically - and the default must land on the real one rather than a duplicate.
    /// </summary>
    [Theory]
    [InlineData("Sample Type")]
    [InlineData("sample type")]
    [InlineData("SampleType")]
    public void AnExportedSampleTypeIsNotDuplicatedUnderAnySpelling(string exported)
    {
        var columns = QcGroupColumns.Offer(new[] { "Condition", exported }, new[] { "batch" });

        Assert.Single(columns, c =>
            c.Replace(" ", "").Equals("SampleType", StringComparison.OrdinalIgnoreCase));
        Assert.Equal(columns.IndexOf(exported), QcGroupColumns.DefaultIndex(columns));
    }

    /// <summary>
    /// The Replicates report comes first: it is what the person analyzing the data curated in
    /// Skyline, and a richer source should add columns without reordering it.
    /// </summary>
    [Fact]
    public void TheReplicatesReportKeepsPriorityOverTheRunsOwnColumns()
    {
        var columns = QcGroupColumns.Offer(new[] { "Condition" }, new[] { "batch" });

        Assert.Equal(new[] { QcGroupColumns.SampleType, "Condition", "batch" }, columns);
    }

    /// <summary>A column in both sources is listed once, not twice.</summary>
    [Fact]
    public void AColumnInBothSourcesIsListedOnce()
    {
        var columns = QcGroupColumns.Offer(new[] { "batch" }, new[] { "batch" });

        Assert.Single(columns, c => c == "batch");
    }
}
