using System;
using System.Collections.Generic;
using System.Linq;

namespace SkylinePrism.App;

/// <summary>
/// The columns the QC pane's Group-by dropdown offers, and which of them it starts on.
/// </summary>
/// <remarks>
/// Its own class for the same reason as <see cref="VizNavigation"/>: pure state with no window in
/// it, so it can be tested without a dispatcher. The rule is worth pinning because every way of
/// getting it wrong is silent - a missing column is a grouping the user simply cannot pick, with
/// nothing said. Two cases in particular:
/// <list type="bullet">
/// <item>A run from the CLI has no <c>skyline-reports/</c>, so the Replicates report contributes
/// nothing and <c>batch</c> has to arrive from the run's own <c>sample_metadata.csv</c>. For a
/// batch-correction tool, grouping by batch is the first thing anyone reaches for.</item>
/// <item>A clinical CSV attached in the Differential pane must reach this list too, or attaching one
/// enriches the Volcano's covariates while leaving the PCA unable to color by any of them.</item>
/// </list>
/// </remarks>
internal static class QcGroupColumns
{
    /// <summary>The synthetic column the pane always offers, backed by sample_metadata.csv.</summary>
    public const string SampleType = "Sample Type";

    /// <summary>
    /// Every source's columns in priority order, de-duplicated, with <see cref="SampleType"/>
    /// guaranteed present.
    /// </summary>
    /// <param name="replicateColumns">
    /// Columns from the document's own Replicates report - first because that is what the person
    /// analyzing the data curated in Skyline. Empty on a CLI run.
    /// </param>
    /// <param name="extraColumns">
    /// Columns from the run's <c>sample_metadata.csv</c> and any joined clinical table.
    /// </param>
    public static List<string> Offer(
        IEnumerable<string> replicateColumns, IEnumerable<string> extraColumns)
    {
        var columns = new List<string>();
        foreach (var c in replicateColumns.Concat(extraColumns))
            if (!columns.Contains(c, StringComparer.Ordinal))
                columns.Add(c);
        // Always available: it is the fallback the pane is documented to default to, and it must not
        // disappear just because a directory happened to carry richer annotations.
        if (!columns.Any(IsSampleType))
            columns.Insert(0, SampleType);
        return columns;
    }

    /// <summary>
    /// The index <see cref="Offer"/>'s list should start on - Sample Type where it is, else the
    /// first column. Never -1, which would leave the pane grouping by nothing.
    /// </summary>
    public static int DefaultIndex(IReadOnlyList<string> columns)
    {
        for (var i = 0; i < columns.Count; i++)
            if (IsSampleType(columns[i]))
                return i;
        return 0;
    }

    /// <summary>
    /// Whether a column IS the sample-type grouping under any spelling. Matched ignoring spaces and
    /// case so a Skyline-exported "sample type" does not end up listed beside the synthetic
    /// "Sample Type" as two separate groupings of the same thing.
    /// </summary>
    private static bool IsSampleType(string column) =>
        column.Replace(" ", "").Equals("SampleType", StringComparison.OrdinalIgnoreCase);
}
