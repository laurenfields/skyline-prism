using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SkylinePrism.Core.IO;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>Feature level a differential dataset is loaded at.</summary>
public enum FeatureLevel
{
    Protein,
    Peptide,
}

/// <summary>
/// A loaded PRISM output ready for differential analysis, ported from the explorer's
/// <c>load_prism</c>. Reads the corrected (LINEAR-scale) protein or peptide matrix and
/// <c>sample_metadata.csv</c> from a run's output directory, identifies the sample columns by matching
/// the metadata's sample ids, and log2-transforms the abundances (non-positive values become NaN). The
/// sample metadata is exposed so a contrast can be built from any annotation column.
/// </summary>
public sealed class DifferentialDataset
{
    private readonly Dictionary<string, string?[]> _metaByColumn;

    private DifferentialDataset(FeatureLevel level, double[,] exprLog2, string[] featureIds,
        string[] featureLabels, string[] sampleIds, string idColumn, string labelColumn,
        IReadOnlyList<string> metadataColumns, Dictionary<string, string?[]> metaByColumn)
    {
        Level = level;
        ExprLog2 = exprLog2;
        FeatureIds = featureIds;
        FeatureLabels = featureLabels;
        SampleIds = sampleIds;
        IdColumn = idColumn;
        LabelColumn = labelColumn;
        MetadataColumns = metadataColumns;
        _metaByColumn = metaByColumn;
    }

    /// <summary>The level this was loaded at.</summary>
    public FeatureLevel Level { get; }

    /// <summary>LOG2 abundance matrix, <c>[feature, sample]</c> (NaN where the linear value was &lt;= 0).</summary>
    public double[,] ExprLog2 { get; }

    /// <summary>Feature identifiers (matrix rows).</summary>
    public string[] FeatureIds { get; }

    /// <summary>Human-readable feature labels (gene names for proteins); parallel to <see cref="FeatureIds"/>.</summary>
    public string[] FeatureLabels { get; }

    /// <summary>Sample identifiers (matrix columns), matching the metadata's sample ids.</summary>
    public string[] SampleIds { get; }

    /// <summary>The parquet column used as the feature id.</summary>
    public string IdColumn { get; }

    /// <summary>The parquet column used as the feature label.</summary>
    public string LabelColumn { get; }

    /// <summary>Metadata columns available for building a contrast (e.g. sample_type, batch, ...).</summary>
    public IReadOnlyList<string> MetadataColumns { get; }

    /// <summary>Metadata values for <paramref name="column"/>, aligned to <see cref="SampleIds"/>.</summary>
    public string?[] MetadataValues(string column) =>
        _metaByColumn.TryGetValue(column, out var v)
            ? v
            : throw new ArgumentException($"Unknown metadata column '{column}'.", nameof(column));

    /// <summary>
    /// Load the corrected matrix and sample metadata from a run's <paramref name="outputDir"/>.
    /// </summary>
    public static DifferentialDataset Load(string outputDir, FeatureLevel level)
    {
        var matrixName = level == FeatureLevel.Protein ? "corrected_proteins.parquet" : "corrected_peptides.parquet";
        var matrixPath = Path.Combine(outputDir, matrixName);
        var metaPath = Path.Combine(outputDir, "sample_metadata.csv");
        if (!File.Exists(matrixPath))
            throw new FileNotFoundException($"Corrected matrix not found: {matrixPath}", matrixPath);
        if (!File.Exists(metaPath))
            throw new FileNotFoundException($"sample_metadata.csv not found: {metaPath}", metaPath);

        var (sampleIdToRow, metaColumns) = ReadSampleMetadata(metaPath);

        var table = ParquetTable.Load(matrixPath);
        var sampleCols = table.ColumnNames.Where(sampleIdToRow.ContainsKey).ToArray();
        if (sampleCols.Length == 0)
            throw new InvalidOperationException(
                "No sample columns in the matrix matched sample_metadata.csv sample ids.");
        var annotCols = table.ColumnNames.Where(c => !sampleIdToRow.ContainsKey(c)).ToList();

        var idColumn = ChooseColumn(table, level == FeatureLevel.Protein
            ? new[] { "protein_group" }
            : new[] { "PeptideModifiedSequenceUnimodIds" }, annotCols);
        var labelColumn = level == FeatureLevel.Protein && table.HasColumn("leading_gene_name")
            ? "leading_gene_name"
            : idColumn;

        var nFeatures = table.RowCount;
        var featureIds = table.GetString(idColumn).Select(s => s ?? string.Empty).ToArray();
        var featureLabels = table.GetString(labelColumn).Select(s => s ?? string.Empty).ToArray();

        var exprLog2 = new double[nFeatures, sampleCols.Length];
        for (var j = 0; j < sampleCols.Length; j++)
        {
            var col = table.GetDouble(sampleCols[j]);
            for (var i = 0; i < nFeatures; i++)
            {
                var v = col[i];
                exprLog2[i, j] = v.HasValue && v.Value > 0.0 ? Math.Log2(v.Value) : double.NaN;
            }
        }

        // Metadata aligned to the sample-column order.
        var metaByColumn = new Dictionary<string, string?[]>(StringComparer.Ordinal);
        for (var m = 0; m < metaColumns.Count; m++)
        {
            var values = new string?[sampleCols.Length];
            for (var j = 0; j < sampleCols.Length; j++)
                values[j] = sampleIdToRow[sampleCols[j]][m];
            metaByColumn[metaColumns[m]] = values;
        }

        return new DifferentialDataset(level, exprLog2, featureIds, featureLabels, sampleCols,
            idColumn, labelColumn, metaColumns, metaByColumn);
    }

    private static string ChooseColumn(ParquetTable table, string[] preferred, List<string> fallback)
    {
        foreach (var name in preferred)
            if (table.HasColumn(name))
                return name;
        if (fallback.Count == 0)
            throw new InvalidOperationException("No annotation column available for the feature id.");
        return fallback[0];
    }

    private static (Dictionary<string, string?[]> ByRow, IReadOnlyList<string> Columns) ReadSampleMetadata(
        string path)
    {
        var lines = File.ReadAllLines(path);
        if (lines.Length == 0)
            throw new InvalidOperationException("sample_metadata.csv is empty.");

        var header = CsvLine.Split(lines[0]);
        var idIdx = CsvLine.IndexOf(header, "sample_id");
        if (idIdx < 0)
            throw new InvalidOperationException("sample_metadata.csv has no 'sample_id' column.");

        // Every column except sample_id, preserving file order.
        var otherIdx = Enumerable.Range(0, header.Length).Where(i => i != idIdx).ToArray();
        var columns = otherIdx.Select(i => header[i]).ToList();

        var byRow = new Dictionary<string, string?[]>(StringComparer.Ordinal);
        for (var r = 1; r < lines.Length; r++)
        {
            if (string.IsNullOrEmpty(lines[r]))
                continue;
            var fields = CsvLine.Split(lines[r]);
            if (idIdx >= fields.Length)
                continue;
            var sampleId = fields[idIdx];
            var values = new string?[otherIdx.Length];
            for (var k = 0; k < otherIdx.Length; k++)
            {
                // PRISM writes "" for a sample with no value; pandas reads that as NaN, and the
                // covariate/finder contracts treat missing as null - so an empty cell is null, not "".
                var raw = otherIdx[k] < fields.Length ? fields[otherIdx[k]] : null;
                values[k] = string.IsNullOrEmpty(raw) ? null : raw;
            }

            byRow[sampleId] = values;
        }

        return (byRow, columns);
    }
}
