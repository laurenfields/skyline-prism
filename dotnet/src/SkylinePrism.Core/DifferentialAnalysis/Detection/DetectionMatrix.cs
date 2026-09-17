using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using DuckDB.NET.Data;
using SkylinePrism.Core.IO;

namespace SkylinePrism.Core.DifferentialAnalysis.Detection;

/// <summary>A binary peptide x sample detection matrix.</summary>
public sealed class DetectionMatrixData
{
    internal DetectionMatrixData(string[] peptideIds, string[] sampleIds, double[,] matrix)
    {
        PeptideIds = peptideIds;
        SampleIds = sampleIds;
        Matrix = matrix;
    }

    /// <summary>Peptide ids, sorted ascending (matrix rows).</summary>
    public string[] PeptideIds { get; }

    /// <summary>Sample ids, sorted ascending (matrix columns).</summary>
    public string[] SampleIds { get; }

    /// <summary>Detection indicator, <c>[peptide, sample]</c>: 1 if detected, else 0.</summary>
    public double[,] Matrix { get; }
}

/// <summary>
/// Builds the binary peptide x sample detection matrix from a PRISM run's transition-level
/// <c>merged_data</c> (a hash-partitioned parquet dataset), ported from the explorer's
/// <c>cryptic_detection_matrix</c> / <c>cryptic_peptide_map</c> / <c>cryptic_short_label</c>. A cell is
/// 1 iff the peptide was genuinely detected in that sample (<c>DetectionQValue</c> present and below the
/// threshold), the on/off signal the dense abundance matrix cannot give. Reads with DuckDB, filtering
/// server-side.
/// </summary>
public static class DetectionMatrix
{
    /// <summary>
    /// Load the detection matrix from the run's output directory or a <c>merged_data</c> root.
    /// <paramref name="term"/> restricts to peptides whose <c>Protein</c> contains it (the cryptic
    /// flag, e.g. "cryptic"); pass null to cover all peptides. A cell is 1 when
    /// <c>DetectionQValue &lt; qThreshold</c> in that sample.
    /// </summary>
    public static DetectionMatrixData Load(string outputDirOrMergedRoot, double qThreshold = 0.01,
        string? term = "cryptic")
    {
        var dataset = OpenDataset(outputDirOrMergedRoot);
        var thr = qThreshold.ToString(CultureInfo.InvariantCulture);
        var sql =
            "SELECT \"PeptideModifiedSequenceUnimodIds\" AS pep, \"Sample ID\" AS samp, " +
            $"MAX(CASE WHEN \"DetectionQValue\" IS NOT NULL AND \"DetectionQValue\" < {thr} " +
            "THEN 1 ELSE 0 END) AS det " +
            $"FROM {MergedParquetReader.Scan(dataset.ScanTarget)} " +
            (term is null
                ? "GROUP BY pep, samp"
                : $"WHERE \"Protein\" ILIKE '%{Esc(term)}%' GROUP BY pep, samp");

        var triples = new List<(string Pep, string Samp, int Det)>();
        using (var conn = OpenBounded(dataset))
        using (var cmd = DuckDbTuning.StreamingCommand(conn, sql))
        using (var reader = cmd.ExecuteReader())
        {
            while (reader.Read())
                triples.Add((reader.GetString(0), reader.GetString(1),
                    Convert.ToInt32(reader.GetValue(2), CultureInfo.InvariantCulture)));
        }

        var peptides = triples.Select(t => t.Pep).Distinct().OrderBy(p => p, StringComparer.Ordinal).ToArray();
        var samples = triples.Select(t => t.Samp).Distinct().OrderBy(s => s, StringComparer.Ordinal).ToArray();
        var pepIdx = peptides.Select((p, i) => (p, i)).ToDictionary(x => x.p, x => x.i);
        var sampIdx = samples.Select((s, i) => (s, i)).ToDictionary(x => x.s, x => x.i);

        var matrix = new double[peptides.Length, samples.Length];
        foreach (var (pep, samp, det) in triples)
            matrix[pepIdx[pep], sampIdx[samp]] = det;

        return new DetectionMatrixData(peptides, samples, matrix);
    }

    /// <summary>
    /// Map each peptide whose <c>Protein</c> contains <paramref name="term"/> to a representative such
    /// protein string (first occurrence wins). Cryptic peptides carry the flag in their Protein field,
    /// which the protein rollup discards, so this recovers it from the transition table.
    /// </summary>
    public static IReadOnlyDictionary<string, string> CrypticPeptideMap(string outputDirOrMergedRoot,
        string term = "cryptic")
    {
        var dataset = OpenDataset(outputDirOrMergedRoot);
        var sql =
            "SELECT DISTINCT \"PeptideModifiedSequenceUnimodIds\" AS pep, \"Protein\" AS prot " +
            $"FROM {MergedParquetReader.Scan(dataset.ScanTarget)} " +
            $"WHERE \"Protein\" ILIKE '%{Esc(term)}%'";

        var map = new Dictionary<string, string>();
        using var conn = OpenBounded(dataset);
        using var cmd = DuckDbTuning.StreamingCommand(conn, sql);
        using var reader = cmd.ExecuteReader();
        while (reader.Read())
        {
            var pep = reader.GetString(0);
            if (!map.ContainsKey(pep)) // setdefault: first occurrence wins
                map[pep] = reader.GetString(1);
        }

        return map;
    }

    /// <summary>
    /// Shorten a cryptic protein string to a readable name (e.g. "S35U4_HUMAN"): the first
    /// pipe-delimited part ending in "_HUMAN", else the second part, else the string itself.
    /// </summary>
    public static string CrypticShortLabel(string proteinString)
    {
        var parts = (proteinString ?? string.Empty).Split('|');
        foreach (var p in parts)
            if (p.EndsWith("_HUMAN", StringComparison.Ordinal))
                return p;
        return parts.Length > 1 ? parts[1] : proteinString ?? string.Empty;
    }

    private static MergedDataset OpenDataset(string outputDirOrMergedRoot)
    {
        // Accept either a run output directory (holding a merged_data child) or a merged_data root
        // (a directory of _pep_bucket partitions, or a legacy single file) passed directly.
        var root = MergedDataset.Locate(outputDirOrMergedRoot) ?? outputDirOrMergedRoot;
        return MergedDataset.Open(root);
    }

    private static DuckDBConnection OpenBounded(MergedDataset dataset)
    {
        var conn = new DuckDBConnection("Data Source=:memory:");
        conn.Open();
        DuckDbTuning.Apply(conn, DuckDbMerge.AutoMemoryBudgetMb(),
            DuckDbMerge.ResolveTempDirectory(dataset.Root));
        return conn;
    }

    private static string Esc(string s) => s.Replace("'", "''");
}
