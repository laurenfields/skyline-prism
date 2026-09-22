using SkylinePrism.Core.IO;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using SkylinePrism.Cli;
using Xunit;

namespace SkylinePrism.Tests.Cli;

/// <summary>
/// End-to-end CLI coverage: drives Program.Main in-process against the mini fixture and asserts exit
/// codes + output files. Covers Program.cs (arg parsing + dispatch), the run/merge/qc/compare/
/// config-template commands, and RollupComparison - none of which the unit suite otherwise touches.
/// </summary>
[Collection("cli")]
public class CliIntegrationTests
{
    private static readonly object ConsoleLock = new();

    private static (int Code, string Output) Invoke(params string[] args)
    {
        lock (ConsoleLock)
        {
            var origOut = Console.Out;
            var origErr = Console.Error;
            using var sw = new StringWriter();
            Console.SetOut(sw);
            Console.SetError(sw);
            try
            {
                return (Program.Main(args), sw.ToString());
            }
            finally
            {
                Console.SetOut(origOut);
                Console.SetError(origErr);
            }
        }
    }

    private static string Fixture(params string[] parts)
        => Path.Combine(new[] { AppContext.BaseDirectory, "fixtures" }.Concat(parts).ToArray());

    private static readonly string Input1 = Fixture("mini", "merge", "mini_plate1.csv");
    private static readonly string Input2 = Fixture("mini", "merge", "mini_plate2.csv");
    private static readonly string Config = Fixture("mini", "e2e-sum", "config.yaml");

    private static string TempDir()
    {
        var d = Path.Combine(Path.GetTempPath(), "prism_cli_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(d);
        return d;
    }

    private static void Cleanup(string dir)
    {
        try { Directory.Delete(dir, recursive: true); }
        catch (IOException) { /* best-effort */ }
    }

    private static int Run(string outDir)
        => Invoke("run", "-i", Input1, Input2, "-o", outDir, "-c", Config).Code;

    [Fact]
    public void Run_ProducesExpectedOutputs()
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));
            // qc_report.html is gated by qc_report.enabled (off in this config); the qc command
            // test covers report generation.
            foreach (var f in new[]
            {
                "corrected_peptides.parquet", "corrected_proteins.parquet",
                "peptides_log2_internal.parquet", "proteins_raw.parquet",
                "protein_groups.csv", "parameters.json", "sample_metadata.csv",
            })
                Assert.True(File.Exists(Path.Combine(outDir, f)), $"missing output: {f}");
            Assert.Contains(Directory.GetFiles(outDir), p => Path.GetFileName(p).StartsWith("prism_run_"));
        }
        finally { Cleanup(outDir); }
    }

    [Fact]
    public void Merge_ProducesParquet()
    {
        var outDir = TempDir();
        // A ".parquet" -o is still accepted (it was the contract before the merge partitioned), but the
        // data lands in a directory of that name without the extension.
        var requested = Path.Combine(outDir, "merged.parquet");
        var actual = Path.Combine(outDir, "merged");
        try
        {
            var (code, output) = Invoke("merge", Input1, Input2, "-o", requested);
            Assert.Equal(0, code);
            Assert.False(File.Exists(requested));
            Assert.True(MergedDataset.Exists(actual), "merge produced no dataset");
            Assert.NotEmpty(MergedDataset.Open(actual).Partitions);
            Assert.Contains("rows", output);
        }
        finally { Cleanup(outDir); }
    }

    [Fact]
    public void Compare_ProducesReport()
    {
        var a = TempDir();
        var b = TempDir();
        var report = Path.Combine(b, "compare.html");
        try
        {
            Assert.Equal(0, Run(a));
            Assert.Equal(0, Run(b));
            var (code, _) = Invoke("compare", "-1", a, "-2", b, "-o", report, "-s", "all", "-n", "5");
            Assert.Equal(0, code);
            Assert.True(File.Exists(report));
            Assert.Contains("Rollup Comparison", File.ReadAllText(report));
        }
        finally { Cleanup(a); Cleanup(b); }
    }

    [Fact]
    public void Qc_RegeneratesReport()
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));
            var html = Path.Combine(outDir, "qc_report.html");
            File.Delete(html);
            Assert.Equal(0, Invoke("qc", "-d", outDir).Code);
            Assert.True(File.Exists(html));
        }
        finally { Cleanup(outDir); }
    }

    [Fact]
    public void ConfigTemplate_FullAndMinimal_WriteFiles()
    {
        var outDir = TempDir();
        try
        {
            var full = Path.Combine(outDir, "full.yaml");
            var min = Path.Combine(outDir, "min.yaml");
            Assert.Equal(0, Invoke("config-template", "-o", full).Code);
            Assert.Equal(0, Invoke("config-template", "--minimal", "-o", min).Code);
            Assert.Contains("transition_rollup", File.ReadAllText(full));
            // Minimal is a strict subset of the full template.
            Assert.True(new FileInfo(min).Length < new FileInfo(full).Length);
        }
        finally { Cleanup(outDir); }
    }

    [Fact]
    public void Version_PrintsVersion()
    {
        var (code, output) = Invoke("--version");
        Assert.Equal(0, code);
        Assert.Contains("prism", output);
    }

    [Theory]
    [InlineData(new[] { "bogus" }, 2)]                       // unknown command
    [InlineData(new[] { "run", "-o", "x" }, 2)]              // run without -i
    [InlineData(new[] { "merge", "-o", "x.parquet" }, 2)]    // merge without inputs
    [InlineData(new[] { "qc" }, 2)]                          // qc without -d
    [InlineData(new[] { "compare", "-1", "a" }, 2)]          // compare without -2
    public void InvalidInvocations_ReturnUsageCode(string[] args, int expected)
        => Assert.Equal(expected, Invoke(args).Code);

    [Fact]
    public void NoArgs_PrintsUsage()
    {
        var (code, output) = Invoke();
        Assert.Equal(0, code);
        Assert.Contains("Usage", output);
    }

    [Theory]
    [InlineData("run")]
    [InlineData("merge")]
    [InlineData("qc")]
    [InlineData("compare")]
    [InlineData("config-template")]
    public void CommandHelp_ReturnsZeroWithCommandUsage(string cmd)
    {
        var (code, output) = Invoke(cmd, "--help");
        Assert.Equal(0, code);
        Assert.Contains($"Usage: prism {cmd}", output);
    }

    /// <summary>
    /// An input edited in place is caught, even though not one setting moved.
    /// </summary>
    /// <remarks>
    /// This is the direction a config comparison cannot see at all, and it is the one that matters
    /// most: the YAML is identical, so comparing it says "nothing will change" while the merge and
    /// everything below it is about to be recomputed from different data. The merge is checked the
    /// way the pipeline checks it - against its own sidecar beside <c>merged_data</c>, which stamps
    /// each input's path, size and write time.
    /// </remarks>
    [Fact]
    public void Run_WarnsWhenAnInputChangedUnderUnchangedSettings()
    {
        var outDir = TempDir();
        var inputDir = TempDir();
        try
        {
            var a = Path.Combine(inputDir, "plate1.csv");
            var b = Path.Combine(inputDir, "plate2.csv");
            File.Copy(Input1, a);
            File.Copy(Input2, b);

            Assert.Equal(0, Invoke("run", "-i", a, b, "-o", outDir, "-c", Config).Code);

            // Same settings, same inputs: silent.
            var repeat = Invoke("run", "-i", a, b, "-o", outDir, "-c", Config);
            Assert.Equal(0, repeat.Code);
            Assert.DoesNotContain("already holds results", repeat.Output, StringComparison.Ordinal);

            // The same bytes, written again - which is what an export re-run looks like. Nothing in
            // the config has moved; only the file's stamp has.
            File.WriteAllBytes(a, File.ReadAllBytes(a));
            File.SetLastWriteTimeUtc(a, DateTime.UtcNow.AddSeconds(5));

            var changed = Invoke("run", "-i", a, b, "-o", outDir, "-c", Config);
            Assert.Equal(0, changed.Code);
            Assert.Contains("already holds results", changed.Output, StringComparison.Ordinal);

            // Named as what it is. Not one setting moved, so reporting this as "different settings"
            // would send the reader to a config diff that shows nothing at all.
            Assert.Contains("input files that have changed", changed.Output, StringComparison.Ordinal);
            Assert.DoesNotContain("different settings", changed.Output, StringComparison.Ordinal);
        }
        finally
        {
            Cleanup(outDir);
            Cleanup(inputDir);
        }
    }

    /// <summary>
    /// `prism differential` end to end: run the pipeline, then ask a contrast of its output.
    /// </summary>
    /// <remarks>
    /// One pipeline run serves every assertion below - it is about a second on the mini fixture, but
    /// a contrast is microseconds, so paying for it once per fact would be most of the test time.
    /// </remarks>
    [Fact]
    public void Differential_TestsAContrastAgainstAFinishedRun()
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));

            var csv = Path.Combine(outDir, "differential.csv");
            var (code, output) = Invoke(
                "differential", "-d", outDir, "--group-by", "sample_type", "-a", "qc", "-b", "experimental");

            Assert.Equal(0, code);
            Assert.True(File.Exists(csv), "the default results path is inside the output directory");

            // The status line names the method that produced the numbers - with a menu this size it
            // is the only thing that makes a saved console log interpretable.
            Assert.Contains("moderated t (intensity-trend prior), unpaired", output, StringComparison.Ordinal);
            Assert.Contains("sample_type = experimental vs qc", output, StringComparison.Ordinal);
            Assert.Contains("Benjamini-Hochberg", output, StringComparison.Ordinal);

            var lines = File.ReadAllLines(csv);

            // The header outlives the shell it was produced in, so it carries the contrast direction.
            Assert.StartsWith("# contrast: sample_type = experimental vs qc", lines[0], StringComparison.Ordinal);
            Assert.Contains("positive log2FC is higher in experimental", lines[0], StringComparison.Ordinal);
            Assert.StartsWith("# method:", lines[1], StringComparison.Ordinal);
            // The rule is recorded, and recorded as NOT having filtered the rows - every tested
            // feature is in the file, so a reader must not take the header as a description of
            // which rows survived.
            Assert.Contains(lines, l => l.StartsWith("# hit rule", StringComparison.Ordinal)
                                        && l.Contains("NOT filtered", StringComparison.Ordinal));

            // Found by its content rather than its index, so adding another comment line is not a
            // test failure - the index is what broke when the hit rule was added.
            var headerIndex = Array.FindIndex(lines, l => l.StartsWith("feature_id,", StringComparison.Ordinal));
            Assert.Equal(
                "feature_id,label,gene,protein,accession,log2fc,fc,ave_expr,statistic,"
                + "p_value,adj_p_value,mean_a,mean_b",
                lines[headerIndex]);

            var rows = lines.Skip(headerIndex + 1).Where(l => l.Length > 0).ToList();
            Assert.NotEmpty(rows);

            // Rows come out most significant first, and the p-values are real numbers in invariant
            // culture - a decimal comma would silently shift every column one to the right.
            var pValues = rows
                .Select(r => double.Parse(r.Split(',')[9], CultureInfo.InvariantCulture))
                .ToList();
            Assert.Equal(pValues.OrderBy(v => v), pValues);
            Assert.All(pValues, v => Assert.InRange(v, 0.0, 1.0));
        }
        finally
        {
            Cleanup(outDir);
        }
    }

    /// <summary>
    /// Each arm is the UNION of its levels, so a level added to an arm can only grow it.
    /// </summary>
    [Fact]
    public void Differential_PoolsSeveralLevelsIntoOneArm()
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));

            var (oneCode, one) = Invoke(
                "differential", "-d", outDir, "-g", "sample_type", "-a", "qc", "-b", "experimental",
                "-o", Path.Combine(outDir, "one.csv"));
            var (bothCode, both) = Invoke(
                "differential", "-d", outDir, "-g", "sample_type", "-a", "qc,reference",
                "-b", "experimental", "-o", Path.Combine(outDir, "both.csv"));

            Assert.Equal(0, oneCode);
            Assert.Equal(0, bothCode);
            Assert.Contains("= experimental vs qc", one, StringComparison.Ordinal);
            Assert.Contains("= experimental vs qc + reference", both, StringComparison.Ordinal);

            // Same arm B, larger arm A: the A count must rise and B must not move.
            var (aOne, bOne) = ArmCounts(one);
            var (aBoth, bBoth) = ArmCounts(both);
            Assert.True(aBoth > aOne, $"pooling should grow arm A ({aOne} -> {aBoth})");
            Assert.Equal(bOne, bBoth);
        }
        finally
        {
            Cleanup(outDir);
        }
    }

    /// <summary>The "n = A vs B" the status line reports.</summary>
    private static (int A, int B) ArmCounts(string output)
    {
        var m = Regex.Match(output, @"n = (\d+) vs (\d+)");
        Assert.True(m.Success, "the status line should report both arm sizes");
        return (int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture),
                int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture));
    }

    /// <summary>
    /// The refusals. Each one is a thing a user can type that would otherwise produce a plausible
    /// answer to a different question than the one asked.
    /// </summary>
    [Theory]
    // A level on both sides puts the same samples on both sides of the contrast.
    [InlineData("both arms", "-g", "sample_type", "-a", "qc,experimental", "-b", "experimental")]
    // A typo in a level name, named by the arm it emptied.
    [InlineData("Arm B matched no samples", "-g", "sample_type", "-a", "qc", "-b", "nosuchlevel")]
    // A column that is not in the metadata at all.
    [InlineData("No metadata column", "-g", "nosuchcolumn", "-a", "qc", "-b", "experimental")]
    // Paired without the column that says which samples are a pair.
    [InlineData("needs --pair-by", "-g", "sample_type", "-a", "qc", "-b", "experimental", "--design", "paired")]
    // ... and the reverse, which would otherwise report an unpaired result for a paired-looking command.
    [InlineData("--pair-by needs --design paired", "-g", "sample_type", "-a", "qc", "-b", "experimental",
        "--pair-by", "batch")]
    // A covariate handed to a test with no design matrix to put it in.
    [InlineData("--adjust-for needs --test moderated", "-g", "sample_type", "-a", "qc", "-b", "experimental",
        "--test", "welch", "--adjust-for", "batch")]
    // Adjusting for the contrast itself leaves nothing to test.
    [InlineData("is the contrast itself", "-g", "sample_type", "-a", "qc", "-b", "experimental",
        "--adjust-for", "sample_type")]
    [InlineData("Unknown --test", "-g", "sample_type", "-a", "qc", "-b", "experimental", "--test", "ttest")]
    [InlineData("Unknown --correction", "-g", "sample_type", "-a", "qc", "-b", "experimental",
        "--correction", "fdr")]
    [InlineData("--level must be", "-g", "sample_type", "-a", "qc", "-b", "experimental", "--level", "gene")]
    public void Differential_RefusesAndSaysWhy(string expected, params string[] args)
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));

            var (code, output) = Invoke(new[] { "differential", "-d", outDir }.Concat(args).ToArray());

            Assert.NotEqual(0, code);
            Assert.Contains(expected, output, StringComparison.Ordinal);
            Assert.False(File.Exists(Path.Combine(outDir, "differential.csv")),
                "a refused contrast must not leave a results file behind");
        }
        finally
        {
            Cleanup(outDir);
        }
    }

    [Fact]
    public void Differential_WithoutTheRequiredFlags_PrintsUsage()
    {
        var (code, output) = Invoke("differential", "-d", "nowhere");

        Assert.Equal(2, code);
        Assert.Contains("Usage: prism differential", output, StringComparison.Ordinal);
    }

    [Fact]
    public void Differential_IsListedAndDocumented()
    {
        var usage = Invoke("--help").Output;
        Assert.Contains("differential", usage, StringComparison.Ordinal);

        var help = Invoke("differential", "--help").Output;
        Assert.Contains("--group-by", help, StringComparison.Ordinal);
        Assert.Contains("--prior-from-controls", help, StringComparison.Ordinal);
        Assert.Contains("--correction", help, StringComparison.Ordinal);
    }

    /// <summary>
    /// Running onto a directory that already holds results says so - and only when the results would
    /// actually differ, which is what keeps the warning worth reading.
    /// </summary>
    /// <remarks>
    /// The decision itself is covered by <c>ExistingResultsTests</c>; this drives the real command
    /// three times so the wiring is proven end to end rather than by inspection. A run is about a
    /// second on the mini fixture.
    /// </remarks>
    [Fact]
    public void Run_WarnsOnlyWhenItWouldReplaceDifferentResults()
    {
        var outDir = TempDir();
        try
        {
            Assert.Equal(0, Run(outDir));

            var (repeatCode, repeatOutput) = Invoke("run", "-i", Input1, Input2, "-o", outDir, "-c", Config);
            Assert.Equal(0, repeatCode);
            Assert.DoesNotContain("already holds results", repeatOutput, StringComparison.Ordinal);

            // One setting changed, which moves the transition rollup and everything downstream of it.
            var changed = Path.Combine(outDir, "changed-config.yaml");
            File.WriteAllText(
                changed, File.ReadAllText(Config).Replace("min_transitions: 1", "min_transitions: 2"));

            var (changedCode, changedOutput) =
                Invoke("run", "-i", Input1, Input2, "-o", outDir, "-c", changed);
            Assert.Equal(0, changedCode);
            Assert.Contains("WARNING:", changedOutput, StringComparison.Ordinal);
            Assert.Contains("already holds results", changedOutput, StringComparison.Ordinal);
            Assert.Contains("corrected_peptides.parquet", changedOutput, StringComparison.Ordinal);
        }
        finally
        {
            Cleanup(outDir);
        }
    }
}
