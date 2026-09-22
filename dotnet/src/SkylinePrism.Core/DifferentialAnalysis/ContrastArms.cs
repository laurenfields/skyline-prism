using System;
using System.Collections.Generic;
using System.Linq;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// The sample columns of a two-arm contrast, or the reason there is no contrast to run.
/// </summary>
/// <remarks>
/// <see cref="Error"/> is written for a user, because both callers show it: the Differential pane in
/// its status line, the CLI on stderr before exiting non-zero.
/// </remarks>
public sealed record ArmColumns(IReadOnlyList<int> A, IReadOnlyList<int> B, string? Error)
{
    /// <summary>Whether a contrast was resolved. False means <see cref="Error"/> says why not.</summary>
    public bool Ok => Error is null;
}

/// <summary>
/// Turns a metadata column plus the levels picked for each arm into the two column-index lists a
/// contrast runs over. Each arm is a UNION of levels, so "A = Control + Mild" against "B = Severe"
/// is one contrast rather than three.
/// </summary>
/// <remarks>
/// In Core, and not in the window that grew it, because it decides which samples land on which side
/// of a test: the GUI and <c>prism differential</c> must resolve the arms identically or the same
/// selection gives two different answers, with nothing to show which one was asked for.
/// </remarks>
public static class ContrastArms
{
    /// <summary>
    /// Resolve the arms. <paramref name="values"/> is the group-by column's value per sample column
    /// (null where the sample has none), matching the matrix's column order.
    /// </summary>
    public static ArmColumns Resolve(
        string?[] values, IEnumerable<string> aLevels, IEnumerable<string> bLevels)
    {
        var aSet = aLevels.ToList();
        var bSet = bLevels.ToList();
        if (aSet.Count == 0 || bSet.Count == 0)
            return new ArmColumns(Array.Empty<int>(), Array.Empty<int>(),
                "Pick at least one value for each arm.");

        // A level in both arms would put the same samples on both sides, which is not a contrast.
        // Refused rather than dropped from one side, because WHICH side it was dropped from would
        // change the answer, and nothing in the output would record the choice.
        var both = aSet.Intersect(bSet, StringComparer.Ordinal).ToList();
        if (both.Count > 0)
            return new ArmColumns(Array.Empty<int>(), Array.Empty<int>(),
                $"'{string.Join("', '", both)}' is in both arms; a value can only be on one side.");

        var inA = new HashSet<string>(aSet, StringComparer.Ordinal);
        var inB = new HashSet<string>(bSet, StringComparer.Ordinal);
        var groupA = new List<int>();
        var groupB = new List<int>();
        for (var j = 0; j < values.Length; j++)
        {
            if (values[j] is not { } v)
                continue;
            if (inA.Contains(v))
                groupA.Add(j);
            else if (inB.Contains(v))
                groupB.Add(j);
        }

        // Named separately: "B matched nothing" is a typo in one level, and saying which arm is
        // empty is the difference between fixing it and re-reading the whole command line.
        if (groupA.Count == 0 || groupB.Count == 0)
        {
            var empty = groupA.Count == 0 && groupB.Count == 0
                ? "Neither arm matched any sample."
                : groupA.Count == 0 ? "Arm A matched no samples." : "Arm B matched no samples.";
            return new ArmColumns(Array.Empty<int>(), Array.Empty<int>(), empty);
        }

        return new ArmColumns(groupA, groupB, null);
    }

    /// <summary>How the arms read in a status line or a results header ("Control + Mild").</summary>
    public static string Describe(IEnumerable<string> levels) => string.Join(" + ", levels);
}

/// <summary>
/// The sample-type values that mean "this is a control injection", in both vocabularies: the
/// Replicates report carries Skyline's spellings, while PRISM's synthetic Sample Type column uses
/// its own mapped names.
/// </summary>
public static class ControlSampleTypes
{
    /// <summary>The control spellings, across both sources.</summary>
    public static readonly string[] Values = { "Standard", "Quality Control", "QC", "reference", "qc" };

    /// <summary>Comparer used throughout: annotation spellings vary in case between sources.</summary>
    public static StringComparer Comparer => StringComparer.OrdinalIgnoreCase;

    /// <summary>Whether a sample-type value denotes a control sample.</summary>
    public static bool IsControl(string? value) => value is not null && Values.Contains(value, Comparer);

    /// <summary>
    /// The control columns, grouped by type, for fitting the variance prior on controls rather than
    /// on the contrast groups. Null when there are too few to fit anything.
    /// </summary>
    /// <remarks>
    /// One group per control TYPE, not one pooled set: the variance is computed within a group, so
    /// pooling QC and reference - different materials, injected at different amounts - would count
    /// the systematic gap between them as measurement noise.
    ///
    /// <para>Two is the minimum a variance can be computed from at all, so fewer than that is not a
    /// thin prior but no prior - the caller turns the option off rather than falling back
    /// silently.</para>
    /// </remarks>
    public static IReadOnlyList<IReadOnlyList<int>>? PriorGroups(string?[] sampleTypes)
    {
        var groups = Enumerable.Range(0, sampleTypes.Length)
            .Where(i => IsControl(sampleTypes[i]))
            // Grouped by the SAME comparer that selected them. Ordinal here while IsControl is
            // case-insensitive splits one control type across its spellings: a cohort merged from
            // two documents writing "QC" and "qc" - both in Values - yields two groups of one,
            // the two-replicate floor drops both, and the prior silently falls back to the design
            // groups. That is the most consequential default in the pane, reverting without a word.
            .GroupBy(i => sampleTypes[i], Comparer)
            .Select(g => (IReadOnlyList<int>)g.ToList())
            .Where(g => g.Count >= 2)
            .ToList();

        return groups.Count > 0 ? groups : null;
    }
}
