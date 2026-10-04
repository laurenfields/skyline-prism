using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Which subjects a paired analysis actually runs over, and whether it says so accurately.
/// </summary>
/// <remarks>
/// The counts in these messages are the only place a reader learns that part of their cohort took no
/// part in the contrast, so an inflated or missing count is a quietly wrong statement about the
/// result rather than a cosmetic defect.
/// </remarks>
public class PairedSamplesTests
{
    private static string?[] Subjects(params string?[] labels) => labels;

    [Fact]
    public void MatchesEachSubjectsTwoSamples()
    {
        // Columns 0..2 are arm A, 3..5 arm B, deliberately in a different subject order.
        var subjects = Subjects("s1", "s2", "s3", "s3", "s1", "s2");

        var (pairs, messages) = PairedSamples.Resolve(subjects, new[] { 0, 1, 2 }, new[] { 3, 4, 5 });

        Assert.Equal(3, pairs.Count);
        Assert.Empty(messages);
        // Matched by label, not by position.
        Assert.Equal(new[] { (0, 4), (1, 5), (2, 3) },
            pairs.Select(p => (p.AColumn, p.BColumn)).ToArray());
    }

    /// <summary>
    /// A subject with two samples in one arm is counted ONCE, not once per column it occupies. The
    /// resolution loop runs over columns, so the obvious guard - "have I already paired this
    /// subject" - only skips the ones that succeeded, and reported a single ambiguous subject as two.
    /// </summary>
    [Fact]
    public void AnAmbiguousSubjectIsReportedOnce()
    {
        // s1 appears twice in arm A: which of its samples pairs with B's is not in the metadata.
        var subjects = Subjects("s1", "s1", "s2", "s1", "s2");

        var (pairs, messages) = PairedSamples.Resolve(subjects, new[] { 0, 1, 2 }, new[] { 3, 4 });

        Assert.Single(pairs);
        Assert.Equal("s2", pairs[0].Subject);

        var ambiguous = Assert.Single(messages, m => m.Contains("ambiguous"));
        Assert.Contains("1 subject(s)", ambiguous);
        Assert.DoesNotContain("2 subject(s)", ambiguous);
        // And the preview names it once.
        Assert.Equal(1, ambiguous.Split("s1").Length - 1);
    }

    /// <summary>A subject in one arm only is counted once too, from whichever side it is missing.</summary>
    [Fact]
    public void AnUnmatchedSubjectIsReportedOnceFromEitherSide()
    {
        // aOnly is twice in A and absent from B; bOnly is only in B.
        var subjects = Subjects("both", "aOnly", "aOnly", "both", "bOnly");

        var (pairs, messages) = PairedSamples.Resolve(subjects, new[] { 0, 1, 2 }, new[] { 3, 4 });

        Assert.Single(pairs);
        var unmatched = Assert.Single(messages, m => m.Contains("only one arm"));
        Assert.Contains("2 subject(s)", unmatched); // aOnly and bOnly, each once
        Assert.Equal(1, unmatched.Split("aOnly").Length - 1);
    }

    /// <summary>
    /// A sample with no value in the pairing column is counted separately, because the cause is a gap
    /// in the metadata rather than an unbalanced design - and the fix is different.
    /// </summary>
    [Fact]
    public void SamplesWithNoSubjectValueAreCountedSeparately()
    {
        var subjects = Subjects("s1", null, "", "s1", "s2");

        var (pairs, messages) = PairedSamples.Resolve(subjects, new[] { 0, 1, 2 }, new[] { 3, 4 });

        Assert.Single(pairs);
        var gap = Assert.Single(messages, m => m.Contains("no value in the pairing column"));
        Assert.Contains("2 sample(s)", gap);
    }

    [Fact]
    public void NothingMatches_IsEmptyRatherThanAnException()
    {
        var subjects = Subjects("a", "b");

        var (pairs, messages) = PairedSamples.Resolve(subjects, new[] { 0 }, new[] { 1 });

        Assert.Empty(pairs);
        Assert.NotEmpty(messages);
    }
}
