using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// The name-and-count preview shared by the pairing messages and the Markers "not detected" note.
/// </summary>
public class NamePreviewTests
{
    [Fact]
    public void FewerNamesThanTheCap_AreAllShown_WithNoSuffix()
    {
        // The common case must read naturally, and must never say "+0 more".
        Assert.Equal("A, B", NamePreview.Of(new[] { "A", "B" }, show: 3));
        Assert.Equal("A, B, C", NamePreview.Of(new[] { "A", "B", "C" }, show: 3));
    }

    [Fact]
    public void MoreNamesThanTheCap_AreTruncated_AndTheRestCounted()
        => Assert.Equal("A, B, C, +2 more", NamePreview.Of(new[] { "A", "B", "C", "D", "E" }, show: 3));

    [Fact]
    public void TheCountIsOfWhatWasLeftOut_NotOfTheWholeList()
    {
        // 1,273 not-detected members with a cap of 12 leaves 1,261 - the number that made this
        // helper necessary.
        var names = Enumerable.Range(0, 1273).Select(i => $"G{i}").ToArray();

        Assert.EndsWith(", +1261 more", NamePreview.Of(names, show: 12));
    }

    [Fact]
    public void AnEmptyList_IsAnEmptyString()
        => Assert.Equal("", NamePreview.Of(Array.Empty<string>(), show: 3));

    [Fact]
    public void ANonPositiveCap_IsRefused()
        => Assert.Throws<ArgumentOutOfRangeException>(() => NamePreview.Of(new[] { "A" }, show: 0));
}
