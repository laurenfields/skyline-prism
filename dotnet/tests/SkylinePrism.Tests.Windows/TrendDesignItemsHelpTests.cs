using System;
using System.IO;
using Xunit;

namespace SkylinePrism.Tests.Windows;

/// <summary>
/// The two Trend design entries are GRAYED, not hidden, when there is no axis to fit a trend against,
/// and their tooltip is what says why. The same two rules as <see cref="DisabledControlHelpTests"/>:
/// the tooltip must survive being disabled, and it must not be fixed in the markup, because the code
/// writes the reason into it.
/// </summary>
public class TrendDesignItemsHelpTests
{
    private static string Xaml => File.ReadAllText(Path.GetFullPath(Path.Combine(
        AppContext.BaseDirectory, "..", "..", "..", "..", "..",
        "src", "SkylinePrism.App", "MainWindow.xaml")));

    [Theory]
    [InlineData("DiffDesignTrendItem")]
    [InlineData("DiffDesignTrendSubjectItem")]
    public void ATrendEntry_ShowsItsTooltipWhileDisabled_AndTakesItFromCode(string name)
    {
        var element = ElementMarkup(name);

        Assert.Contains("ToolTipService.ShowOnDisabled=\"True\"", element, StringComparison.Ordinal);
        Assert.DoesNotContain("ToolTip=\"", element, StringComparison.Ordinal);
    }

    /// <summary>The markup of one named element's opening tag.</summary>
    private static string ElementMarkup(string name)
    {
        var xaml = Xaml;
        var start = xaml.IndexOf($"x:Name=\"{name}\"", StringComparison.Ordinal);
        Assert.True(start >= 0, $"{name} was not found in MainWindow.xaml - was it renamed?");
        var open = xaml.LastIndexOf('<', start);
        var end = xaml.IndexOf('>', start);
        Assert.True(open >= 0 && end > open, $"could not delimit the {name} element");
        return xaml[open..end];
    }
}
