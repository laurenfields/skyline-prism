using System;
using System.IO;
using System.Text.RegularExpressions;
using Xunit;

namespace SkylinePrism.Tests.Windows;

/// <summary>
/// A control this window GREYS rather than hides has to be able to say why.
/// </summary>
/// <remarks>
/// WPF suppresses tooltips on disabled elements unless <c>ToolTipService.ShowOnDisabled</c> is set,
/// so a greyed control's explanation is invisible at exactly the moment it is the only thing the
/// reader needs. That is a one-attribute mistake to make and an invisible one to find - a reviewer
/// sees a tooltip in the markup and assumes it shows - so it is pinned here rather than trusted.
/// </remarks>
public class DisabledControlHelpTests
{
    private static string Xaml => File.ReadAllText(Path.GetFullPath(Path.Combine(
        AppContext.BaseDirectory, "..", "..", "..", "..", "..",
        "src", "SkylinePrism.App", "MainWindow.xaml")));

    /// <summary>
    /// The prior-source checkbox is greyed when the run has no control replicates - deliberately,
    /// rather than hidden, so a reader learns the option exists and that setting sample types in
    /// Skyline would enable it. That only works if the tooltip survives being disabled.
    /// </summary>
    [Fact]
    public void ThePriorSourceCheckbox_ShowsItsTooltipWhileDisabled()
    {
        var element = ElementMarkup("DiffPriorFromControlsCheck");

        Assert.Contains("ToolTipService.ShowOnDisabled=\"True\"", element, StringComparison.Ordinal);
    }

    /// <summary>
    /// Its help text is set in code, not markup, because the disabled case prefixes it with the
    /// reason. A static ToolTip in the XAML would silently win over that.
    /// </summary>
    [Fact]
    public void ThePriorSourceCheckbox_HasNoStaticTooltipToOverrideTheDynamicOne()
    {
        var element = ElementMarkup("DiffPriorFromControlsCheck");

        Assert.DoesNotContain("ToolTip=\"", element, StringComparison.Ordinal);
    }

    /// <summary>The markup of one named element, from its opening tag to the end of that tag.</summary>
    private static string ElementMarkup(string name)
    {
        var xaml = Xaml;
        var start = xaml.IndexOf($"x:Name=\"{name}\"", StringComparison.Ordinal);
        Assert.True(start >= 0, $"{name} was not found in MainWindow.xaml - was it renamed?");

        // Back up to the opening '<' of this element, then forward to the tag's close.
        var open = xaml.LastIndexOf('<', start);
        var end = xaml.IndexOf('>', start);
        Assert.True(open >= 0 && end > open, $"could not delimit the {name} element");
        return xaml[open..end];
    }
}
