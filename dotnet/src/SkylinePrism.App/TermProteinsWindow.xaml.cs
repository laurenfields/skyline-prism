using System.Collections.Generic;
using System.Globalization;
using System.Windows;

namespace SkylinePrism.App;

/// <summary>
/// A per-term detail popup: the significant proteins an enrichment term contains, opened by clicking a
/// term in the Enrichment table. Reused across clicks (one window, redrawn).
/// </summary>
public partial class TermProteinsWindow : Window
{
    /// <summary>One protein row shown for a term.</summary>
    public sealed record TermProteinRow(string Protein, string Gene, double Log2FC, double AdjP);

    public TermProteinsWindow()
    {
        InitializeComponent();
    }

    /// <summary>Show (or refresh) the proteins for one term.</summary>
    public void ShowTerm(string title, int count, double pValue, IReadOnlyList<TermProteinRow> proteins)
    {
        var inv = CultureInfo.InvariantCulture;
        HeaderText.Text = title;
        SubText.Text = $"{count} significant protein(s) in this term   -   term p = {pValue.ToString("0.##e0", inv)}";
        Grid.ItemsSource = proteins;
    }
}
