using System;
using System.Linq;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Tests.TestSupport;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// What a feature is known by - the hover readout on the Volcano, and the names used to find it in a
/// Skyline document. The rule that matters throughout is that a missing identity column yields
/// NOTHING rather than a placeholder: a caller looking an element up by name must never be handed a
/// protein-group id or a peptide sequence to search for.
/// </summary>
public class FeatureIdentityTests
{
    private static string MiniOutput => Fixtures.Path2("mini", "e2e-sum", "output");

    private static FeatureIdentity Identity(
        string id, string label, string[]? groups = null, string[]? accessions = null,
        string[]? names = null, string[]? genes = null) =>
        new(id, label,
            groups ?? Array.Empty<string>(),
            accessions ?? Array.Empty<string>(),
            names ?? Array.Empty<string>(),
            genes ?? Array.Empty<string>());

    [Fact]
    public void Describe_AddsTheGeneAndProteinBehindTheLabel()
    {
        var d = Identity("PG0001", "ALB",
            groups: new[] { "PG0001" },
            names: new[] { "sp|P02768|ALBU_HUMAN" },
            genes: new[] { "ALB" });

        // The gene is the label here, so it is not repeated; the protein name still earns its place.
        Assert.Equal("ALB  (sp|P02768|ALBU_HUMAN)", d.Describe());
    }

    /// <summary>
    /// A peptide's label is its sequence, which says nothing about what protein it came from - the
    /// whole reason hovering a Volcano point needs this.
    /// </summary>
    [Fact]
    public void Describe_APeptideCarriesItsGeneAndProtein()
    {
        var d = Identity("LVNELTEFAK", "LVNELTEFAK",
            groups: new[] { "PG0001" },
            names: new[] { "sp|P02768|ALBU_HUMAN" },
            genes: new[] { "ALB" });

        Assert.Equal("LVNELTEFAK  (ALB - sp|P02768|ALBU_HUMAN)", d.Describe());
    }

    [Fact]
    public void Describe_SaysWhenAPeptideIsSharedAcrossGroups()
    {
        var d = Identity("SHAREDPEPK", "SHAREDPEPK",
            groups: new[] { "PG0002", "PG0007" },
            names: new[] { "sp|P68871|HBB_HUMAN", "sp|Q9Y6K9|NEMO_HUMAN" },
            genes: new[] { "HBB", "IKBKG" });

        var text = d.Describe();

        Assert.True(d.IsShared);
        Assert.Contains("shared across 2 protein groups", text);
        Assert.StartsWith("SHAREDPEPK", text);
    }

    /// <summary>
    /// With no identity columns there is nothing to add, and nothing is invented - the label stands
    /// alone rather than being decorated with the id it already is.
    /// </summary>
    [Fact]
    public void Describe_WithNoIdentityColumns_IsJustTheLabel()
    {
        var d = Identity("PEPTIDEK", "PEPTIDEK");

        Assert.Equal("PEPTIDEK", d.Describe());
        Assert.False(d.IsShared);
    }

    [Fact]
    public void Describe_FallsBackToTheAccessionWhenThereIsNoProteinName()
    {
        var d = Identity("PG0003", "PEPTIDEK",
            groups: new[] { "PG0003" },
            accessions: new[] { "P12345" });

        Assert.Equal("PEPTIDEK  (P12345)", d.Describe());
    }

    /// <summary>
    /// The protein matrix carries the full identity, so a protein feature resolves to a real gene
    /// and protein name off the committed fixture.
    /// </summary>
    [Fact]
    public void IdentityOf_Protein_ReadsTheGroupColumnsFromTheMatrix()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);

        var d = ds.IdentityOf(0);

        Assert.Equal(ds.FeatureIds[0], d.FeatureId);
        Assert.NotEmpty(d.Genes);
        Assert.NotEmpty(d.ProteinNames);
        Assert.All(d.Genes, g => Assert.False(string.IsNullOrEmpty(g)));
    }

    /// <summary>
    /// This fixture predates the group columns on corrected_peptides, so every list must come back
    /// empty - never seeded with the sequence, which is the failure this whole type guards against.
    /// </summary>
    [Fact]
    public void IdentityOf_Peptide_OnAMatrixWithoutGroupColumns_HasNoNames()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Peptide);

        var d = ds.IdentityOf(0);

        Assert.Equal(ds.FeatureIds[0], d.FeatureId);
        Assert.Empty(d.ProteinGroups);
        Assert.Empty(d.Accessions);
        Assert.Empty(d.ProteinNames);
        Assert.Empty(d.Genes);
        // And Describe still has something useful to say.
        Assert.Equal(d.Label, d.Describe());
    }

    [Fact]
    public void IdentityOf_AnUnknownFeatureId_IsNull()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);

        Assert.Null(ds.IdentityOf("not-a-feature"));
        Assert.NotNull(ds.IdentityOf(ds.FeatureIds[0]));
    }

    /// <summary>
    /// A shared peptide's four group columns are ";"-separated and index-aligned, so they must be
    /// split apart rather than used whole - "PG0002;PG0007" matches no node in a document, and the
    /// point of keeping all of them is that a group whose leading protein is missing can still be
    /// found by another member.
    /// </summary>
    [Fact]
    public void IdentityOf_SplitsSemicolonSeparatedGroupsAndKeepsThemAligned()
    {
        var ds = DifferentialDataset.Load(MiniOutput, FeatureLevel.Protein);

        for (var i = 0; i < ds.FeatureIds.Length; i++)
        {
            var d = ds.IdentityOf(i);
            Assert.DoesNotContain(d.ProteinGroups, g => g.Contains(';', StringComparison.Ordinal));
            Assert.DoesNotContain(d.Genes, g => g.Contains(';', StringComparison.Ordinal));
            // Index-aligned: a name per group, where names are carried at all.
            if (d.ProteinNames.Count > 0)
                Assert.Equal(d.ProteinGroups.Count, d.ProteinNames.Count);
        }
    }
}
