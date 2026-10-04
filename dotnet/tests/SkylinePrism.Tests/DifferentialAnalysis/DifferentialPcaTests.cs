using System;
using System.Linq;
using SkylinePrism.Core.Numerics;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Parity tests for the differential explorer's PCA configuration against the Python explorer's
/// <c>compute_pca</c> (numpy SVD). Variance ratios are compared exactly; scores are compared by
/// absolute value because the sign of a principal component is arbitrary (a component and its
/// negation are equivalent).
///
/// <para>These ran against a separate <c>DifferentialPca</c> until the two PCA implementations were
/// merged into <see cref="Pca"/>. The expected values are unchanged, which is the point of keeping
/// them: the merge had to leave every number where it was, and a rewritten expectation would not
/// have shown that.</para>
/// </summary>
public class DifferentialPcaTests
{
    private const int Fn = 10;
    private const int Sn = 6;

    private static double[,] BuildMatrix()
    {
        var m = new double[Fn, Sn];
        for (var f = 0; f < Fn; f++)
        for (var s = 0; s < Sn; s++)
        {
            var v = 10.0 + ((f * 2 + s) % 5) / 2.0;
            if (s >= 3)
                v += f % 3;
            m[f, s] = v;
        }

        return m;
    }

    private static readonly string[] Ids = Enumerable.Range(0, Sn).Select(s => $"s{s}").ToArray();
    private static readonly int[] AllSamples = Enumerable.Range(0, Sn).ToArray();

    /// <summary>The differential pane's PCA settings: centered but unscaled, complete-case, 6
    /// components, and failing loudly rather than returning zeros.</summary>
    private static PcaOptions Options(string[] ids, int[] samples) => new()
    {
        Components = 6,
        Scaling = PcaScaling.CenterOnly,
        Missing = PcaMissingPolicy.CompleteCase,
        SampleColumns = samples,
        SampleIds = ids,
        RequireSufficientData = true,
    };

    [Fact]
    public void Compute_MatchesExplorer()
    {
        var res = Pca.Fit(BuildMatrix(), Options(Ids, AllSamples));

        Assert.Equal(10, res.NFeaturesUsed);
        Assert.Equal(new[] { "s0", "s1", "s2", "s3", "s4", "s5" }, res.SampleIds);

        var vr = new[]
        {
            0.6159451748402038, 0.21865889212828, 0.07288629737609328,
            0.06336746036742594, 0.02914217528799707,
        };
        for (var j = 0; j < vr.Length; j++)
            Assert.Equal(vr[j], res.VarianceRatio[j], 9);
        Assert.True(Math.Abs(res.VarianceRatio[5]) < 1e-9); // degenerate last component

        // |PC1| per sample (Python signs: s0..s2 positive, s3..s5 negative).
        var pc1 = new[]
        {
            1.6407831333439498, 3.0717601340573677, 2.402011336841677,
            2.592986409380687, 2.689809988979348, 1.831758205882959,
        };
        for (var i = 0; i < Sn; i++)
            Assert.Equal(pc1[i], Math.Abs(res.Scores[i, 0]), 9);

        // A couple of |PC2| values.
        Assert.Equal(1.7677669529663693, Math.Abs(res.Scores[0, 1]), 9);
        Assert.Equal(1.7677669529663704, Math.Abs(res.Scores[2, 1]), 9);
    }

    [Fact]
    public void Compute_MoreSamplesThanFeatures_LimitsComponentCount()
    {
        // 3 features x 6 samples: numpy's thin SVD yields min(nSamples, nFeatures) = 3 components,
        // not 6 padded with zero-rank eigenpairs.
        var m = new double[3, 6];
        for (var f = 0; f < 3; f++)
        for (var s = 0; s < 6; s++)
            m[f, s] = 10.0 + ((f * 2 + s) % 5) / 2.0;

        var ids = Enumerable.Range(0, 6).Select(s => $"s{s}").ToArray();
        var res = Pca.Fit(m, Options(ids, Enumerable.Range(0, 6).ToArray()));

        Assert.Equal(3, res.NFeaturesUsed);
        Assert.Equal(3, res.VarianceRatio.Length);
        Assert.Equal(3, res.Scores.GetLength(1));
    }

    [Fact]
    public void Compute_RejectsTooFewSamples()
    {
        Assert.Throws<ArgumentException>(() =>
            Pca.Fit(BuildMatrix(), Options(Ids, new[] { 0 })));
    }

    [Fact]
    public void Compute_RejectsTooFewCompleteFeatures()
    {
        // Knock out all but one feature by putting a NaN in a selected column of each.
        var m = BuildMatrix();
        for (var f = 1; f < Fn; f++)
            m[f, 0] = double.NaN;
        Assert.Throws<ArgumentException>(() => Pca.Fit(m, Options(Ids, AllSamples)));
    }
}
