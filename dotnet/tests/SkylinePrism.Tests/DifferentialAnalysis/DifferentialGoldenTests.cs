using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;
using SkylinePrism.Core.DifferentialAnalysis;
using SkylinePrism.Core.DifferentialAnalysis.Detection;
using SkylinePrism.Core.Numerics;
using SkylinePrism.Tests.TestSupport;
using Xunit;

namespace SkylinePrism.Tests.DifferentialAnalysis;

/// <summary>
/// Reads one golden fixture from <c>dotnet/tests/fixtures/differential/</c>.
///
/// <para>Every floating-point value in those files is a JSON <b>string</b> holding Python's
/// shortest round-trip repr, not a JSON number. That is deliberate: bare <c>NaN</c>/<c>Infinity</c>
/// are not valid JSON and <c>System.Text.Json</c> rejects them, and a value routed through both
/// encoders' number paths is only approximately preserved - which is no basis for an assertion at
/// 1e-15. Through a string, the 64 bits Python wrote are the 64 bits this reads.</para>
/// </summary>
internal static class Golden
{
    public static JsonElement Load(string file)
    {
        var path = Path.Combine(Fixtures.Root, "differential", file);
        Assert.True(File.Exists(path),
            $"Missing golden fixture {path}. Regenerate with "
            + "'uv run dotnet/tests/fixtures/differential/generate.py' from the repository root.");
        using var doc = JsonDocument.Parse(File.ReadAllText(path));
        return doc.RootElement.Clone();
    }

    public static IEnumerable<JsonElement> Cases(string file)
        => Load(file).GetProperty("cases").EnumerateArray();

    /// <summary>All case names in a fixture, as xUnit theory data.</summary>
    public static IEnumerable<object[]> CaseNames(string file)
        => Cases(file).Select(c => new object[] { c.GetProperty("name").GetString()! });

    public static JsonElement Case(string file, string name)
        => Cases(file).Single(c => c.GetProperty("name").GetString() == name);

    /// <summary>
    /// One golden value. Python's repr spells the non-finite doubles <c>inf</c>, <c>-inf</c> and
    /// <c>nan</c>, none of which .NET's parser accepts, so they are mapped rather than being
    /// rewritten in the fixture - a golden should read the way the tool that wrote it writes.
    /// </summary>
    public static double Num(JsonElement e)
    {
        var s = e.GetString()!;
        return s switch
        {
            "inf" => double.PositiveInfinity,
            "-inf" => double.NegativeInfinity,
            "nan" => double.NaN,
            _ => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture),
        };
    }

    public static double Num(JsonElement parent, string prop) => Num(parent.GetProperty(prop));

    public static double[] Vec(JsonElement e)
        => e.EnumerateArray().Select(Num).ToArray();

    public static double[] Vec(JsonElement parent, string prop) => Vec(parent.GetProperty(prop));

    public static double[,] Mat(JsonElement e)
    {
        var rows = e.EnumerateArray().Select(r => Vec(r)).ToArray();
        if (rows.Length == 0)
            return new double[0, 0];
        var m = new double[rows.Length, rows[0].Length];
        for (var i = 0; i < rows.Length; i++)
            for (var j = 0; j < rows[0].Length; j++)
                m[i, j] = rows[i][j];
        return m;
    }

    public static double[,] Mat(JsonElement parent, string prop) => Mat(parent.GetProperty(prop));

    /// <summary>
    /// Relative-error assertion that treats an exact-zero expectation as an absolute one, and
    /// requires the two sides to agree on NaN and on infinities rather than silently passing when
    /// both are non-finite in different ways.
    /// </summary>
    public static void Close(double expected, double actual, double rtol, string what)
    {
        if (double.IsNaN(expected) || double.IsNaN(actual))
        {
            Assert.True(double.IsNaN(expected) && double.IsNaN(actual),
                $"{what}: expected {expected:R}, actual {actual:R} (NaN mismatch)");
            return;
        }

        if (double.IsInfinity(expected) || double.IsInfinity(actual))
        {
            Assert.True(expected.Equals(actual),
                $"{what}: expected {expected:R}, actual {actual:R} (infinity mismatch)");
            return;
        }

        var denom = Math.Abs(expected);
        var err = denom > 0 ? Math.Abs(actual - expected) / denom : Math.Abs(actual - expected);
        Assert.True(err <= rtol,
            $"{what}: expected {expected:R}, actual {actual:R}, rel err {err:R} > {rtol:R}");
    }

    public static void CloseAll(double[] expected, double[] actual, double rtol, string what)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (var i = 0; i < expected.Length; i++)
            Close(expected[i], actual[i], rtol, $"{what}[{i}]");
    }
}

/// <summary>
/// Cross-implementation parity for <see cref="SkylinePrism.Core.DifferentialAnalysis"/> against
/// committed goldens generated from scipy, statsmodels and inmoose
/// (<c>dotnet/tests/fixtures/differential/generate.py</c>).
///
/// <para>This is a different kind of test from the hand-written parity tests beside it. Those carry
/// a handful of reference values transcribed into the source; these read whole vectors and matrices
/// from a file that a script regenerates, so coverage is a question of editing the generator rather
/// than of pasting more constants, and the provenance of every number is the script rather than a
/// comment claiming what produced it.</para>
///
/// <para>The tolerances are not uniform, and the differences are the point:</para>
/// <list type="bullet">
/// <item><b>1e-15 (bit-exact for practical purposes)</b> for closed-form arithmetic on the same
/// inputs - BH, Fisher, polygamma, OLS.</item>
/// <item><b>1e-9</b> where an iteration converges to a tolerance - squeezeVar's Newton solve for the
/// prior df, and everything downstream of it.</item>
/// <item><b>1e-4 on Firth coefficients</b>, because the reference is a derivative-free optimizer:
/// it pins the log-likelihood far more tightly than the coefficients, and on a flat penalized
/// likelihood the coefficients are genuinely not determined to more than that. The log-likelihood
/// is asserted at 1e-9 in the same test, which is where the evidence actually is.</item>
/// </list>
/// </summary>
public class DifferentialGoldenTests
{
    // ------------------------------------------------------------------------------------
    // Benjamini-Hochberg  vs  statsmodels multipletests(method="fdr_bh")
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> FdrCases() => Golden.CaseNames("fdr.json");

    [Theory]
    [MemberData(nameof(FdrCases))]
    public void BenjaminiHochberg_MatchesStatsmodels(string name)
    {
        var c = Golden.Case("fdr.json", name);
        var actual = Fdr.BenjaminiHochberg(Golden.Vec(c, "p"));
        Golden.CloseAll(Golden.Vec(c, "expected"), actual, 1e-15, $"fdr/{name}");
    }

    // ------------------------------------------------------------------------------------
    // Polygamma  vs  scipy.special.polygamma
    // ------------------------------------------------------------------------------------

    [Fact]
    public void Trigamma_And_Tetragamma_MatchScipy()
    {
        var g = Golden.Load("polygamma.json");
        var xs = Golden.Vec(g, "x");
        var tri = Golden.Vec(g, "trigamma");
        var tet = Golden.Vec(g, "tetragamma");

        for (var i = 0; i < xs.Length; i++)
        {
            Golden.Close(tri[i], EmpiricalBayes.Trigamma(xs[i]), 1e-13, $"trigamma({xs[i]:R})");
            Golden.Close(tet[i], EmpiricalBayes.Tetragamma(xs[i]), 1e-13, $"tetragamma({xs[i]:R})");
        }
    }

    [Fact]
    public void TrigammaInverse_MatchesScipyInvertedNumerically()
    {
        var g = Golden.Load("polygamma.json");
        var xs = Golden.Vec(g, "trigamma_inverse_x");
        var expected = Golden.Vec(g, "trigamma_inverse_expected");

        // 1e-7 rather than 1e-13: limma's trigammaInverse stops when the Newton step is small
        // relative to y (-dif/y < 1e-8), not when the residual is, so the two agree to about that
        // and no further. The quantity that matters downstream - the prior df - inherits exactly
        // this, which is why squeezeVar below is asserted at 1e-9 and not tighter.
        for (var i = 0; i < xs.Length; i++)
            Golden.Close(expected[i], EmpiricalBayes.TrigammaInverse(xs[i]), 1e-7,
                $"trigammaInverse({xs[i]:R})");
    }

    // ------------------------------------------------------------------------------------
    // squeezeVar  vs  inmoose.limma.squeezeVar
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> SqueezeVarCases() => Golden.CaseNames("squeezevar.json");

    [Theory]
    [MemberData(nameof(SqueezeVarCases))]
    public void SqueezeVar_MatchesInmoose(string name)
    {
        var c = Golden.Case("squeezevar.json", name);
        var variances = Golden.Vec(c, "variances");
        var df = Golden.Num(c, "df_residual");
        var covariateEl = c.GetProperty("covariate");

        var res = covariateEl.ValueKind == JsonValueKind.Null
            ? EmpiricalBayes.SqueezeVarGlobal(variances, df)
            : EmpiricalBayes.SqueezeVarTrend(variances, df, Golden.Vec(covariateEl));

        Golden.Close(Golden.Num(c, "df_prior"), res.DfPrior, 1e-9, $"{name}/df_prior");
        Golden.CloseAll(Golden.Vec(c, "var_prior"), res.VarPrior, 1e-9, $"{name}/var_prior");
        Golden.CloseAll(Golden.Vec(c, "var_post"), res.VarPost, 1e-9, $"{name}/var_post");
    }

    // ------------------------------------------------------------------------------------
    // Natural spline basis  vs  inmoose.utils.splines.ns
    //
    // Compared through the orthogonal projector onto the column span, because the basis itself is
    // fixed only up to the rotation the constraint QR happens to produce - two correct
    // implementations can return different matrices with the same span, and the span is the whole
    // of what the trend fit uses.
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> SplineCases() => Golden.CaseNames("spline.json");

    [Theory]
    [MemberData(nameof(SplineCases))]
    public void NaturalSplineBasis_SpansTheSameSubspaceAsInmoose(string name)
    {
        var c = Golden.Case("spline.json", name);
        var x = Golden.Vec(c, "x");
        var df = c.GetProperty("df").GetInt32();
        var expectedHat = Golden.Mat(c, "hat");

        var basis = NaturalSplineBasis.Build(x, df, includeIntercept: true);
        Assert.Equal(c.GetProperty("n_basis_columns").GetInt32(), basis.GetLength(1));

        var q = DenseMatrix.OfArray(basis).QR(QRMethod.Thin).Q;
        var hat = q * q.Transpose();

        var n = x.Length;
        for (var i = 0; i < n; i++)
            for (var j = 0; j < n; j++)
                // Absolute, not relative: a projector entry is legitimately near zero off the
                // diagonal, and a relative test there would be comparing rounding noise.
                Assert.True(Math.Abs(expectedHat[i, j] - hat[i, j]) < 1e-12,
                    $"spline/{name} hat[{i},{j}]: expected {expectedHat[i, j]:R}, got {hat[i, j]:R}");
    }

    // ------------------------------------------------------------------------------------
    // lmFit  vs  numpy.linalg.lstsq
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> LmFitCases() => Golden.CaseNames("lmfit.json");

    [Theory]
    [MemberData(nameof(LmFitCases))]
    public void LinearModelFit_MatchesNumpy(string name)
    {
        var c = Golden.Case("lmfit.json", name);
        var expr = Golden.Mat(c, "expr");
        var design = Golden.Mat(c, "design");

        var fit = LinearModel.Fit(expr, design);

        // PRISM solves through a thin QR; the reference forms inv(X'X) explicitly, which squares
        // the condition number. On a well-conditioned design that difference is invisible, and the
        // assertion is correspondingly tight. The ill_conditioned case exists precisely to make it
        // visible - its covariate almost tracks the group - and there the two legitimately part
        // company around 1e-11. Loosening only that case keeps the strong claim where it is
        // earned, and records that PRISM's path is the more accurate of the two rather than
        // hiding the gap under a blanket tolerance.
        var tol = name == "ill_conditioned" ? 1e-9 : 1e-12;

        Golden.Close(Golden.Num(c, "df_residual"), fit.DfResidual, 1e-15, $"{name}/df_residual");
        Golden.CloseAll(Golden.Vec(c, "sigma"), fit.Sigma, tol, $"{name}/sigma");
        Golden.CloseAll(Golden.Vec(c, "amean"), fit.Amean, 1e-14, $"{name}/amean");
        Golden.CloseAll(Golden.Vec(c, "stdev_unscaled"), fit.StdevUnscaled, tol,
            $"{name}/stdev_unscaled");

        var expectedCoef = Golden.Mat(c, "coefficients");
        for (var i = 0; i < expectedCoef.GetLength(0); i++)
            for (var j = 0; j < expectedCoef.GetLength(1); j++)
                Golden.Close(expectedCoef[i, j], fit.Coefficients[i, j], tol,
                    $"{name}/coef[{i},{j}]");
    }

    // ------------------------------------------------------------------------------------
    // End-to-end moderated t  vs  lstsq + squeezeVar + scipy t + statsmodels BH
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> ModeratedTCases() => Golden.CaseNames("moderated_t.json");

    [Theory]
    [MemberData(nameof(ModeratedTCases))]
    public void Differential_MatchesTheLimmaComposition(string name)
    {
        var c = Golden.Case("moderated_t.json", name);
        var expr = Golden.Mat(c, "expr");
        var nA = c.GetProperty("n_a").GetInt32();
        var nB = c.GetProperty("n_b").GetInt32();
        var trend = c.GetProperty("trend").GetBoolean();

        var nFeatures = expr.GetLength(0);
        var featureIds = Enumerable.Range(0, nFeatures).Select(i => $"F{i}").ToArray();
        var groupA = Enumerable.Range(0, nA).ToArray();
        var groupB = Enumerable.Range(nA, nB).ToArray();

        // Rebuild the covariates as PRISM takes them - raw and uncentered - so that the design
        // construction (centering, dummy coding, level ordering) is under test too, rather than
        // being handed the finished design columns.
        var covariates = new List<Covariate>();
        foreach (var cov in c.GetProperty("covariates").EnumerateArray())
        {
            var covName = cov.GetProperty("name").GetString()!;
            var kind = cov.GetProperty("kind").GetString();
            if (kind == "numeric")
                covariates.Add(new NumericCovariate(covName, Golden.Vec(cov, "values")));
            else
                covariates.Add(new CategoricalCovariate(covName,
                    cov.GetProperty("values").EnumerateArray().Select(v => v.GetString()).ToArray()));
        }

        var res = Differential.Run(expr, featureIds, groupA, groupB, minPerGroup: 2,
            covariates: covariates.Count > 0 ? covariates : null, trend: trend);

        Assert.Equal(nFeatures, res.NFeaturesTested);
        Golden.Close(Golden.Num(c, "df_residual"), res.DfResidual, 1e-15, $"{name}/df_residual");
        Golden.Close(Golden.Num(c, "df_prior"), res.DfPrior, 1e-9, $"{name}/df_prior");

        // Rows come back sorted by p-value, so map them home by feature id before comparing.
        var byId = res.Rows.ToDictionary(r => r.FeatureId, StringComparer.Ordinal);
        var logfc = Golden.Vec(c, "logfc");
        var amean = Golden.Vec(c, "amean");
        var t = Golden.Vec(c, "t");
        var p = Golden.Vec(c, "p");
        var adjP = Golden.Vec(c, "adj_p");

        for (var i = 0; i < nFeatures; i++)
        {
            var row = byId[$"F{i}"];
            Golden.Close(logfc[i], row.LogFc, 1e-11, $"{name}/F{i}/logFC");
            Golden.Close(amean[i], row.AveExpr, 1e-14, $"{name}/F{i}/AveExpr");
            Golden.Close(t[i], row.T, 1e-9, $"{name}/F{i}/t");
            Golden.Close(p[i], row.PValue, 1e-9, $"{name}/F{i}/p");
            Golden.Close(adjP[i], row.AdjPValue, 1e-9, $"{name}/F{i}/adj.P");
        }
    }

    // ------------------------------------------------------------------------------------
    // Fisher exact  vs  scipy.stats.fisher_exact
    // ------------------------------------------------------------------------------------

    [Fact]
    public void FisherExact_MatchesScipy()
    {
        foreach (var c in Golden.Cases("fisher.json"))
        {
            int a = c.GetProperty("a").GetInt32(), b = c.GetProperty("b").GetInt32();
            int d = c.GetProperty("c").GetInt32(), e = c.GetProperty("d").GetInt32();
            Golden.Close(Golden.Num(c, "p"), FisherExact.TwoSidedP(a, b, d, e), 1e-12,
                $"fisher({a},{b},{d},{e})");
        }
    }

    // ------------------------------------------------------------------------------------
    // Firth  vs  a derivative-free maximization of the same penalized log-likelihood
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> FirthCases() => Golden.CaseNames("firth.json");

    [Theory]
    [MemberData(nameof(FirthCases))]
    public void FirthLogit_ReachesTheSamePenalizedOptimum(string name)
    {
        var c = Golden.Case("firth.json", name);
        // tol far below the 1e-7 default on purpose. At the default the iteration stops when the
        // penalized log-likelihood moves by less than 1e-7, so comparing against a reference that
        // ran to convergence would be measuring that stopping rule rather than the estimator, and
        // would pass or fail on where the last accepted step happened to land.
        var res = FirthLogit.Fit(Golden.Mat(c, "x"), Golden.Vec(c, "y"), maxIter: 500, tol: 1e-14);

        Assert.True(res.Converged, $"firth/{name} did not converge");

        // The log-likelihood is what a direct-search optimum actually determines, so it carries the
        // assertion. Nelder-Mead and Powell agree with the Newton estimator to ~1e-9 here.
        Golden.Close(Golden.Num(c, "penalized_loglik"), res.PenalizedLogLik, 1e-9,
            $"firth/{name}/loglik");

        // The coefficients are a much weaker claim, and deliberately so: near a flat optimum the
        // argument moves far for a log-likelihood change of 1e-12. Asserting them tightly would be
        // asserting the optimizer's stopping point, not the estimate.
        var expectedBeta = Golden.Vec(c, "beta");
        Assert.Equal(expectedBeta.Length, res.Beta.Length);
        for (var i = 0; i < expectedBeta.Length; i++)
        {
            var err = Math.Abs(res.Beta[i] - expectedBeta[i])
                      / Math.Max(1.0, Math.Abs(expectedBeta[i]));
            Assert.True(err < 1e-4,
                $"firth/{name}/beta[{i}]: expected {expectedBeta[i]:R}, got {res.Beta[i]:R} "
                + $"(scaled err {err:R})");
        }
    }

    // ------------------------------------------------------------------------------------
    // Sample PCA  vs  numpy.linalg.svd
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> PcaCases() => Golden.CaseNames("pca.json");

    [Theory]
    [MemberData(nameof(PcaCases))]
    public void Pca_CenterOnlyCompleteCase_MatchesNumpySvd(string name)
    {
        var c = Golden.Case("pca.json", name);
        var expr = Golden.Mat(c, "expr");
        var nSamples = expr.GetLength(1);
        var sampleIds = Enumerable.Range(0, nSamples).Select(i => $"S{i}").ToArray();
        var nComponents = c.GetProperty("n_components").GetInt32();

        var res = Pca.Fit(expr, new PcaOptions
        {
            Components = nComponents,
            Scaling = PcaScaling.CenterOnly,
            Missing = PcaMissingPolicy.CompleteCase,
            SampleColumns = Enumerable.Range(0, nSamples).ToArray(),
            SampleIds = sampleIds,
            RequireSufficientData = true,
        });

        Assert.Equal(c.GetProperty("n_features_used").GetInt32(), res.NFeaturesUsed);
        Golden.CloseAll(Golden.Vec(c, "variance_ratio"), res.VarianceRatio, 1e-10,
            $"pca/{name}/variance_ratio");

        // Component signs are arbitrary (an eigenvector and its negation are the same component),
        // so the golden holds magnitudes and so does the comparison.
        var expectedAbs = Golden.Mat(c, "abs_scores");
        Assert.Equal(expectedAbs.GetLength(1), res.Scores.GetLength(1));
        for (var i = 0; i < expectedAbs.GetLength(0); i++)
            for (var j = 0; j < expectedAbs.GetLength(1); j++)
                Golden.Close(expectedAbs[i, j], Math.Abs(res.Scores[i, j]), 1e-9,
                    $"pca/{name}/|score[{i},{j}]|");
    }

    // ------------------------------------------------------------------------------------
    // Covariate-adjusted detection  vs  the Firth-penalized LRT computed independently
    // ------------------------------------------------------------------------------------

    public static IEnumerable<object[]> DetectionCases() => Golden.CaseNames("detection_lrt.json");

    [Theory]
    [MemberData(nameof(DetectionCases))]
    public void DetectionGlm_MatchesThePenalizedLrt(string name)
    {
        var c = Golden.Case("detection_lrt.json", name);
        var xFull = Golden.Mat(c, "x_full");
        var y = Golden.Vec(c, "y");
        var nS = y.Length;

        // Rebuild the call from the design: column 1 is the group indicator (so its 0s and 1s give
        // the two column lists) and any further column is a numeric covariate. The fixture's
        // covariate is already centered, and PRISM centers again - which is a no-op on an
        // already-centered column and leaves the design identical either way.
        var groupA = Enumerable.Range(0, nS).Where(i => xFull[i, 1] < 0.5).ToArray();
        var groupB = Enumerable.Range(0, nS).Where(i => xFull[i, 1] >= 0.5).ToArray();

        var covariates = new List<Covariate>();
        for (var col = 2; col < xFull.GetLength(1); col++)
        {
            var values = new double[nS];
            for (var i = 0; i < nS; i++)
                values[i] = xFull[i, col];
            covariates.Add(new NumericCovariate($"cov{col}", values));
        }

        var detection = new double[1, nS];
        for (var i = 0; i < nS; i++)
            detection[0, i] = y[i];

        var res = DetectionGlm.Run(detection, new[] { "PEP1" }, groupA, groupB,
            covariates.Count > 0 ? covariates : null);

        Assert.True(res.Identifiable, $"detection/{name}: design reported unidentifiable");
        var row = Assert.Single(res.Rows);

        // 1e-5, and the reason is worth stating: DetectionGlm calls FirthLogit.Fit at its default
        // tol of 1e-7 ON THE LOG-LIKELIHOOD, then forms the statistic by DIFFERENCING two such
        // log-likelihoods. Each is converged to ~1e-7, so their difference carries ~2e-7 of
        // absolute error, and the p-value inherits it. That is a property of the estimator's
        // stopping rule, not of the arithmetic here, and it is the floor on how precisely this
        // test can be asserted.
        Golden.Close(Golden.Num(c, "p"), row.P, 1e-5, $"detection/{name}/p");
    }

    // ------------------------------------------------------------------------------------
    // A divergence from inmoose, pinned deliberately rather than left to be discovered
    // ------------------------------------------------------------------------------------

    /// <summary>
    /// With 3 to 5 features the trend prior's spline df comes out at exactly 2 - one basis column
    /// per boundary knot and no interior knot - and <b>inmoose 0.9.1's <c>ns()</c> raises</b> on
    /// that basis rather than returning it. R's <c>splines::ns</c> builds it without complaint, and
    /// so does PRISM.
    ///
    /// <para>So there is no inmoose value to pin this against, and it is absent from
    /// <c>squeezevar.json</c> for that reason. What is asserted here is the part that does not need
    /// a reference: PRISM must produce a finite, usable prior instead of throwing, and must not
    /// quietly fall through to the global prior - which would be a different estimator giving
    /// different numbers under the same setting.</para>
    /// </summary>
    [Fact]
    public void SqueezeVarTrend_SplineDfTwo_IsSupportedWhereInmooseRaises()
    {
        var variances = new[] { 0.2, 0.5, 0.1, 0.8, 0.3 };
        var covariate = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

        var trend = EmpiricalBayes.SqueezeVarTrend(variances, 5.0, covariate);

        Assert.Equal(variances.Length, trend.VarPrior.Length);
        Assert.All(trend.VarPrior, v => Assert.True(v > 0 && double.IsFinite(v),
            $"expected a finite positive prior scale, got {v:R}"));
        Assert.All(trend.VarPost, v => Assert.True(v > 0 && double.IsFinite(v),
            $"expected a finite positive posterior variance, got {v:R}"));

        // The trend prior varies across features; the global prior is one number repeated. If the
        // spline had silently collapsed to the global path this would be a single distinct value.
        Assert.True(trend.VarPrior.Distinct().Count() > 1,
            "prior scale is constant across features - the spline path did not run");
    }
    /// <summary>
    /// The lab's intensity-trend variance prior, against the toolkit that defines it.
    /// </summary>
    /// <remarks>
    /// This is the one golden in the directory whose reference is another MacCoss Lab tool rather
    /// than a third-party library, and deliberately so: the estimator is not a published formula
    /// with an independent implementation to check against - it IS
    /// <c>proteomics-toolkit</c>'s <c>moderation="intensity_trend"</c>, and reproducing that is the
    /// whole requirement. It stands to PRISM as inmoose does for squeezeVar.
    /// <para>Note what is asserted: the per-feature prior SCALE only. The prior degrees of freedom
    /// are not part of this estimator - they stay at the global squeezeVar value - and that is
    /// precisely what distinguishes it from limma's <c>trend=TRUE</c>, which re-estimates both.</para>
    /// </remarks>
    [Theory]
    [MemberData(nameof(IntensityTrendCases))]
    public void IntensityTrendPrior_MatchesTheToolkit(string name)
    {
        var c = Golden.Case("intensity_trend.json", name);
        var exprLog2 = Golden.Mat(c, "expr_log2");
        var nA = c.GetProperty("n_a").GetInt32();
        var nB = c.GetProperty("n_b").GetInt32();
        var expected = Golden.Vec(c, "expected");

        var groups = new IReadOnlyList<int>[]
        {
            Enumerable.Range(0, nA).ToArray(),
            Enumerable.Range(nA, nB).ToArray(),
        };

        var actual = VariancePriors.IntensityTrend(exprLog2, groups);

        Assert.NotNull(actual);
        // 1e-9: the chain is a LOWESS fit, an interpolation and a delta-method division, all in
        // double precision with no iterative solve - so this is ordinary floating-point agreement,
        // not the looser tolerance squeezeVar's Newton step forces elsewhere in this file.
        Golden.CloseAll(expected, actual!, 1e-9, $"{name} prior scale");
    }

    public static IEnumerable<object[]> IntensityTrendCases() => Golden.CaseNames("intensity_trend.json");

    /// <summary>
    /// Welch, Student and Mann-Whitney against scipy, run through the public
    /// <see cref="Differential.Run"/> entry point so the dispatch is covered too.
    /// </summary>
    /// <remarks>
    /// Mann-Whitney is pinned to scipy's <c>method="asymptotic"</c>, not its <c>"auto"</c> default:
    /// PRISM implements only the normal approximation, and <c>auto</c> switches to the exact
    /// permutation distribution when the larger sample is 8 or fewer with no ties. Pinning
    /// <c>auto</c> would assert a branch PRISM does not have.
    /// </remarks>
    [Theory]
    [MemberData(nameof(SimpleTestCases))]
    public void SimpleTests_MatchScipy(string name)
    {
        var c = Golden.Case("simple_tests.json", name);
        var a = Golden.Vec(c, "a");
        var b = Golden.Vec(c, "b");

        // One feature, samples laid out as [A..., B...].
        var expr = new double[1, a.Length + b.Length];
        for (var i = 0; i < a.Length; i++)
            expr[0, i] = a[i];
        for (var i = 0; i < b.Length; i++)
            expr[0, a.Length + i] = b[i];
        var groupA = Enumerable.Range(0, a.Length).ToArray();
        var groupB = Enumerable.Range(a.Length, b.Length).ToArray();
        var ids = new[] { "f0" };

        DifferentialRow Run(DifferentialTest test) => Differential.Run(
            expr, ids, groupA, groupB,
            new DifferentialOptions { Test = test, Correction = MultipleTesting.None }).Rows[0];

        var welch = Run(DifferentialTest.WelchT);
        Golden.Close(Golden.Num(c, "welch_t"), welch.T, 1e-12, $"{name} welch t");
        Golden.Close(Golden.Num(c, "welch_p"), welch.PValue, 1e-12, $"{name} welch p");
        Golden.Close(Golden.Num(c, "logfc"), welch.LogFc, 1e-12, $"{name} logFC");

        var student = Run(DifferentialTest.StudentT);
        Golden.Close(Golden.Num(c, "student_t"), student.T, 1e-12, $"{name} student t");
        Golden.Close(Golden.Num(c, "student_p"), student.PValue, 1e-12, $"{name} student p");

        var mw = Run(DifferentialTest.MannWhitney);
        Golden.Close(Golden.Num(c, "mw_u"), mw.T, 1e-12, $"{name} mann-whitney U");
        Golden.Close(Golden.Num(c, "mw_p"), mw.PValue, 1e-12, $"{name} mann-whitney p");
        // A rank test reports a median shift, not a mean difference.
        Golden.Close(Golden.Num(c, "median_diff"), mw.LogFc, 1e-12, $"{name} median diff");
    }

    public static IEnumerable<object[]> SimpleTestCases() => Golden.CaseNames("simple_tests.json");

    /// <summary>
    /// The multiple-testing methods beside BH, against statsmodels. Every case is NaN-free, where
    /// PRISM's own NaN policy and statsmodels' agree; <c>FdrTests</c> pins where they do not.
    /// </summary>
    [Theory]
    [MemberData(nameof(CorrectionCases))]
    public void Corrections_MatchStatsmodels(string name)
    {
        var c = Golden.Case("corrections.json", name);
        var p = Golden.Vec(c, "p");

        Golden.CloseAll(Golden.Vec(c, "by"),
            Fdr.Adjust(p, MultipleTesting.BenjaminiYekutieli), 1e-12, $"{name} BY");
        Golden.CloseAll(Golden.Vec(c, "bonferroni"),
            Fdr.Adjust(p, MultipleTesting.Bonferroni), 1e-12, $"{name} bonferroni");
        Golden.CloseAll(Golden.Vec(c, "holm"),
            Fdr.Adjust(p, MultipleTesting.Holm), 1e-12, $"{name} holm");
    }

    public static IEnumerable<object[]> CorrectionCases() => Golden.CaseNames("corrections.json");
}
