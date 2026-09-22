#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy==2.5.3",
#     "scipy==1.18.1",
#     "statsmodels==0.15.0",
#     "inmoose==0.9.1",
# ]
# ///
"""Generate the differential-analysis golden fixtures from the reference implementations.

Unlike the `refanchored/` and `sva/` fixtures, these are NOT frozen: every reference here is a
maintained third-party library on PyPI, so the goldens can be regenerated at any time and the
dependency versions are pinned in the PEP 723 header above.

Which library is the reference for which quantity:

| PRISM type                        | reference                                              |
|-----------------------------------|--------------------------------------------------------|
| `Fdr.BenjaminiHochberg`           | `statsmodels.stats.multitest.multipletests("fdr_bh")`   |
| `EmpiricalBayes.Trigamma`         | `scipy.special.polygamma(1, .)`                         |
| `EmpiricalBayes.Tetragamma`       | `scipy.special.polygamma(2, .)`                         |
| `EmpiricalBayes.TrigammaInverse`  | `scipy.special.polygamma(1, .)` inverted (round trip)   |
| `EmpiricalBayes.SqueezeVar*`      | `inmoose.limma.squeezeVar`                              |
| `NaturalSplineBasis.Build`        | `inmoose.utils.splines.ns`                              |
| `LinearModel.Fit`                 | `numpy.linalg.lstsq` + the textbook OLS formulas        |
| `Differential.Run`                | the four above, composed the way limma composes them    |
| `Detection.FisherExact`           | `scipy.stats.fisher_exact(alternative="two-sided")`     |
| `Detection.FirthLogit`            | `scipy.optimize` on the penalized log-likelihood        |
| `DifferentialPca.Compute`         | `numpy.linalg.svd(full_matrices=False)`                 |
| `Detection.DetectionGlm`          | penalized LRT: `scipy.optimize` twice + `scipy.stats.chi2` |

Nothing here imports PRISM. The point of a golden is that it was produced without reference to the
code under test, so a shared mistake cannot cancel out.

`FirthLogit` is the one entry with no library implementation to call. Rather than pin it to the
sibling Python implementation it was ported from - which would only prove the two agree - the
reference maximizes the Jeffreys-penalized log-likelihood directly with a derivative-free
optimizer. Different objective formulation, different algorithm, same fixed point.

Run from the repository root:

    uv run dotnet/tests/fixtures/differential/generate.py

or, without uv, in an environment holding the pinned versions above:

    python dotnet/tests/fixtures/differential/generate.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.special
import scipy.stats
from inmoose.limma import squeezeVar
from inmoose.utils.splines import ns
from statsmodels.stats.multitest import multipletests

OUT = Path("dotnet/tests/fixtures/differential")


# --------------------------------------------------------------------------------------------
# Encoding
#
# Every floating-point value is written as a STRING holding Python's shortest round-trip repr.
# Two reasons, both of which bit a golden elsewhere in this repo:
#   * bare NaN / Infinity are not valid JSON, and System.Text.Json rejects them outright;
#   * a float written through the JSON number path is at the mercy of both encoders' shortest-repr
#     rules, and a golden that is only nearly bit-exact cannot be asserted at 1e-15.
# Strings sidestep both: Python writes repr(), .NET reads double.Parse(InvariantCulture), and the
# value that comes back is the same 64 bits that went in.
# --------------------------------------------------------------------------------------------


def num(x) -> str:
    """One float as its exact round-trip string."""
    return repr(float(x))


def vec(a) -> list[str]:
    return [num(v) for v in np.asarray(a, dtype=float).ravel()]


def mat(a) -> list[list[str]]:
    a = np.asarray(a, dtype=float)
    return [[num(v) for v in row] for row in a]


def write(name: str, payload: dict) -> None:
    path = OUT / name
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


class Rng:
    """A tiny reproducible LCG.

    numpy's Generator is reproducible too, but its stream is a guarantee about numpy, not about
    this fixture; a bare LCG keeps the inputs readable and pins them to this file.
    """

    def __init__(self, seed: int) -> None:
        self._s = seed | 1

    def next_double(self) -> float:
        self._s = (self._s * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return ((self._s >> 11) & ((1 << 53) - 1)) / float(1 << 53)

    def normal(self) -> float:
        # Box-Muller; one draw per call is wasteful and entirely beside the point here.
        u1 = max(self.next_double(), 1e-12)
        u2 = self.next_double()
        return float(np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2))

    def matrix(self, rows: int, cols: int, loc: float = 10.0, scale: float = 1.0) -> np.ndarray:
        return np.array(
            [[loc + scale * self.normal() for _ in range(cols)] for _ in range(rows)],
            dtype=float,
        )


# --------------------------------------------------------------------------------------------
# Benjamini-Hochberg
# --------------------------------------------------------------------------------------------


def gen_fdr() -> None:
    rng = Rng(11)
    cases = []

    def add(name: str, p, note: str) -> None:
        p = np.asarray(p, dtype=float)
        expected = multipletests(p, method="fdr_bh")[1]
        cases.append({"name": name, "note": note, "p": vec(p), "expected": vec(expected)})

    add("uniform_steps", [0.01, 0.02, 0.03, 0.04, 0.05], "every step-up value equal; all take the max")
    add("one_strong_hit", [1e-8, 0.2, 0.4, 0.6, 0.8, 0.9], "a single dominant hit")
    add("ties", [0.04, 0.04, 0.04, 0.2, 0.2, 0.9], "tied p-values must take an identical q")
    add("all_large", [0.6, 0.7, 0.8, 0.9, 0.95], "everything clips to 1")
    add("single", [0.031], "m = 1: the adjusted value is the raw one")
    add("descending", [0.9, 0.5, 0.2, 0.05, 0.001], "input order is not sorted order")
    add("many", [rng.next_double() ** 3 for _ in range(200)], "200 values, heavily skewed to small p")
    add("with_exact_zero", [0.0, 0.0, 0.3, 0.7], "exact zeros stay zero")

    write(
        "fdr.json",
        {
            "reference": "statsmodels.stats.multitest.multipletests(method='fdr_bh')",
            "note": (
                "PRISM's BenjaminiHochberg additionally passes NaN through and excludes it from m, "
                "where statsmodels returns all-NaN. That behavior is PRISM's own and is pinned by "
                "FdrTests, not here - every case below is NaN-free, where the two agree exactly."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Polygamma
# --------------------------------------------------------------------------------------------


def gen_polygamma() -> None:
    xs = [
        0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0,
        # Either side of the asymptotic-series threshold (30), where the recurrence hands over.
        29.0, 29.999, 30.0, 30.001, 31.0, 50.0, 100.0, 500.0, 5000.0, 1e6,
        # Half-integers are what fitFDist actually evaluates: trigamma(df/2).
        2.5, 3.5, 25.5, 60.5,
    ]
    trigamma = [scipy.special.polygamma(1, x) for x in xs]
    tetragamma = [scipy.special.polygamma(2, x) for x in xs]

    # trigammaInverse is pinned by round-trip: solve trigamma(y) = x, then check trigamma(y) == x.
    # Inverting through scipy rather than through limma's own Newton keeps the reference independent.
    inv_inputs = [1e-8, 1e-6, 1e-4, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3, 1e6, 1e7, 1e8, 1e10]
    inv_expected = []
    for x in inv_inputs:
        if x > 1e7:
            inv_expected.append(1.0 / np.sqrt(x))
        elif x < 1e-6:
            inv_expected.append(1.0 / x)
        else:
            root = scipy.optimize.brentq(
                lambda y: float(scipy.special.polygamma(1, y)) - x,
                1e-12, 1e12, xtol=1e-300, rtol=8.9e-16, maxiter=500,
            )
            inv_expected.append(root)

    write(
        "polygamma.json",
        {
            "reference": "scipy.special.polygamma(1|2, .); trigammaInverse via scipy.optimize.brentq",
            "note": (
                "The two smallest/largest trigammaInverse inputs take limma's closed-form branches "
                "(1/x and 1/sqrt(x)) rather than Newton, and are reproduced here the same way."
            ),
            "x": vec(xs),
            "trigamma": vec(trigamma),
            "tetragamma": vec(tetragamma),
            "trigamma_inverse_x": vec(inv_inputs),
            "trigamma_inverse_expected": vec(inv_expected),
        },
    )


# --------------------------------------------------------------------------------------------
# squeezeVar (global and intensity-trend priors)
# --------------------------------------------------------------------------------------------


def gen_squeezevar() -> None:
    rng = Rng(23)
    cases = []

    def add(name: str, variances, df: float, covariate=None, note: str = "") -> None:
        v = np.asarray(variances, dtype=float)
        cov = None if covariate is None else np.asarray(covariate, dtype=float)
        res = squeezeVar(v, df, covariate=cov)
        cases.append(
            {
                "name": name,
                "note": note,
                "variances": vec(v),
                "df_residual": num(df),
                "covariate": None if cov is None else vec(cov),
                "df_prior": num(np.atleast_1d(res["df_prior"])[0]),
                "var_prior": vec(np.atleast_1d(res["var_prior"])),
                "var_post": vec(np.atleast_1d(res["var_post"])),
            }
        )

    # --- global prior -------------------------------------------------------------------
    add(
        "global_finite_prior",
        [0.05, 0.5, 0.1, 1.2, 0.3, 0.8, 0.02, 0.15, 2.0, 0.25,
         0.6, 0.09, 0.35, 1.5, 0.04, 0.7, 0.12, 0.9, 0.2, 0.45],
        7.0,
        note="the ordinary case: a finite, small prior df",
    )
    add(
        "global_large_prior_df",
        [0.1, 0.2, 0.4, 0.5],
        5.0,
        note="tightly clustered variances drive df_prior to ~51, the regime the polygamma "
             "asymptotic threshold has to hold 1e-9 on",
    )
    add(
        "global_uniform_variance",
        [0.25] * 12,
        6.0,
        note="identical variances: evar <= 0, so df_prior is infinite and every posterior is the "
             "pooled mean",
    )
    add(
        "global_wide_spread",
        [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0],
        4.0,
        note="seven orders of magnitude; exercises the log-F fit at its limits",
    )
    add(
        "global_many_features",
        [abs(rng.normal()) * 0.4 + 0.02 for _ in range(500)],
        9.0,
        note="500 features, the realistic size",
    )

    # --- intensity-trend prior ----------------------------------------------------------
    n = 60
    amean = np.array([6.0 + 8.0 * i / (n - 1) for i in range(n)])
    # A genuine mean-variance trend: variance falls as intensity rises, which is the whole reason
    # limma-trend exists.
    trend_var = np.array(
        [float(np.exp(-0.35 * (a - 6.0)) * (0.8 + 0.4 * rng.next_double())) for a in amean]
    )
    add("trend_decreasing", trend_var, 8.0, covariate=amean,
        note="variance decays with mean intensity - the limma-trend motivating case")

    n2 = 40
    amean2 = np.array([4.0 + 0.25 * i for i in range(n2)])
    flat_var = np.array([0.3 + 0.05 * rng.normal() for _ in range(n2)])
    flat_var = np.abs(flat_var)
    add("trend_flat", flat_var, 6.0, covariate=amean2,
        note="no real trend: the spline should flatten and land near the global prior")

    # splinedf = 1 + (n>=3) + (n>=6) + (n>=30), capped at the number of distinct covariate values.
    add("trend_three_distinct", [0.2, 0.5, 0.1, 0.8, 0.3, 0.4], 5.0,
        covariate=[1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
        note="only 3 distinct covariate values: splinedf caps at 3")
    add("trend_six_even", [0.2, 0.5, 0.1, 0.8, 0.3, 0.4], 5.0,
        covariate=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        note="n=6 is exactly where the (nok >= 6) term turns splinedf from 2 to 3")

    # NOT covered here, deliberately: splinedf == 2, which happens for 3 <= n <= 5 features.
    # inmoose 0.9.1's ns() raises on a zero-interior-knot basis (see README), so there is no
    # reference value to pin against. R's splines::ns handles it, and so does PRISM; that case is
    # pinned as PRISM-only behavior in EmpiricalBayesTests instead.

    write(
        "squeezevar.json",
        {
            "reference": "inmoose.limma.squeezeVar (inmoose 0.9.1, a port of limma 3.55.1)",
            "note": (
                "var_prior is length 1 for the global prior and length n for the trend prior, "
                "exactly as inmoose returns it."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Natural spline basis
# --------------------------------------------------------------------------------------------


def gen_spline() -> None:
    """Pin the spline by its COLUMN SPAN, not its entries.

    R's `ns` (and inmoose's port) fixes the basis only up to the orthogonal rotation the QR of the
    natural-boundary constraint happens to produce, so two correct implementations can return
    different matrices. What the trend fit depends on - and all it depends on - is the span, and
    the orthogonal projector onto the span, H = B (B'B)^-1 B', is invariant to that rotation.
    Comparing H is therefore both the strictest test available and the only one that is not
    asserting an arbitrary choice.
    """
    rng = Rng(37)
    cases = []

    def add(name: str, x, df: int, note: str) -> None:
        x = np.asarray(x, dtype=float)
        basis = np.asarray(ns(x, df=df, include_intercept=True).basis, dtype=float)
        q, _ = np.linalg.qr(basis)
        hat = q @ q.T
        cases.append(
            {
                "name": name,
                "note": note,
                "x": vec(x),
                "df": df,
                "n_basis_columns": int(basis.shape[1]),
                "hat": mat(hat),
            }
        )

    add("even_grid_df4", [4.0 + 0.5 * i for i in range(18)], 4, "evenly spaced, the default df")
    add("even_grid_df3", [4.0 + 0.5 * i for i in range(18)], 3, "one interior knot fewer")
    add("skewed_df4", [float(10.0 * rng.next_double() ** 3) for _ in range(20)], 4,
        "clustered near the lower boundary, so the interior knots are not evenly spaced")
    add("with_duplicates_df3", [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0], 3,
        "repeated covariate values; percentile knots can coincide with data points")
    add("five_points_df3", [1.0, 2.0, 3.0, 4.0, 5.0], 3, "the smallest basis with an interior knot")

    # df=2 (zero interior knots) is absent because inmoose 0.9.1's ns() raises on it. See README.

    write(
        "spline.json",
        {
            "reference": "inmoose.utils.splines.ns(x, df, include_intercept=True), inmoose 0.9.1",
            "note": (
                "`hat` is the orthogonal projector onto the basis column span, B (B'B)^-1 B'. The "
                "basis itself is defined only up to an orthogonal rotation, so its entries are NOT "
                "a parity target; the span is."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# lmFit
# --------------------------------------------------------------------------------------------


def gen_lmfit() -> None:
    rng = Rng(51)
    cases = []

    def add(name: str, expr, design, note: str) -> None:
        expr = np.asarray(expr, dtype=float)      # features x samples
        design = np.asarray(design, dtype=float)  # samples x coef
        n, p = design.shape
        beta, _, _, _ = np.linalg.lstsq(design, expr.T, rcond=None)
        beta = beta.T
        resid = expr - beta @ design.T
        df_res = float(n - p)
        sigma = np.sqrt((resid**2).sum(axis=1) / df_res)
        stdev_unscaled = np.sqrt(np.diag(np.linalg.inv(design.T @ design)))
        cases.append(
            {
                "name": name,
                "note": note,
                "expr": mat(expr),
                "design": mat(design),
                "df_residual": num(df_res),
                "coefficients": mat(beta),
                "sigma": vec(sigma),
                "amean": vec(expr.mean(axis=1)),
                "stdev_unscaled": vec(stdev_unscaled),
            }
        )

    # Two groups of 4, design [intercept, groupB].
    design2 = np.column_stack([np.ones(8), np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)])
    add("two_groups", rng.matrix(15, 8), design2, "the plain two-group contrast")

    # Two groups plus a centered numeric covariate and a dummy-coded categorical one.
    age = np.array([61.0, 55.0, 70.0, 48.0, 66.0, 59.0, 73.0, 52.0])
    sex = np.array([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0])
    design_cov = np.column_stack([design2, age - age.mean(), sex])
    add("with_covariates", rng.matrix(12, 8), design_cov,
        "intercept + group + centered numeric + dummy categorical")

    # An ill-conditioned design: the covariate nearly tracks the group. This is where forming the
    # normal equations instead of a QR would visibly lose digits.
    near = np.array([0.0, 0.0, 0.0, 0.001, 1.0, 1.0, 1.0, 0.999])
    design_ill = np.column_stack([design2, near])
    add("ill_conditioned", rng.matrix(10, 8), design_ill,
        "covariate almost collinear with group; QR vs normal equations diverge here")

    # Unbalanced groups, larger n.
    grp = np.array([0.0] * 5 + [1.0] * 12)
    design_unbal = np.column_stack([np.ones(17), grp])
    add("unbalanced", rng.matrix(30, 17), design_unbal, "5 vs 12")

    write(
        "lmfit.json",
        {
            "reference": "numpy.linalg.lstsq plus the textbook OLS formulas (limma lmFit)",
            "note": (
                "`expr` is features x samples (PRISM's orientation) and `design` is samples x coef. "
                "stdev_unscaled is sqrt(diag((X'X)^-1)) and depends only on the design."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# End-to-end moderated t
# --------------------------------------------------------------------------------------------


def moderated_t(expr: np.ndarray, design: np.ndarray, coef_idx: int, trend: bool):
    """limma's lmFit -> squeezeVar -> moderated t -> BH, composed from the reference libraries.

    This is the composition the limma vignette describes, written out rather than called through
    a wrapper so that every step is visibly one of the pinned references.
    """
    n, p = design.shape
    beta, _, _, _ = np.linalg.lstsq(design, expr.T, rcond=None)
    beta = beta.T
    resid = expr - beta @ design.T
    df_res = float(n - p)
    sigma2 = (resid**2).sum(axis=1) / df_res
    stdev_unscaled = np.sqrt(np.diag(np.linalg.inv(design.T @ design)))
    amean = expr.mean(axis=1)

    sv = squeezeVar(sigma2, df_res, covariate=amean if trend else None)
    df_prior = float(np.atleast_1d(sv["df_prior"])[0])
    var_post = np.atleast_1d(sv["var_post"])

    df_total = df_res + df_prior
    t = beta[:, coef_idx] / (stdev_unscaled[coef_idx] * np.sqrt(var_post))
    if np.isinf(df_total):
        pval = 2.0 * scipy.stats.norm.cdf(-np.abs(t))
    else:
        pval = 2.0 * scipy.stats.t.cdf(-np.abs(t), df=df_total)
    adj = multipletests(pval, method="fdr_bh")[1]
    return {
        "logfc": beta[:, coef_idx],
        "amean": amean,
        "t": t,
        "p": pval,
        "adj_p": adj,
        "df_residual": df_res,
        "df_prior": df_prior,
    }


def gen_moderated_t() -> None:
    rng = Rng(67)
    cases = []

    def add(name: str, expr, n_a: int, design, trend: bool, note: str,
            covariates=None) -> None:
        expr = np.asarray(expr, dtype=float)
        design = np.asarray(design, dtype=float)
        res = moderated_t(expr, design, coef_idx=1, trend=trend)
        n_b = expr.shape[1] - n_a
        cases.append(
            {
                "name": name,
                "note": note,
                "expr": mat(expr),
                "n_a": n_a,
                "n_b": n_b,
                "design": mat(design),
                # The covariates as PRISM receives them - RAW, uncentered, in [A..., B...] order.
                # PRISM centers numerics and dummy-codes categoricals itself, so handing it the
                # design columns directly would test the arithmetic twice and the design-building
                # code not at all.
                "covariates": covariates or [],
                "trend": trend,
                "df_residual": num(res["df_residual"]),
                "df_prior": num(res["df_prior"]),
                "logfc": vec(res["logfc"]),
                "amean": vec(res["amean"]),
                "t": vec(res["t"]),
                "p": vec(res["p"]),
                "adj_p": vec(res["adj_p"]),
            }
        )

    # --- no covariates, heterogeneous variance (finite prior df) ------------------------
    n_a, n_b = 6, 7
    grp = np.array([0.0] * n_a + [1.0] * n_b)
    design = np.column_stack([np.ones(n_a + n_b), grp])
    # Per-feature noise scales must genuinely differ, or every residual variance comes out the same,
    # evar <= 0, and the prior df is infinite - which is a branch worth testing but not the typical
    # one. Spreading the scale over an order of magnitude gives the finite-prior case.
    expr = np.array(
        [[18.0 + (0.15 + 1.5 * (i / 39.0) ** 2) * rng.normal() for _ in range(n_a + n_b)]
         for i in range(40)]
    )
    # Put a real effect into the first five features so the volcano is not all noise.
    expr[:5, n_a:] += 1.7
    add("no_covariates", expr, n_a, design, False,
        "6 vs 7, heterogeneous per-feature noise, a real effect in the first 5 features")

    # --- uniform variance drives the prior df to infinity --------------------------------
    base = np.array([[10.0, 10.5, 11.0, 12.0, 12.5, 13.0]] * 8)
    add("infinite_prior_df", base, 3, np.column_stack([np.ones(6), [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]),
        False, "every feature the same row: evar <= 0, df_prior is infinite, p comes from the normal")

    # --- intensity trend ------------------------------------------------------------------
    n_a2, n_b2 = 8, 8
    grp2 = np.array([0.0] * n_a2 + [1.0] * n_b2)
    design2 = np.column_stack([np.ones(16), grp2])
    # Abundance spans a wide dynamic range and the noise shrinks with it: limma-trend's premise.
    rows = []
    for i in range(50):
        level = 8.0 + 12.0 * i / 49.0
        noise = float(np.exp(-0.12 * (level - 8.0)))
        rows.append([level + noise * rng.normal() for _ in range(16)])
    expr_tr = np.array(rows)
    expr_tr[:6, n_a2:] += 0.9
    add("trend_prior", expr_tr, n_a2, design2, True,
        "wide dynamic range with intensity-dependent noise; the trend prior is the point")
    add("trend_prior_off", expr_tr, n_a2, design2, False,
        "the same matrix with the global prior, so the two priors can be told apart")

    # --- covariate adjustment ---------------------------------------------------------------
    n_a3, n_b3 = 7, 7
    grp3 = np.array([0.0] * n_a3 + [1.0] * n_b3)
    age3 = np.array([61.0, 55.0, 70.0, 48.0, 66.0, 59.0, 73.0,
                     52.0, 64.0, 58.0, 71.0, 49.0, 67.0, 60.0])
    # Sorted levels, first dropped: F is the reference level, so the dummy marks M.
    sex3 = np.array([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
                     0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0])
    design3 = np.column_stack([np.ones(14), grp3, age3 - age3.mean(), sex3])
    expr3 = rng.matrix(35, 14, loc=15.0, scale=0.8)
    # Make age genuinely predictive, so dropping the covariate would change the answer.
    expr3 += 0.03 * (age3 - age3.mean())[None, :]
    expr3[:4, n_a3:] += 1.2
    add("with_covariates", expr3, n_a3, design3, False,
        "centered age + sex dummy; age has a real effect so the adjustment matters",
        covariates=[
            {"name": "age", "kind": "numeric", "values": vec(age3)},
            # Sorted levels are [F, M] and the first is dropped, so the dummy marks M - which is
            # exactly the 1s in sex3.
            {"name": "sex", "kind": "categorical",
             "values": ["M" if v == 1.0 else "F" for v in sex3]},
        ])

    # --- a single feature, the smallest possible run -------------------------------------
    add("single_feature", np.array([[10.0, 10.2, 9.8, 12.0, 12.4, 11.6]]), 3,
        np.column_stack([np.ones(6), [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]), False,
        "n_features = 1: fitFDist takes its nok == 1 branch (df_prior 0)")

    write(
        "moderated_t.json",
        {
            "reference": (
                "numpy.linalg.lstsq + inmoose.limma.squeezeVar + scipy.stats.t.cdf + "
                "statsmodels multipletests('fdr_bh'), composed as limma composes them"
            ),
            "note": (
                "`design` is given explicitly and its column 1 is the group indicator, so the C# "
                "side must reproduce it from (n_a, n_b) and the covariates: intercept, groupB, "
                "centered numeric covariates, then dummy-coded categoricals with the first sorted "
                "level dropped. Samples are ordered [A..., B...]."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Fisher exact
# --------------------------------------------------------------------------------------------


def gen_fisher() -> None:
    tables = [
        (8, 2, 1, 9), (0, 10, 10, 0), (5, 5, 5, 5), (1, 0, 0, 1),
        (0, 0, 4, 6), (12, 3, 4, 11), (30, 10, 12, 28), (1, 19, 18, 2),
        (7, 0, 0, 7), (2, 8, 8, 2), (100, 50, 60, 90), (3, 1, 1, 3),
        (0, 0, 0, 0), (25, 0, 0, 25), (14, 6, 9, 11),
    ]
    cases = []
    for a, b, c, d in tables:
        _, p = scipy.stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        cases.append({"a": a, "b": b, "c": c, "d": d, "p": num(p)})

    write(
        "fisher.json",
        {
            "reference": "scipy.stats.fisher_exact(alternative='two-sided')",
            "note": "Table is [[a, b], [c, d]]: detected/not-detected by group.",
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Firth-penalized logistic regression
# --------------------------------------------------------------------------------------------


def firth_reference(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Maximize l(b) + 0.5 log det(X' W X) with a derivative-free optimizer.

    Deliberately not Newton-with-hat-matrix-adjustment, which is what the implementation under
    test does: an independent algorithm reaching the same stationary point is evidence; the same
    algorithm reaching it twice is not.
    """

    def penalized(b: np.ndarray) -> float:
        eta = x @ b
        # log(1 + exp(eta)) computed stably; eta can be large under separation, which is the case
        # Firth exists for.
        log1pexp = np.logaddexp(0.0, eta)
        ll = float(np.sum(y * eta - log1pexp))
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        info = x.T @ (w[:, None] * x)
        sign, logdet = np.linalg.slogdet(info)
        if sign <= 0:
            return -np.inf
        return ll + 0.5 * logdet

    def objective(b: np.ndarray) -> float:
        v = penalized(b)
        return np.inf if not np.isfinite(v) else -v

    best = np.zeros(x.shape[1])
    # Nelder-Mead first (no derivatives, no curvature assumptions), then Powell to polish. Both are
    # direct-search methods, so neither reuses the Fisher information the estimator under test
    # builds its step from.
    for method, opts in (
        ("Nelder-Mead", {"xatol": 1e-14, "fatol": 1e-14, "maxiter": 100000, "maxfev": 100000}),
        ("Powell", {"xtol": 1e-14, "ftol": 1e-14, "maxiter": 100000, "maxfev": 100000}),
        ("Nelder-Mead", {"xatol": 1e-15, "fatol": 1e-15, "maxiter": 100000, "maxfev": 100000}),
    ):
        res = scipy.optimize.minimize(objective, best, method=method, options=opts)
        best = res.x
    return best, penalized(best)


def gen_firth() -> None:
    cases = []

    def add(name: str, x, y, note: str) -> None:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        beta, ll = firth_reference(x, y)
        cases.append(
            {"name": name, "note": note, "x": mat(x), "y": vec(y),
             "beta": vec(beta), "penalized_loglik": num(ll)}
        )

    add(
        "separable",
        [[1.0, -2.0], [1.0, -1.0], [1.0, -0.5], [1.0, 0.5], [1.0, 1.0], [1.0, 2.0]],
        [0, 0, 0, 1, 1, 1],
        "perfect separation: the unpenalized MLE diverges, Firth stays finite",
    )
    add(
        "mixed",
        [[1.0, -2.0], [1.0, -1.0], [1.0, 0.0], [1.0, 1.0],
         [1.0, 2.0], [1.0, 3.0], [1.0, 4.0], [1.0, 5.0]],
        [0, 1, 0, 0, 1, 1, 0, 1],
        "overlapping classes: an ordinary well-posed fit",
    )
    add(
        "two_predictors",
        [[1.0, 0.0, 1.2], [1.0, 1.0, -0.4], [1.0, 0.0, 0.7], [1.0, 1.0, 2.1],
         [1.0, 0.0, -1.1], [1.0, 1.0, 0.3], [1.0, 0.0, 1.9], [1.0, 1.0, -0.8],
         [1.0, 0.0, 0.1], [1.0, 1.0, 1.4]],
        [0, 1, 0, 1, 0, 0, 1, 0, 1, 1],
        "a group indicator plus a continuous covariate, the adjusted-detection design",
    )
    add(
        "quasi_separation",
        [[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0], [1.0, 5.0], [1.0, 5.0]],
        [0, 0, 0, 1, 1, 0],
        "one point on the wrong side of an otherwise clean split",
    )
    add(
        "all_one_class",
        [[1.0, -1.0], [1.0, 0.0], [1.0, 1.0], [1.0, 2.0]],
        [1, 1, 1, 1],
        "no contrast at all: the intercept runs away and only the penalty holds it",
    )

    write(
        "firth.json",
        {
            "reference": (
                "scipy.optimize (Nelder-Mead then Powell) maximizing "
                "l(b) + 0.5*log det(X' diag(p(1-p)) X)"
            ),
            "note": (
                "A direct-search optimum is accurate in the OBJECTIVE to roughly the tolerance "
                "asked for, but the coefficients themselves are only as well determined as the "
                "curvature allows - a flat penalized likelihood (all_one_class) pins beta far more "
                "loosely than it pins the log-likelihood. Assert the log-likelihood tightly and the "
                "coefficients loosely; that is what the fixture is able to support."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Sample-space PCA
# --------------------------------------------------------------------------------------------


def gen_pca() -> None:
    rng = Rng(83)
    cases = []

    def add(name: str, expr, note: str, n_components: int = 6) -> None:
        expr = np.asarray(expr, dtype=float)  # features x samples
        n_samples = expr.shape[1]
        # Complete-case, then center each feature across the selected samples (no scaling).
        complete = expr[~np.isnan(expr).any(axis=1), :]
        centered = complete - complete.mean(axis=1, keepdims=True)
        # numpy's SVD of the sample x feature matrix; scores are U * S.
        u, s, _ = np.linalg.svd(centered.T, full_matrices=False)
        k = min(n_components, n_samples, complete.shape[0])
        scores = u[:, :k] * s[:k]
        total = float((s**2).sum())
        ratio = (s[:k] ** 2) / total
        cases.append(
            {
                "name": name,
                "note": note,
                "expr": mat(expr),
                "n_components": n_components,
                "n_features_used": int(complete.shape[0]),
                # Signs are arbitrary per component, so the fixture carries |scores| and the test
                # compares magnitudes. A sign flip is not a defect; a magnitude change is.
                "abs_scores": mat(np.abs(scores)),
                "variance_ratio": vec(ratio),
            }
        )

    add("dense", rng.matrix(60, 9, loc=14.0, scale=1.5), "no missing values")

    # A matrix with structure: two clusters of samples.
    struct = rng.matrix(80, 10, loc=12.0, scale=0.4)
    struct[:, 5:] += 2.0
    add("two_clusters", struct, "PC1 should separate the two groups cleanly")

    # Missing values: whole features drop out under the complete-case rule.
    with_nan = rng.matrix(40, 8, loc=11.0, scale=1.0)
    with_nan[3, 2] = np.nan
    with_nan[17, 0] = np.nan
    with_nan[17, 7] = np.nan
    with_nan[39, 4] = np.nan
    add("with_missing", with_nan, "4 features carry a NaN and are dropped entirely")

    # Fewer features than samples: the component count is capped by rank, not by n_components.
    add("more_samples_than_features", rng.matrix(4, 9, loc=10.0, scale=1.0),
        "only 4 features, so at most 4 components exist")

    write(
        "pca.json",
        {
            "reference": "numpy.linalg.svd(full_matrices=False) on the centered complete-case matrix",
            "note": (
                "`expr` is features x samples and every column is selected. Component signs are "
                "arbitrary, so the golden holds |scores|."
            ),
            "cases": cases,
        },
    )


# --------------------------------------------------------------------------------------------
# Detection GLM
# --------------------------------------------------------------------------------------------


def gen_detection_lrt() -> None:
    """The Firth-penalized likelihood-ratio test behind the covariate-adjusted Detection view.

    `DetectionGlm` does NOT fit an unpenalized GLM, so statsmodels' `GLM(family=Binomial())` is the
    wrong reference for it: it fits Firth twice - once on the full design, once with the group
    column removed - and refers 2*(ll_full - ll_reduced) to chi2 with 1 df. The penalized
    log-likelihood is what firth.json already pins independently; this adds the two-fit arithmetic
    on top of it, with the tail probability from scipy.
    """
    cases = []

    def add(name: str, x_full, y, note: str) -> None:
        x_full = np.asarray(x_full, dtype=float)
        y = np.asarray(y, dtype=float)
        # PRISM drops design column 1 (the group indicator) to form the reduced model.
        x_red = np.delete(x_full, 1, axis=1)
        _, ll_full = firth_reference(x_full, y)
        _, ll_red = firth_reference(x_red, y)
        lrt = 2.0 * (ll_full - ll_red)
        pval = float(scipy.stats.chi2.sf(max(lrt, 0.0), 1))
        cases.append(
            {
                "name": name,
                "note": note,
                "x_full": mat(x_full),
                "y": vec(y),
                "ll_full": num(ll_full),
                "ll_reduced": num(ll_red),
                "lrt": num(lrt),
                "p": num(pval),
            }
        )

    n = 24
    grp = np.array([0.0] * 12 + [1.0] * 12)
    age = np.array([50.0 + 1.7 * i for i in range(n)])
    age_c = age - age.mean()

    y1 = np.array([0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                   1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1], dtype=float)
    add("group_effect_with_covariate", np.column_stack([np.ones(n), grp, age_c]), y1,
        "detection commoner in group B, adjusted for a centered numeric covariate")

    y2 = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
                   1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0], dtype=float)
    add("no_group_effect", np.column_stack([np.ones(n), grp, age_c]), y2,
        "detection unrelated to group; the LRT should be near zero")

    y3 = np.array([0] * 12 + [1] * 12, dtype=float)
    add("complete_separation", np.column_stack([np.ones(n), grp, age_c]), y3,
        "detected in every B sample and no A sample - the case an unpenalized GLM cannot fit at all")

    add("group_only", np.column_stack([np.ones(n), grp]), y1,
        "no covariate: the reduced model is intercept-only")

    write(
        "detection_lrt.json",
        {
            "reference": (
                "scipy.optimize on the Jeffreys-penalized log-likelihood (twice) + "
                "scipy.stats.chi2.sf(., 1)"
            ),
            "note": (
                "Design column 1 is the group indicator and is the column the reduced model drops. "
                "Assert `p` and `lrt`, which are well determined; the individual log-likelihoods "
                "are reported for diagnosis. statsmodels' unpenalized GLM is deliberately NOT the "
                "reference here - it diverges on complete_separation, which is the case this path "
                "exists to handle."
            ),
            "cases": cases,
        },
    )


def main() -> None:
    if not OUT.is_dir():
        raise SystemExit(f"run from the repository root: {OUT} not found")
    gen_fdr()
    gen_polygamma()
    gen_squeezevar()
    gen_spline()
    gen_lmfit()
    gen_moderated_t()
    gen_fisher()
    gen_firth()
    gen_pca()
    gen_detection_lrt()


if __name__ == "__main__":
    main()
