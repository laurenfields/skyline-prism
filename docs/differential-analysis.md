# Differential analysis

The Skyline tool's **Differential** pane runs a limma empirical-Bayes moderated-t contrast on a
PRISM output directory: protein or peptide level, two groups from any metadata or clinical column,
optional covariate adjustment, and optionally the intensity-trend variance prior (limma's
`trend=TRUE`). The statistics live in `SkylinePrism.Core.DifferentialAnalysis` and are reachable
without the GUI.

## What each quantity is held to

Every component with a published reference implementation is pinned to it by committed goldens in
`dotnet/tests/fixtures/differential/`, regenerable with
`uv run dotnet/tests/fixtures/differential/generate.py`. See that directory's README for the
per-quantity reference table and for the tolerances, which are not uniform and are not arbitrary.

In summary: BH, Fisher exact, the polygamma functions and OLS agree with scipy/statsmodels/numpy to
around 1e-13 or better; `squeezeVar` and everything downstream of it agrees with inmoose to 1e-9,
which is the precision limma's own Newton solve for the prior degrees of freedom delivers.

## PRISM and `proteomics-toolkit` are not interchangeable

The lab's [`proteomics-toolkit`](https://github.com/uw-maccosslab/proteomics-toolkit) implements the
same idea in `run_moderated_linear_model`. It is a **third** implementation, not a second copy of
this one, and the two agree in one mode and not in the other. Both were run on the same inputs (the
`moderated_t.json` fixture cases) against the same limma reference:

| toolkit `moderation` | agreement with limma / PRISM |
|---|---|
| `"limma"` (global prior) | **~1e-13 on logFC, t and P.Value** - the three implementations are the same estimator |
| `"intensity_trend"` (the toolkit's **default**) | **median 3-7% on P.Value, up to 178% on individual features**; hit lists at p < 0.05 differed by 1-2 features out of 6-8 |

The global-prior agreement is the reassuring half: three independent implementations of
Smyth (2004) landing on the same numbers is good evidence all three are right.

The trend disagreement is not a bug in either tool. They are **different estimators that share a
name**:

- **limma `trend=TRUE`, which PRISM implements** (`EmpiricalBayes.SqueezeVarTrend`): fitFDist with
  a covariate. A natural cubic spline of `log(s^2)` on mean **log2** expression, with the spline df
  chosen as `1 + (n>=3) + (n>=6) + (n>=30)`. The prior scale **and** the prior degrees of freedom
  are both re-estimated from that spline fit.
- **toolkit `moderation="intensity_trend"`** (`_fit_intensity_trend_prior`): a LOWESS of
  `log(within-group variance)` on `log(within-group mean)` computed on **raw, pre-log** intensities,
  one point per (feature, group); converted back to log space by the delta method
  (`var_log ~ var_raw / mean_raw^2`) and combined across groups as a sample-size-weighted mean. The
  prior degrees of freedom stay at the **global** `d0` - only the per-feature scale is replaced.

So they differ in the smoother (spline vs LOWESS), in the space the trend is fit in (log2
abundance vs raw intensity), in what contributes a point (a feature vs a feature-group pair), and
in whether the prior df is re-estimated. Expecting them to agree to more than a few percent would
be expecting a coincidence.

Two practical consequences:

- **Do not treat a result from one as reproducing a result from the other** when the trend prior is
  on. Quote which tool and which moderation produced a hit list.
- **`moderation="limma"` is the setting to compare across the two tools.** If a comparison is the
  point, use it in both.

The toolkit's own global prior also differs from limma's `fitFDist` in two small ways that do not
matter on dense proteomics data and could on sparse: it drops zero and negative residual variances
from the prior fit where limma floors them at `1e-5 * median` and keeps them, and it has no
single-feature branch (limma returns `df_prior = 0` there).

## Cross-checking a PRISM result against the toolkit

```python
cfg.moderation = "limma"   # not the default; the comparable setting
```

and read the PRISM side from the same matrix PRISM reported on
(`corrected_proteins.parquet` / `corrected_peptides.parquet`, which are **linear** - log2 them
first; see the scale conventions in `CLAUDE.md`).

## Detection (peptide on/off)

The Detection view tests whether a peptide is observed at different rates between the groups. With
no covariates that is Fisher's exact test, pinned to `scipy.stats.fisher_exact`. With covariates it
is a **Firth-penalized** likelihood-ratio test - the full design against the same design with the
group column removed, referred to chi-square on 1 df.

Firth rather than an ordinary logistic regression because the interesting case is exactly the one
an unpenalized fit cannot express: a peptide detected in every sample of one group and none of the
other has no finite maximum-likelihood estimate. No standard library implements Firth, so the
golden reference maximizes the same penalized log-likelihood with a derivative-free optimizer
instead - a different algorithm reaching the same fixed point.
