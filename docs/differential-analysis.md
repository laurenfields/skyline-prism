# Differential analysis

The Skyline tool's **Differential** pane runs a limma empirical-Bayes moderated-t contrast on a
PRISM output directory: protein or peptide level, two groups from any metadata or clinical column,
optional covariate adjustment, and a choice of empirical-Bayes variance prior. The statistics live
in `SkylinePrism.Core.DifferentialAnalysis` and are reachable without the GUI.

The **default prior is the lab's own**, `VariancePriors.IntensityTrend` - the same estimator as
`proteomics-toolkit`'s `moderation="intensity_trend"`, and pinned to it by
`intensity_trend.json`. limma's `trend=TRUE` is still available as **limma-trend**, and the two are
not the same thing; the section below is about exactly that.

## What each quantity is held to

Every component with a published reference implementation is pinned to it by committed goldens in
`dotnet/tests/fixtures/differential/`, regenerable with
`uv run dotnet/tests/fixtures/differential/generate.py`. See that directory's README for the
per-quantity reference table and for the tolerances, which are not uniform and are not arbitrary.

In summary: BH, Fisher exact, the polygamma functions and OLS agree with scipy/statsmodels/numpy to
around 1e-13 or better; `squeezeVar` and everything downstream of it agrees with inmoose to 1e-9,
which is the precision limma's own Newton solve for the prior degrees of freedom delivers.

## What the pane can run

Three orthogonal choices - a **design**, a **test**, and (for the moderated t alone) a **variance
prior** - which is how `proteomics-toolkit` models it too.

| Design | Test | Reference |
|---|---|---|
| unpaired | moderated t | lstsq + `squeezeVar`, composed as limma composes them |
| unpaired | Welch t / Student t | `scipy.stats.ttest_ind(equal_var=False/True)` |
| unpaired | Mann-Whitney U | `scipy.stats.mannwhitneyu(method="asymptotic")` |
| paired | moderated t | the same composition over `[1, group, subject dummies]` |
| paired | paired t | `scipy.stats.ttest_rel` |
| paired | Wilcoxon signed-rank | `scipy.stats.wilcoxon(method="asymptotic")` |

Detection (the peptide on/off view) follows the design too: Fisher exact unpaired, McNemar's exact
test paired (`statsmodels.stats.contingency_tables.mcnemar(exact=True)`), and the Firth-penalized LRT
when covariates are set — the last being unpaired whatever the design, since a paired adjusted binary
outcome needs conditional logistic regression, which is not implemented.

Variance priors: **global** (Smyth 2004), **intensity trend** (the default - see below), **limma-trend**
(`trend=TRUE`), and **peptide count** (DEqMS, protein level only, pinned to the toolkit's
`_fit_count_dependent_prior`). Any prior that fits a per-feature scale leaves the prior *degrees of
freedom* at the global value; only limma-trend re-estimates both. The prior can optionally be fitted on
the run's QC and reference replicates instead of on the contrast groups.

Multiple-testing corrections: Benjamini-Hochberg, Benjamini-Yekutieli, Holm, Bonferroni, none - all
pinned to `statsmodels multipletests`, and all keeping PRISM's own NaN policy (a NaN passes through and
does not count toward *m*, where statsmodels returns all-NaN).

### Two things the rank tests do not do

**Neither implements an exact permutation distribution.** scipy's `method="auto"` uses one below n = 9
for Mann-Whitney and up to n = 50 for Wilcoxon; PRISM always uses the normal approximation with a tie
correction, so on a small group its p-value differs from a default-argument scipy run. The goldens pin
`method="asymptotic"` rather than encode a branch PRISM does not have.

**scipy's two rank tests disagree about the continuity correction** - `mannwhitneyu` applies one by
default, `wilcoxon` does not - and PRISM follows each one's own convention. Copying either to the other
shifts every p-value in that test by 10-20%, which is how this was found.

### Paired is a fixed-effect subject block

`[1, group, subject dummies]`, matching `statistical_analysis.py:1321`. It is not a random intercept
and not a mixed model; those are a different estimator and are not implemented. What the block buys is
that each subject's overall level leaves the residual, so a within-subject shift is tested against
within-subject noise instead of against the spread between people.

## From the command line

The same analysis, against the same output directory, without Skyline or Windows:

```bash
prism differential -d output/ --group-by condition -a Control -b Disease
```

`prism differential` takes the whole menu above as flags - `--level`, `--design`, `--pair-by`,
`--test`, `--prior`, `--prior-from-controls`, `--adjust-for`, `--correction` - and writes a results
CSV (`differential.csv` in the output directory unless `-o` says otherwise) whose header records the
contrast, its direction and the method that produced it. `prism differential --help` lists every
flag with its default.

Each arm takes several levels, and the arm is their **union**, so a three-level column can be
collapsed into a two-group contrast in one command:

```bash
prism differential -d output/ --group-by stage -a Control Mild -b Severe --adjust-for sex,age
```

A level named on both sides is refused rather than dropped from one, because which side it was
dropped from would change the answer and nothing in the output would record the choice.

The pane and the command resolve their arms through the same `ContrastArms` in Core, and run the
same `Differential.Run`, so a contrast set up by clicking and one typed out mean the same samples
and give the same numbers - checked on a 192-sample cohort, where the two agree bit for bit on
log2FC, p and adjusted p.

One caveat worth knowing before comparing two runs of your own: **which arm is A and which is B is
not a pure sign flip.** Swapping them reverses the sign of log2FC as you would expect, but it also
changes the order values are accumulated in, so the magnitudes move in the last digit or two
(around 1e-13 relative on that cohort). That is floating-point summation order, not a difference in
what was computed - but it means two results are only comparable digit-for-digit if the arms were
given the same way round.

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

- **limma `trend=TRUE`** (`EmpiricalBayes.SqueezeVarTrend`, offered as **limma-trend**): fitFDist with
  a covariate. A natural cubic spline of `log(s^2)` on mean **log2** expression, with the spline df
  chosen as `1 + (n>=3) + (n>=6) + (n>=30)`. The prior scale **and** the prior degrees of freedom
  are both re-estimated from that spline fit.
- **toolkit `moderation="intensity_trend"`** (`_fit_intensity_trend_prior`; PRISM's
  `VariancePriors.IntensityTrend`, and the **default**): a LOWESS of
  `log(within-group variance)` on `log(within-group mean)` computed on **raw, pre-log** intensities,
  one point per (feature, group); converted back to log space by the delta method
  (`var_log ~ var_raw / mean_raw^2`) and combined across groups as a sample-size-weighted mean. The
  prior degrees of freedom stay at the **global** `d0` - only the per-feature scale is replaced.

So they differ in the smoother (spline vs LOWESS), in the space the trend is fit in (log2
abundance vs raw intensity), in what contributes a point (a feature vs a feature-group pair), and
in whether the prior df is re-estimated. Expecting them to agree to more than a few percent would
be expecting a coincidence.

PRISM now implements both, so the disagreement is selectable rather than baked in - but it is still
a disagreement, and the two remain different estimators that happen to share a name. Two practical
consequences:

- **Do not treat a result from one as reproducing a result from the other** when the trend prior is
  on. Quote which tool and which moderation produced a hit list.
- **`moderation="limma"` is the setting to compare across the two tools.** If a comparison is the
  point, use it in both - PRISM's equivalent is the **Global** prior.
- **To reproduce a toolkit run at its own default**, leave PRISM on **Intensity trend**: that is the
  same estimator, held to the toolkit by a committed golden.

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
