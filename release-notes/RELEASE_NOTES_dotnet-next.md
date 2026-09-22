# Skyline-PRISM (C#) dotnet-vNEXT Release Notes

Working draft for the next C# (.NET) release. Append entries as they land on the development branch;
rename to `RELEASE_NOTES_dotnet-v{version}.md` at release time - the release workflow publishes this file
as the GitHub Release description and fails if it is missing.

## New Features

- **Differential analysis pane.** A new "Differential" visualization pane runs a limma-style moderated-t
  contrast over the corrected matrix (protein or peptide level) entirely locally. Pick a grouping column
  and two values (A vs B), optionally adjust for covariates, and toggle the intensity-trend variance
  prior (limma-trend). Views: a **Volcano** plot with a ranked hit table, a sample **PCA**, a peptide
  **Detection**-frequency test (Fisher exact, or Firth-penalized logistic regression when covariates are
  set), and functional **Enrichment** of the significant hits via g:Profiler. The statistics live in
  `SkylinePrism.Core.DifferentialAnalysis` and are validated to 1e-9 against the reference implementation.
- **Per-feature boxplots.** Clicking a point on the Volcano plot - or a row in the hit table - opens a
  detail window with that feature's log2 abundance split into a boxplot per contrast group, with jittered
  per-sample points and the log2FC / adj.P in the title.
- **Attach an external clinical CSV.** A top-level "Clinical CSV" input (beside the metadata report)
  joins an external clinical metadata table to the samples by an auto-detected key column (best value
  match, preferring near one-to-one keys). Joined columns become available for grouping and covariate
  adjustment across the analysis panes.
- **Markers pane.** A new "Markers" visualization pane evaluates a protein panel against the corrected
  matrix: a row z-scored log2 heatmap of the panel's members across a chosen grouping column (group means,
  toggle to per-sample) plus a per-group boxplot of each sample's mean marker z-score, and a "found
  N/total members / not detected" report. Panels are the same protein lists the Dynamic Range plot uses
  (your own plus the shipped panels), and several can be ticked to union into one heatmap.

- **One PCA, not two.** The QC scatter and the Differential pane's sample plot were separate
  implementations; they are now one `Core/Numerics/Pca.cs` with options for what actually differed
  - standardize vs center-only, impute-to-feature-mean vs complete-case, two components vs k with
  their variance ratios, and all samples vs a chosen subset. The merged engine keeps the QC path's
  blocked Gram accumulation, so the differential configuration no longer runs an unblocked
  O(samples^2) loop per feature.
- **The QC PCA can be grouped by the run's own metadata and by a clinical CSV.** Its "Group by"
  list used to come only from an exported Skyline Replicates report, so a run produced by the CLI
  offered nothing but Sample Type - and a batch-correction tool could not colour its PCA by batch.
  It now also reads the run's `sample_metadata.csv` (which supplies `batch`) and any clinical CSV
  attached in the Differential pane. The QC PCA also gained a component-pair selector up to PC6,
  with each component's variance-explained on its axis.
- **The Differential pane's PCA sub-view is gone**, now that the QC pane's PCA can group by the same
  columns. It existed largely to work around the gap above. Note the surviving plot uses the QC
  pane's long-standing settings - features standardized, a missing cell imputed to its feature mean
  - where the removed one centered without scaling and dropped any feature with a gap, so the two
  do not draw identical point positions on a matrix with missing values.
- **Committed goldens for the differential statistics.** `dotnet/tests/fixtures/differential/`
  holds reference values generated from scipy, statsmodels and inmoose by a checked-in script, and
  `DifferentialGoldenTests` holds PRISM to them: BH, the polygamma functions, `lmFit`, `squeezeVar`
  (global and intensity-trend), the natural-spline basis, the end-to-end moderated t, Fisher exact,
  Firth, the penalized detection LRT and the sample PCA. Unlike the other fixture directories these
  are regenerable - `uv run dotnet/tests/fixtures/differential/generate.py`.
- **`docs/differential-analysis.md`** records what each quantity is held to, and where PRISM and the
  lab's `proteomics-toolkit` agree (`moderation="limma"`, to ~1e-13) and do not
  (`moderation="intensity_trend"`, a different estimator: median 3-7% on p-values, and hit lists
  that differ).

## Bug Fixes

- The QC PCA no longer materializes a transposed copy of the abundance matrix before fitting - a
  full second copy of the largest object in the pipeline (5.7 GB on a 100-document peptide matrix,
  on the large object heap), built only to be read one feature at a time, which is the layout it
  started in. It now calls the `Fit2DOfFeaturesBySamples` overload that exists to avoid exactly that.
- The differential loader now aligns a run's sample columns to `sample_metadata.csv` by the bare replicate
  name (the part before `__@__`) when no column matches by the full `sample_id`, so a run whose corrected
  matrix and metadata were written with different document/batch stems still loads.

## Performance

## Breaking Changes
