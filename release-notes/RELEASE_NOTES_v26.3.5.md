# Skyline-PRISM v26.3.5 Release Notes

Patch release: correctness fix for Tukey median polish (matches R `stats::medpolish` and Tukey 1977 reference; thanks to Nick Shulman). Also fixes long-rollup parent-process OOM via streaming output writes, eliminates `BrokenProcessPool` from fork+COW worker OOMs via the spawn start method, suppresses expected `All-NaN slice encountered` warnings in library-assist rollup, and reverts the default `transition_rollup.method` to `sum`.

## New Features

### Default `transition_rollup.method` reverted to `sum`

Both `prism config-template` and `prism config-template --minimal` now emit `transition_rollup.method: "sum"` as the default. The v26.3.3 release had switched the default to `library_assist`, which required users to also provide a spectral library before they could run anything; `sum` is a usable starting point with no extra files needed.

The `library_assist` method is still fully supported. The minimal template's section header documents how to switch (replace the `transition_rollup` block with the four-line `library_assist` recipe and set `library_path`); the full template's method comment now lists `library_assist` with a note that `library_assist.library_path` must be set if you choose it. The "BEFORE RUNNING - edit two paths" header callout in both templates is now a single-path callout: only `parsimony.fasta_path` is required up front.

- **Files modified**: `skyline_prism/cli.py` (template functions `get_minimal_config_template` and `get_full_config_template`)

## Bug Fixes

### Fixed: Tukey median polish recentering (correctness fix, thanks to Nick Shulman, [PR #5](https://github.com/maccoss/skyline-prism/pull/5))

`tukey_median_polish` (`skyline_prism/rollup.py`) now matches the canonical Tukey 1977 / R `stats::medpolish` definition. The previous implementation only re-centered each iteration's per-row update (`row_effect_update = row_medians - median(row_medians)`), but a sum of median-zero vectors does not have median zero in general, so the cumulative `row_effects` and `col_effects` drifted from 0 across iterations. Without the Tukey normalization, the partition between `overall`, `row_effects`, and `col_effects` is ambiguous, and `overall + col_effects[j]` (the per-sample output PRISM consumes) was biased by an unspecified per-matrix constant equal to `median(row_effects)` at convergence. On a real Skyline test dataset (33 peptides × 3 samples), this bias produced log2 offsets up to 0.624 (~1.54x linear) compared to R's medpolish on the same input. After the fix the worst-case offset drops to 0.054 log2 (3.7% linear), entirely from convergence/tie-breaking noise rather than the partition bias.

The fix re-centers `row_effects` and `col_effects` against their own running medians at each iteration and transfers the offset into `overall`, exactly as R's reference implementation does.

**Practical impact on existing PRISM analyses** (the bias is a *constant* additive shift across all samples within one decomposition, with a *different* constant per protein/peptide):

- Sample-vs-sample log2 fold change per protein → **unaffected** (the bias cancels)
- Linear-scale CV across replicates per protein → **unaffected** (multiplicative `2^m_P` factors out)
- Differential expression / t-tests → **unaffected**
- Cross-protein absolute abundance, iBAQ, absolute quantification → **were biased**, now correct

The bug applied to:

- **peptide → protein rollup** with `protein_rollup.method: "median_polish"` (in `rollup.rollup_protein_matrix` and `rollup.rollup_to_proteins`)
- **transition → peptide rollup** with `transition_rollup.method: "median_polish"` (in `chunked_processing._process_single_peptide`)

Both sites consume `result.col_effects` directly and inherit Nick's fix automatically.

The **library-assisted rollup** (`spectral_library.library_median_polish_rollup_vectorized`) is a different algorithm and was **not** affected: it fixes row effects to `log(library)` from the spectral library prior, computes per-sample `beta_s` via a single direct `np.nanmedian`, and never iterates a row+column sweep whose accumulators could drift.

Nick's PR #5 added two regression tests in `tests/test_rollup.py`:

- **`test_tukey_normalization_at_convergence`** — asserts `median(row_effects)` and `median(col_effects - overall)` are both ~0 at convergence.
- **`test_matches_r_medpolish_reference`** — compares `overall` and per-sample `col_effects` against R `medpolish` output captured for a 5×3 matrix.

This release adds three additional regression tests in the same class:

- **`test_within_protein_fold_changes_match_r`** — locks in the practical invariant that proteomics differential-expression analyses depend on: sample-vs-sample log2 fold changes within one protein match R medpolish on the same input.
- **`test_decomposition_reconstructs_input_at_convergence`** — structural correctness: the additive model must reconstruct the input (`X[i,j] == row_effects[i] + col_effects[j] + residuals[i,j]`, where PRISM's `col_effects` already includes `overall`).
- **`test_handles_nan_values_in_matrix`** — recentering still produces a valid decomposition with the Tukey normalization holding (~0 effect medians) when the matrix has missing observations, and the reconstruction is exact at non-NaN cells.

- **Files modified**: `skyline_prism/rollup.py`, `tests/test_rollup.py`

### Fixed: parent-process OOM (`Killed`) on long rollup runs from in-memory result accumulation

After ~45,000 peptides processed (typically several hours into a large run, most reproducibly with the `consensus` rollup method), the parent `prism` process was OOM-killed by the kernel with no Python traceback (just `Killed` in the terminal). The previous `BrokenProcessPool` fix (spawn workers below) had eliminated worker-side OOMs, but the parent process itself was still accumulating the entire run's results in three Python lists before writing them at the end:

- `peptide_rows`: ~1 GB by 45k peptides × 250 samples
- `residual_rows` (when `save_residuals` is on): ~250 MB
- `consensus_diag_rows` (one entry per peptide × transition × sample for the `consensus` method with diagnostics): ~6+ GB by 45k peptides; ~10 GB by run end on 70k peptides × 8 transitions × 250 samples = 140M rows

Fix: stream all three accumulators to disk during the run. `_flush_buffers()` is called after each parallel batch (every ~1000 peptides) and after every `FLUSH_THRESHOLD` (2000) peptides in the single-threaded path:

- `peptide_rows` and `residual_rows` are flushed to their final parquet via per-batch `pq.ParquetWriter.write_table()` calls. The parquet schemas are pre-defined from the full sample list so per-batch DataFrames can be reindexed and cast consistently across writes (handles peptides that are missing some samples).
- `consensus_diag_rows` is streamed to an intermediate parquet during the run. At the end, the intermediate is sorted by `abs_residual DESC` via DuckDB (with a pandas fallback) and written to the user-configured CSV path. The intermediate parquet is then deleted. This avoids holding 140M rows in memory while preserving the "most problematic first" CSV output ordering.

Memory now stays bounded by `FLUSH_THRESHOLD × per-row size` (~hundreds of MB) regardless of total peptide count, so multi-hour 70k-peptide consensus runs no longer OOM the parent.

- **Files modified**: `skyline_prism/chunked_processing.py`

### Fixed: `BrokenProcessPool` on long parallel rollup runs (Linux fork+COW OOM)

After ~45 minutes / tens of thousands of peptides processed, the chunked transition rollup could fail with `concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly`. Most often observed with the `consensus` rollup method on cohorts of 200+ samples.

Root cause: `ProcessPoolExecutor` was created and torn down inside the per-batch loop (every ~1000 peptides) with the default `fork` start method on Linux. Each fork inherited the parent's accumulated `peptide_rows` / `residual_rows` / `consensus_diag_rows` via copy-on-write. As the parent's result state grew (GBs after tens of thousands of peptides), forking workers from that snapshot pushed total memory over the OS limit and the kernel OOM-killer terminated a worker mid-task, surfacing as `BrokenProcessPool`.

Fix: All three `ProcessPoolExecutor(...)` sites in `skyline_prism/chunked_processing.py` now pass `mp_context=mp.get_context("spawn")`. Spawned workers start as fresh Python interpreters with no inherited parent memory, so each worker only carries the chunk it is actually processing. Per-pool startup cost is ~150 ms; on a multi-hour run this adds up to a few seconds — negligible compared to the per-peptide rollup work.

- **Files modified**: `skyline_prism/chunked_processing.py`

### Suppressed expected `All-NaN slice encountered` warnings in library-assisted rollup

Library-assisted rollup (`library_median_polish_rollup_vectorized` in `skyline_prism/spectral_library.py`) called `np.nanmedian` per-sample to compute the per-sample log-scale `beta_s`. For any sample whose fragments are all NaN — a legitimate state when every fragment for that sample was missing or excluded by outlier filtering — numpy emits `RuntimeWarning: All-NaN slice encountered` and returns NaN. The NaN is the intended sentinel and is handled correctly downstream, so the warnings were pure noise and would mask other diagnostic output during test runs.

Fix: Both `np.nanmedian` call sites (lines ~929 and ~985) are now wrapped in a narrowly-scoped `warnings.catch_warnings()` + `filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)`. The filter is targeted: only this specific warning at these specific calls is suppressed; any other `RuntimeWarning` from numpy would still surface.

Test suite now reports zero warnings (was: `342 passed, 2 warnings`).

- **Files modified**: `skyline_prism/spectral_library.py`

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
