# Skyline-PRISM vNEXT Release Notes

Working draft. Rename to `RELEASE_NOTES_v{version}.md` at release time.

## New Features

<!-- none yet -->

## Bug Fixes

### Fixed: `BrokenProcessPool` on long parallel rollup runs (Linux fork+COW OOM)

After ~45 minutes / tens of thousands of peptides processed, the chunked transition rollup could fail with `concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly`. Most often observed with the `consensus` rollup method on cohorts of 200+ samples.

Root cause: `ProcessPoolExecutor` was created and torn down inside the per-batch loop (every ~1000 peptides) with the default `fork` start method on Linux. Each fork inherited the parent's accumulated `peptide_rows` / `residual_rows` / `consensus_diag_rows` via copy-on-write. As the parent's result state grew (GBs after tens of thousands of peptides), forking workers from that snapshot pushed total memory over the OS limit and the kernel OOM-killer terminated a worker mid-task, surfacing as `BrokenProcessPool`.

Fix: All three `ProcessPoolExecutor(...)` sites in `skyline_prism/chunked_processing.py` now pass `mp_context=mp.get_context("spawn")`. Spawned workers start as fresh Python interpreters with no inherited parent memory, so each worker only carries the chunk it is actually processing. Per-pool startup cost is ~150 ms; on a multi-hour run this adds up to a few seconds — negligible compared to the per-peptide rollup work.

- **Files modified**: `skyline_prism/chunked_processing.py`

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
