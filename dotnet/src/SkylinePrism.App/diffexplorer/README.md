# Bundled Differential Explorer payload

This folder is the payload launched by the **"Differential Analysis →"** button
(`OnOpenDifferential` in `MainWindow.xaml.cs`). It is copied next to the exe at
build time and staged into the tool `.zip`.

Contents:

- `prism_diff_explorer.py`, `cryptic_gene_map.csv`, `requirements.txt` — **vendored**
  from https://github.com/laurenfields/prism-diff-explorer (single-module Streamlit
  app). Current copy is from commit **`eafba2f`**.
- `uv.exe` — the [uv](https://docs.astral.sh/uv/) launcher (Astral). **Not committed**
  (`.gitignore`d — ~75 MB); the release/packaging step fetches it into this folder
  before zipping, and a developer copies it here for local runs.

The button runs, with `PDX_PRISM_DIR=<PRISM output dir>` in the environment:

```
uv run --python 3.12 --with-requirements requirements.txt --no-project \
  streamlit run prism_diff_explorer.py --server.headless=true --server.port=<free> --server.fileWatcherType=none
```

`uv` provisions an isolated Python 3.12 + the pinned deps on first use (no system
Python / conda needed); first launch ~1 min, cached launches ~seconds.

## Updating the vendored app

Re-copy the three files from the prism-diff-explorer repo and update the commit
hash above:

```
cp ../../../../../prism-diff-explorer/{prism_diff_explorer.py,cryptic_gene_map.csv,requirements.txt} .
```
