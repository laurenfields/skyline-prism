"""
PRISM Differential Explorer
===========================

Interactive differential-abundance analysis for PRISM output directories.

Load a PRISM `output_dir/` (the folder containing `corrected_proteins.parquet`,
`corrected_peptides.parquet`, and `sample_metadata.csv`), optionally join a
clinical-metadata CSV, define two groups (by a metadata column or by hand), and
run a limma-style empirical-Bayes moderated t-test with optional covariates. Get
an interactive volcano plot, a filterable table of hits, and per-feature plots.

Statistics
----------
Honours the project methodology rules: limma's empirical-Bayes *moderated*
t-test (variance shrinkage via ``inmoose.limma.squeezeVar``) with
Benjamini-Hochberg FDR -- NOT a per-feature Welch t-test. Arbitrary covariates
(batch, cohort/Study_Name, Sex, Age, PMI, ...) can be added to the linear model;
the app reports the group x covariate balance so you can see whether a contrast
is confounded before trusting it.

Caveat surfaced in the UI: ``corrected_*.parquet`` is dense. The test is on
integrated-window intensity, which conflates detection with abundance; cells for
truly-absent features are integrated baseline noise (real but low-information).
For genuine on/off detection, recover it from transition-level ``merged_data``.

Run
---
    conda activate lf01
    streamlit run prism_diff_explorer.py

Core functions (`load_prism`, `attach_clinical`, `differential`) are
Streamlit-free and reusable from a notebook.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

import inmoose.limma as L

# --------------------------------------------------------------------------- #
# Core (Streamlit-free) logic
# --------------------------------------------------------------------------- #

_PROTEIN_ID = "protein_group"
_PROTEIN_LABEL = "leading_gene_name"
_PEPTIDE_ID = "PeptideModifiedSequenceUnimodIds"

# Canonical extracellular-vesicle markers (gene symbols), grouped by MISEV-style
# category. Positive markers should be enriched in EV preparations; the
# "Contamination" set should be depleted (co-isolated plasma / ER / nuclear).
EV_MARKERS: dict[str, list[str]] = {
    "Tetraspanins": ["CD9", "CD63", "CD81", "CD82", "CD37", "CD53"],
    "ESCRT / biogenesis": ["TSG101", "PDCD6IP", "SDCBP", "CHMP4B", "VPS4A", "VPS4B",
                            "ARRDC1", "MVB12A"],
    "Membrane / receptors": ["CD47", "BSG", "ITGB1", "ITGA4", "SLC3A2", "TFRC",
                             "LAMP1", "LAMP2", "ATP1A1", "MFGE8", "PTGFRN"],
    "Flotillins": ["FLOT1", "FLOT2"],
    "Rab / trafficking": ["RAB1A", "RAB5A", "RAB5B", "RAB5C", "RAB7A", "RAB11A",
                          "RAB11B", "RAB27A", "RAB27B", "RAB35", "ARF6", "EHD1", "EHD4"],
    "Cytosolic / chaperones": ["HSPA8", "HSP90AA1", "HSP90AB1", "ANXA1", "ANXA2",
                               "ANXA5", "ANXA6", "ANXA11", "ACTB", "GAPDH", "EEF1A1",
                               "YWHAZ", "YWHAE"],
}
EV_NEGATIVE: list[str] = [
    "ALB", "APOA1", "APOA2", "APOB", "APOE", "FGA", "FGB", "FGG", "SERPINA1",
    "TF", "HP", "CALR", "CANX", "HSPA5", "PDIA3",
]

# --- Local group suggestion (no API, no cost) -----------------------------
# Heuristically proposes two-group comparisons from the metadata by classifying
# each categorical value as control-like vs case-like via keyword matching, then
# auto-picking nuisance covariates and flagging imbalanced ones.

_CONTROL_RE = re.compile(
    r"\b(no|non|without|absence|neg|negative|control|ctrl|ctl|con|normal|healthy|hc|co|cn|nc|"
    r"reference|ref|wt|wild|baseline|untreated|vehicle|mock|sham|low|unaffected|none|benign)\b")
_CASE_RE = re.compile(
    r"\b(disease|patient|case|dementia|ad|add|adnc|tumor|tumour|cancer|malignant|mutant|"
    r"mutation|positive|pos|affected|sporadic|autosomal|psen|carrier|high|treated|stimulated|"
    r"infected|severe|stage|pd|ftd|als|dlb|psp|cbd|mci|lbd|tdp|park)\b")
# column names that are nuisance/technical rather than the biological contrast
_COVARIATE_NAME_RE = re.compile(r"sex|gender|study|cohort|site|batch|plate|apoe|region|age|pmi", re.I)


def _classify_level(value: str) -> str | None:
    """Label a categorical value 'control', 'case', or None (control wins ties)."""
    v = str(value).lower()
    if _CONTROL_RE.search(v):
        return "control"
    if _CASE_RE.search(v):
        return "case"
    return None


def suggest_groups(meta: pd.DataFrame, max_suggestions: int = 6,
                   min_group_n: int = 2) -> list[dict]:
    """Local, deterministic comparison suggester (no external service).

    Returns suggestion dicts with the same shape the UI consumes:
    name, grouping_column, group_a/b_values, label_a/b, suggested_covariates,
    rationale, confound_warning.
    """
    exp = meta[meta["sample_type"] == "experimental"] if "sample_type" in meta.columns else meta
    skip = {"sample_id", "sample", "replicate", "sample_type"}
    cat_cols = [c for c in meta.columns if c not in skip
                and not pd.api.types.is_numeric_dtype(meta[c]) and exp[c].notna().any()]
    num_cols = [c for c in meta.columns if c not in skip and pd.api.types.is_numeric_dtype(meta[c])]

    out: list[dict] = []
    for col in cat_cols:
        counts = exp[col].dropna().astype(str).value_counts()  # descending
        levels = [lv for lv in counts.index if counts[lv] >= min_group_n]
        if not (2 <= len(levels) <= 10):
            continue
        cls = {lv: _classify_level(lv) for lv in levels}
        controls = [lv for lv in levels if cls[lv] == "control"]
        cases = [lv for lv in levels if cls[lv] == "case"]
        if controls and cases:
            a_vals, b_vals, priority = controls, cases, 2.0
        elif len(levels) == 2:
            a_vals, b_vals, priority = [levels[0]], [levels[1]], 1.0  # larger as reference
        else:
            continue
        na, nb = int(counts[a_vals].sum()), int(counts[b_vals].sum())
        if na < min_group_n or nb < min_group_n:
            continue

        la = "Control" if len(a_vals) > 1 else str(a_vals[0])
        lb = "Case" if len(b_vals) > 1 else str(b_vals[0])

        # Only genuine nuisance variables (by name) — NOT disease-severity scores
        # like Braak/CERAD, which are part of the phenotype (colliders), nor other
        # numerics that may be downstream of the disease.
        cov = [c for c in num_cols if _COVARIATE_NAME_RE.search(c)]
        cov += [c for c in cat_cols if c != col and _COVARIATE_NAME_RE.search(c)]
        cov = list(dict.fromkeys(cov))

        # confound check across the proposed groups
        grp = {**{lv: "A" for lv in a_vals}, **{lv: "B" for lv in b_vals}}
        sub = exp[exp[col].astype(str).isin(a_vals + b_vals)].copy()
        sub["_g"] = sub[col].astype(str).map(grp)
        warn = []
        for cv in cov:
            if cv in num_cols:
                med = sub.groupby("_g")[cv].median()
                sd = sub[cv].dropna().std()
                if {"A", "B"} <= set(med.index) and sd and sd > 0:
                    d = abs(med["A"] - med["B"]) / sd
                    if d > 0.8:
                        warn.append(f"{cv} differs ~{d:.1f} SD between groups")
            else:
                ct = pd.crosstab(sub[cv].astype(str), sub["_g"])
                if ct.shape[0] > 1 and ((ct == 0).any(axis=1) & (ct.sum(axis=1) > 0)).any():
                    warn.append(f"{cv} partly confounded with group")

        out.append({
            "name": f"{lb} vs {la} ({col})",
            "grouping_column": col,
            "group_a_values": list(a_vals),
            "group_b_values": list(b_vals),
            "label_a": la, "label_b": lb,
            "suggested_covariates": cov,
            "rationale": f"{col}: {nb} {lb} vs {na} {la} samples"
                         + (" — add the flagged covariates." if warn else "."),
            "confound_warning": "; ".join(warn),
            "_priority": priority + min(na, nb) / 1000.0,
        })

    out.sort(key=lambda s: s["_priority"], reverse=True)
    for s in out:
        s.pop("_priority", None)
    return out[:max_suggestions]


# --- Optional local-LLM backend (Ollama) ----------------------------------
# Same suggestion shape as suggest_groups(), but powered by a local model via
# Ollama's HTTP API. Useful for messy/free-text metadata the heuristic can't
# parse. Fully local: no API key, no cloud. Falls back to the heuristic in the
# UI if the server isn't reachable.

_OLLAMA_HOST = "http://localhost:11434"

_OLLAMA_SYSTEM = (
    "You are a mass-spectrometry proteomics analysis assistant. Given a JSON profile of "
    "sample metadata, propose biologically meaningful two-group comparisons for "
    "differential-abundance testing. Use ONLY the experimental (biological) samples; never "
    "use QC or reference pools. For each comparison pick ONE categorical column, assign a "
    "subset of its EXACT listed values to group A (control/reference) and a disjoint subset "
    "to group B (case/treatment). Suggest nuisance covariates from OTHER columns that could "
    "confound the contrast (age, sex, cohort/study, post-mortem interval, batch) - but NOT "
    "disease-severity scores, which are part of the phenotype. Give a one-line rationale and "
    "a brief confounding caveat. Respond ONLY with JSON matching the schema."
)

# JSON schema for Ollama structured output (mirrors the suggestion dict).
_OLLAMA_SCHEMA = {
    "type": "object",
    "properties": {
        "suggestions": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "grouping_column": {"type": "string"},
                    "group_a_values": {"type": "array", "items": {"type": "string"}},
                    "group_b_values": {"type": "array", "items": {"type": "string"}},
                    "label_a": {"type": "string"},
                    "label_b": {"type": "string"},
                    "suggested_covariates": {"type": "array", "items": {"type": "string"}},
                    "rationale": {"type": "string"},
                    "confound_warning": {"type": "string"},
                },
                "required": ["name", "grouping_column", "group_a_values", "group_b_values",
                             "label_a", "label_b", "suggested_covariates", "rationale",
                             "confound_warning"],
            },
        }
    },
    "required": ["suggestions"],
}


def _profile_metadata(meta: pd.DataFrame, focus_type: str = "experimental",
                      max_levels: int = 20) -> dict:
    """Compact metadata profile (column names + value distributions) for the LLM."""
    prof: dict = {"focus_sample_type": focus_type}
    if "sample_type" in meta.columns:
        prof["sample_type_counts"] = {k: int(v) for k, v in meta["sample_type"].value_counts().items()}
        focus = meta[meta["sample_type"] == focus_type]
    else:
        focus = meta
    prof["n_focus_samples"] = int(len(focus))
    cols: dict = {}
    for c in meta.columns:
        if c in {"sample_id", "sample", "replicate"}:
            continue
        s = focus[c] if c in focus.columns else meta[c]
        if pd.api.types.is_numeric_dtype(s):
            sv = s.dropna()
            cols[c] = {"type": "numeric",
                       "min": float(sv.min()) if len(sv) else None,
                       "median": float(sv.median()) if len(sv) else None,
                       "max": float(sv.max()) if len(sv) else None}
        else:
            vc = s.dropna().astype(str).value_counts()
            cols[c] = {"type": "categorical", "n_levels": int(vc.shape[0]),
                       "values": {k: int(v) for k, v in vc.head(max_levels).items()}}
    prof["columns"] = cols
    return prof


def ollama_status(host: str = _OLLAMA_HOST, timeout: float = 2.0) -> tuple[bool, list[str]]:
    """Return (server_running, [model names]). Never raises."""
    try:
        import requests
        r = requests.get(host.rstrip("/") + "/api/tags", timeout=timeout)
        r.raise_for_status()
        return True, [m["name"] for m in r.json().get("models", [])]
    except Exception:  # noqa: BLE001
        return False, []


def suggest_groups_ollama(meta: pd.DataFrame, model: str, host: str = _OLLAMA_HOST,
                          timeout: float = 300.0) -> list[dict]:
    """Propose comparisons via a local Ollama model (structured JSON output).

    Slow on CPU (~1-2 min for a 3B model, first call also loads it); keep_alive
    keeps the model resident so repeat calls are faster.
    """
    import json

    import requests

    payload = {
        "model": model,
        "stream": False,
        "format": _OLLAMA_SCHEMA,
        "keep_alive": "10m",
        "options": {"temperature": 0, "num_predict": 1500},
        "messages": [
            {"role": "system", "content": _OLLAMA_SYSTEM},
            {"role": "user",
             "content": "Sample metadata profile (JSON):\n"
                        + json.dumps(_profile_metadata(meta), indent=2)},
        ],
    }
    r = requests.post(host.rstrip("/") + "/api/chat", json=payload, timeout=timeout)
    r.raise_for_status()
    content = r.json()["message"]["content"]
    return list(json.loads(content).get("suggestions", []))

# Columns never offered as grouping/covariate choices.
_NON_META_COLS = {"sample_id", "sample", "sample_type", "replicate"}

# A numeric metadata column with at most this many distinct values (e.g. a stage/grade/
# batch code) is also offered as a categorical "Group by" option, its values cast to
# strings. Columns with more distinct values are treated as continuous (covariate/time
# axis) only, so a genuinely continuous variable like Age isn't offered as 40 groups.
_MAX_NUMERIC_CAT_LEVELS = 25


def _as_cat_str(s):
    """String view of a metadata column for categorical grouping. Whole-number numeric
    codes (e.g. a 1.0/2.0 stage read as float because of NaNs) render as '1'/'2', not
    '1.0'. Missing values map to '' (numeric) or 'nan' (text) and are excluded from the
    value list, so they never form a group."""
    import pandas as pd
    if pd.api.types.is_numeric_dtype(s):
        nn = s.dropna()
        whole = len(nn) > 0 and bool((nn == nn.round()).all())
        return s.map(lambda v: "" if pd.isna(v)
                     else (str(int(round(v))) if whole else str(v)))
    return s.astype(str)


def _grouping_columns(meta, meta_cols):
    """Columns offered as 'Group by' categories: every text column, plus low-cardinality
    numeric columns (discrete codes like stage/grade/batch) whose values are treated as
    string categories. Continuous numerics are left out (covariate / time axis only)."""
    import pandas as pd
    out = []
    for c in meta_cols:
        if not pd.api.types.is_numeric_dtype(meta[c]):
            out.append(c)
        elif 2 <= meta[c].dropna().nunique() <= _MAX_NUMERIC_CAT_LEVELS:
            out.append(c)
    return out


def _bh(pvals):
    """Benjamini-Hochberg adjusted p-values for a 1-D array (NaNs preserved)."""
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    adj = np.empty(m)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


# Significance bases for the volcano / hit calls. 'fdr' = genome-wide BH (the adj.P.Val
# limma already computed across ALL tested features); 'family' = BH recomputed within the
# rows passed in (e.g. only cryptic features — the right multiplicity scope for a
# pre-specified cryptic study); 'raw' = uncorrected p-value (exploratory).
def _significance(frame, basis, q_cut, lfc_cut):
    """Return (significant: bool Series, p_threshold: float|None) for a results frame with
    columns P.Value, adj.P.Val, logFC. p_threshold is the raw-P value where the volcano's
    horizontal line goes (BH is monotone in raw P, so it's the largest passing raw P)."""
    import numpy as np
    import pandas as pd
    p = frame["P.Value"]
    if basis == "raw":
        adj = p
    elif basis == "family":
        adj = pd.Series(_bh(p.to_numpy()), index=frame.index)
    else:  # 'fdr' genome-wide
        adj = frame["adj.P.Val"]
    sig = (adj < q_cut) & (frame["logFC"].abs() >= lfc_cut)
    if basis == "raw":
        p_thresh = float(q_cut)
    else:
        passing = p[adj <= q_cut]
        p_thresh = float(passing.max()) if len(passing) else None
    return sig, p_thresh


def _add_volcano_guides(fig, p_thresh, lfc_cut, basis, q_cut):
    """Draw the +/-log2FC verticals and the horizontal significance line (at the raw-P
    threshold, since the y-axis is -log10 raw P) with a label naming the basis."""
    import numpy as np
    if lfc_cut > 0:
        fig.add_vline(x=lfc_cut, line_dash="dot", line_color="grey")
        fig.add_vline(x=-lfc_cut, line_dash="dot", line_color="grey")
    if p_thresh and p_thresh > 0:
        name = {"raw": f"p < {q_cut:g}", "family": f"cryptic-FDR {q_cut:g}",
                "fdr": f"FDR {q_cut:g}"}.get(basis, f"{q_cut:g}")
        extra = "" if basis == "raw" else f"  (p ≤ {p_thresh:.1e})"
        fig.add_hline(y=-np.log10(p_thresh), line_dash="dash", line_color="firebrick",
                      annotation_text=name + extra, annotation_position="top left")


def load_prism(output_dir: str | Path, level: str = "protein") -> dict:
    """Load a PRISM output directory for one feature level.

    Returns a dict: ``annot`` (annotation cols indexed by feature id),
    ``expr_log2`` (features x samples, log2), ``meta`` (sample_metadata indexed
    by sample_id), ``sample_cols``, ``id_col``, ``label_col``, ``level``.
    """
    output_dir = Path(output_dir)
    fname = "corrected_proteins.parquet" if level == "protein" else "corrected_peptides.parquet"
    fpath = output_dir / fname
    if not fpath.exists():
        raise FileNotFoundError(f"Not found: {fpath}")
    meta_path = output_dir / "sample_metadata.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Not found: {meta_path}")

    df = pd.read_parquet(fpath)
    meta = pd.read_csv(meta_path).set_index("sample_id", drop=False)

    sample_ids = set(meta["sample_id"])
    sample_cols = [c for c in df.columns if c in sample_ids]
    if not sample_cols:
        raise ValueError(
            "No columns in the parquet matched sample_metadata.sample_id. "
            "Are the output_dir and level consistent?"
        )
    annot_cols = [c for c in df.columns if c not in sample_ids]

    if level == "protein":
        id_col = _PROTEIN_ID if _PROTEIN_ID in df.columns else annot_cols[0]
        label_col = _PROTEIN_LABEL if _PROTEIN_LABEL in df.columns else id_col
    else:
        id_col = _PEPTIDE_ID if _PEPTIDE_ID in df.columns else annot_cols[0]
        label_col = id_col

    df = df.set_index(id_col, drop=False)
    annot = df[annot_cols].copy()

    expr = df[sample_cols].astype(float)  # corrected_*.parquet is LINEAR scale
    with np.errstate(invalid="ignore", divide="ignore"):
        expr_log2 = np.log2(expr.where(expr > 0))

    return {
        "annot": annot,
        "expr_log2": expr_log2,
        "meta": meta,
        "sample_cols": sample_cols,
        "id_col": id_col,
        "label_col": label_col,
        "level": level,
    }


def _sample_keys(name: str) -> list[str]:
    """Candidate identifier tokens for a sample name, highest-priority first:
    the full name, the trailing '-'/'_'-token, then every token (upper-cased)."""
    n = str(name).strip()
    parts = [p for p in re.split(r"[-_\s/]+", n) if p]
    raw = [n] + ([parts[-1]] if parts else []) + parts
    out, seen = [], set()
    for k in raw:
        ku = k.upper()
        if ku and ku not in seen:
            seen.add(ku)
            out.append(ku)
    return out


def _match_clinical_column(cl: pd.DataFrame, col, sample_tok: dict, name_up: dict) -> dict:
    """Map sample_id -> clinical row index for ONE clinical column, matching by value
    (exact token match, then substring fallback). Used both to auto-detect the key and
    to honour a manually chosen column."""
    ser = cl[col].dropna()
    if ser.empty:
        return {}
    val_to_row: dict[str, object] = {}
    for idx, v in ser.items():
        val_to_row.setdefault(str(v).strip().upper(), idx)
    long_vals = [(val, idx) for val, idx in val_to_row.items() if len(val) >= 3]
    mapping: dict = {}
    for sid, keys in sample_tok.items():
        hit = None
        for k in keys:                       # exact token match
            if k in val_to_row:
                hit = val_to_row[k]
                break
        if hit is None:                       # substring fallback (either direction)
            up = name_up[sid]
            for k in keys:
                if len(k) < 4:
                    continue
                for val, idx in long_vals:
                    if k in val or val in up:
                        hit = idx
                        break
                if hit is not None:
                    break
        if hit is not None:
            mapping[sid] = hit
    return mapping


def infer_clinical_key(meta: pd.DataFrame, cl: pd.DataFrame,
                       min_rate: float = 0.5) -> tuple[str | None, float, dict]:
    """Find the clinical column whose *values* best identify the PRISM samples.

    Column name is irrelevant — matching is by content. For each column it tries
    exact token matches (full sample name, trailing token, any token) then a
    substring fallback, and prefers columns that map ~one clinical row per sample
    (so a low-cardinality column like 'batch' can't win by coincidence).

    Returns (key_column, match_rate, {sample_id: clinical_row_label}).
    """
    sample_tok = {sid: _sample_keys(nm) for sid, nm in zip(meta["sample_id"], meta["sample"])}
    name_up = {sid: str(nm).upper() for sid, nm in zip(meta["sample_id"], meta["sample"])}
    n = len(meta)
    best_col, best_rate, best_map, best_score = None, 0.0, {}, -1.0

    for c in cl.columns:
        mapping = _match_clinical_column(cl, c, sample_tok, name_up)
        matched = len(mapping)
        if not matched:
            continue
        rate = matched / n
        one_to_one = len(set(mapping.values())) / matched  # penalise non-identifier columns
        score = rate * one_to_one
        if score > best_score:
            best_col, best_rate, best_map, best_score = c, rate, mapping, score

    if best_col is None or best_rate < min_rate:
        return None, best_rate, {}
    return best_col, best_rate, best_map


def attach_clinical(meta: pd.DataFrame, clinical_csv: str | Path,
                    key_col: str | None = None) -> tuple[pd.DataFrame, dict]:
    """Join a clinical metadata CSV onto sample metadata by matching the key column's
    values to the PRISM sample names.

    The identifier column in the CSV can be named anything (``replicate``,
    ``Sample Label``, ``UW Test Name``, ...). By default it is auto-detected by value.
    Pass ``key_col`` to force a specific column (overriding a wrong auto-detection).
    Returns (augmented meta indexed by sample_id, info) where
    info = {"match_rate", "key_column", "n_clinical_cols", "manual"}.
    """
    cl = pd.read_csv(clinical_csv)
    manual = bool(key_col) and key_col in cl.columns
    if manual:
        sample_tok = {sid: _sample_keys(nm)
                      for sid, nm in zip(meta["sample_id"], meta["sample"])}
        name_up = {sid: str(nm).upper() for sid, nm in zip(meta["sample_id"], meta["sample"])}
        mapping = _match_clinical_column(cl, key_col, sample_tok, name_up)
        rate = len(mapping) / len(meta) if len(meta) else 0.0
    else:
        key_col, rate, mapping = infer_clinical_key(meta, cl)
    merged = meta.copy()
    info = {"match_rate": float(rate), "key_column": key_col, "manual": manual,
            "n_clinical_cols": int(cl.shape[1] - (1 if key_col else 0))}
    if key_col is None:
        merged = merged.set_index("sample_id", drop=False)
        return merged, info

    rows = [mapping.get(sid) for sid in merged["sample_id"]]
    aligned = cl.reindex(rows)            # None -> NaN row; preserves per-column dtypes
    aligned.index = merged.index
    for c in cl.columns:
        if c == key_col:                  # drop the join key itself (redundant with 'sample')
            continue
        name = c if c not in merged.columns else f"{c}_clin"
        vals = pd.Series(aligned[c].values, index=merged.index)
        # coerce numeric-looking text columns (e.g. Age stored as strings) so they
        # enter models as continuous covariates rather than dozens of dummy levels
        if vals.dtype == object:
            coerced = pd.to_numeric(vals, errors="coerce")
            nonblank = vals.notna() & (vals.astype(str).str.strip() != "")
            # treat as numeric if most non-blank values parse as numbers (the rest,
            # e.g. QC/Pool labels in an Age column, become NaN — those samples aren't
            # experimental anyway). True categoricals (Sex, APOE) parse ~0% -> stay text.
            if nonblank.sum() and coerced[nonblank].notna().mean() >= 0.5:
                vals = coerced
        merged[name] = vals
    merged = merged.set_index("sample_id", drop=False)
    return merged, info


def cryptic_peptide_map(merged_parquet: str | Path, term: str = "cryptic",
                        peptide_col: str = "PeptideModifiedSequenceUnimodIds",
                        protein_col: str = "Protein") -> dict[str, str]:
    """Map each cryptic peptide to a representative cryptic protein string.

    Cryptic peptides are flagged at the transition/peptide level: their
    ``Protein`` field contains the term (e.g. ``CRYPTIC_UNIQUE|...``), which the
    protein rollup discards. Reads only the matching distinct rows from
    ``merged_data.parquet`` via DuckDB (server-side filter, cheap).
    """
    import duckdb

    q = (f'SELECT DISTINCT "{peptide_col}" AS pep, "{protein_col}" AS prot '
         f'FROM read_parquet(?) WHERE "{protein_col}" ILIKE ?')
    con = duckdb.connect()
    try:
        df = con.execute(q, [str(merged_parquet), f"%{term}%"]).df()
    finally:
        con.close()
    out: dict[str, str] = {}
    for pep, prot in zip(df["pep"].astype(str), df["prot"].astype(str)):
        out.setdefault(pep, prot)
    return out


def load_cryptic_gene_map(path: str | Path) -> dict:
    """Load a peptide -> {gene, accession, protein_name, cryptic_string} map.

    Built once from UniProt (see cryptic_gene_map.csv). Lets the cryptic peptides
    show their real gene (MPRIP) instead of an opaque TrEMBL accession, and links
    each to its canonical counterpart for stoichiometry.
    """
    p = Path(path)
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    out: dict = {}
    for _, r in df.iterrows():
        out[str(r["peptide"])] = {
            "gene": str(r.get("gene", "") or ""),
            "accession": str(r.get("cryptic_accession", "") or ""),
            "protein_name": str(r.get("canonical_protein_name", "") or ""),
            "cryptic_string": str(r.get("cryptic_string", "") or ""),
        }
    return out


def cryptic_short_label(protein_string: str) -> str:
    """Shorten a cryptic accession to a readable name, e.g. 'S35U4_HUMAN'."""
    parts = str(protein_string).split("|")
    for p in parts:
        if p.endswith("_HUMAN"):
            return p
    return parts[1] if len(parts) > 1 else str(protein_string)


def cryptic_detection_matrix(merged_parquet: str | Path, term: str | None = "cryptic",
                             q_thresh: float = 0.01,
                             peptide_col: str = "PeptideModifiedSequenceUnimodIds",
                             sample_col: str = "Sample ID") -> pd.DataFrame:
    """Binary detection matrix (peptide x sample_id).

    A cell is 1 if the peptide was genuinely *detected* in that sample
    (``DetectionQValue`` present and < q_thresh), not merely integrated from a
    borrowed RT window. This is the right signal for on/off peptides -- the dense
    Area matrix can't tell a real peak from a borrowed-baseline one. ``term=None``
    covers all peptides; a string restricts to cryptic (Protein ILIKE '%term%').
    """
    import duckdb

    thr = float(q_thresh)
    sel = (f'SELECT "{peptide_col}" AS pep, "{sample_col}" AS samp, '
           f'MAX(CASE WHEN "DetectionQValue" IS NOT NULL AND "DetectionQValue" < {thr} '
           f'THEN 1 ELSE 0 END) AS det FROM read_parquet(?) ')
    if term:
        q, params = sel + 'WHERE "Protein" ILIKE ? GROUP BY pep, samp', [str(merged_parquet), f"%{term}%"]
    else:
        q, params = sel + "GROUP BY pep, samp", [str(merged_parquet)]
    con = duckdb.connect()
    try:
        df = con.execute(q, params).df()
    finally:
        con.close()
    return df.pivot(index="pep", columns="samp", values="det").fillna(0).astype(int)


def detection_test(det_matrix: pd.DataFrame, group_a: list[str],
                   group_b: list[str]) -> pd.DataFrame:
    """Per-feature Fisher exact test of detection rate, group A vs B (BH-FDR)."""
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests

    a_cols = [s for s in group_a if s in det_matrix.columns]
    b_cols = [s for s in group_b if s in det_matrix.columns]
    rows = []
    for pep in det_matrix.index:
        na = int(det_matrix.loc[pep, a_cols].sum())
        nb = int(det_matrix.loc[pep, b_cols].sum())
        _, p = fisher_exact([[na, len(a_cols) - na], [nb, len(b_cols) - nb]])
        rows.append({"feature": pep, "det_A": na, "n_A": len(a_cols),
                     "det_B": nb, "n_B": len(b_cols),
                     "rate_A": na / len(a_cols) if a_cols else np.nan,
                     "rate_B": nb / len(b_cols) if b_cols else np.nan, "p": p})
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = multipletests(out["p"], method="fdr_bh")[1]
        out = out.sort_values("p").reset_index(drop=True)
    return out


def _firth_logit(X: np.ndarray, y: np.ndarray, max_iter: int = 500,
                 tol: float = 1e-7) -> tuple[np.ndarray, float, bool]:
    """Firth-penalized logistic regression (Jeffreys-prior bias reduction).

    Returns (beta, penalized log-likelihood, converged). Gives finite estimates
    under complete/quasi-separation (common at small n), unlike the MLE. Uses
    Newton steps with step-halving for stability.
    """
    n, p = X.shape
    beta = np.zeros(p)
    ll_prev = -np.inf
    converged = False

    def _penll(b):
        eta = X @ b
        pr = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-12, 1 - 1e-12)
        w = pr * (1 - pr)
        xtwx = X.T @ (X * w[:, None])
        try:
            inv = np.linalg.inv(xtwx)
        except np.linalg.LinAlgError:
            inv = np.linalg.pinv(xtwx)
        _, logdet = np.linalg.slogdet(xtwx)
        ll = float(np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr)) + 0.5 * logdet)
        return pr, w, inv, ll

    for it in range(max_iter):
        pr, w, inv, ll = _penll(beta)
        if it > 0 and abs(ll - ll_prev) < tol:   # converge on the penalized log-likelihood
            converged = True
            break
        ll_prev = ll
        h = w * np.einsum("ij,jk,ik->i", X, inv, X)        # hat-matrix diagonal
        step = inv @ (X.T @ (y - pr + h * (0.5 - pr)))     # Firth-modified Newton step
        if not np.all(np.isfinite(step)):
            break
        if np.max(np.abs(step)) > 5:                       # damping against overshoot
            step *= 5.0 / np.max(np.abs(step))
        beta = beta + step
    _, _, _, ll = _penll(beta)
    return beta, ll, converged


def _build_covariate_design(covariate_df: pd.DataFrame | None, samples: list[str],
                            grp: np.ndarray) -> tuple[np.ndarray | None, list[str], list[str]]:
    """Numeric (centered) + categorical (drop-first dummy) covariate matrix for `samples`.

    Drops covariates with missing values in these samples, constant columns, and
    dummy levels collinear with the group vector. Returns (X, names, dropped_msgs).
    """
    if covariate_df is None or covariate_df.shape[1] == 0:
        return None, [], []
    cols, names, dropped = [], [], []
    for c in covariate_df.columns:
        s = covariate_df.reindex(samples)[c]
        if pd.api.types.is_numeric_dtype(s):
            if s.isna().any():
                dropped.append(f"{c} (missing numeric values)")
                continue
            v = s.to_numpy(dtype=float)
            if np.unique(v).size < 2:
                continue
            cols.append((v - v.mean())[:, None])
            names.append(c)
        else:
            sv = s.astype("object").where(s.notna(), "NA").astype(str)  # NaN -> its own level
            d = pd.get_dummies(sv, prefix=c, drop_first=True)
            for col in d.columns:
                cv = d[col].to_numpy(dtype=float)
                if np.unique(cv).size < 2:
                    continue
                # collapse rare levels (singletons) into the reference — they cause
                # quasi-separation, over-parameterise small designs, and aren't identifiable
                if cv.sum() < 3 or cv.sum() > len(cv) - 3:
                    dropped.append(f"{col} (rare level, n<3 — pooled into reference)")
                    continue
                if np.allclose(cv, grp) or np.allclose(cv, 1 - grp):
                    dropped.append(f"{col} (confounded with group)")
                    continue
                cols.append(cv[:, None])
                names.append(col)
    if not cols:
        return None, [], dropped
    return np.hstack(cols), names, dropped


def detection_test_glm(det_matrix: pd.DataFrame, group_a: list[str], group_b: list[str],
                       covariate_df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    """Per-feature Firth logistic regression of detection on group (+ covariates).

    Tests ``detected ~ group [+ cohort/age/...]`` with a penalized likelihood-ratio
    test for the group term (group B = 1, so logOR > 0 = more detected in B).
    Adjusts for confounders (e.g. Study_Name cohort) the Fisher test cannot.
    Returns (results sorted by p with BH-q, info dict).
    """
    from scipy.stats import chi2

    a_cols = [s for s in group_a if s in det_matrix.columns]
    b_cols = [s for s in group_b if s in det_matrix.columns]
    samples = a_cols + b_cols
    grp = np.array([0.0] * len(a_cols) + [1.0] * len(b_cols))
    xcov, cov_names, dropped = _build_covariate_design(covariate_df, samples, grp)
    info = {"covariates_used": cov_names, "dropped": dropped, "n_A": len(a_cols), "n_B": len(b_cols)}

    base = [np.ones(len(samples))] + ([xcov] if xcov is not None else [])
    x_red = np.column_stack(base)
    x_full = np.column_stack([np.ones(len(samples)), grp] + ([xcov] if xcov is not None else []))

    # identifiability: is the group term separable from the covariates, and is the
    # design full-rank for this n? (catches cohort-confounded contrasts like FTD-vs-CO)
    n_s = len(samples)
    r_full = int(np.linalg.matrix_rank(x_full))
    r_red = int(np.linalg.matrix_rank(x_red))
    info["n_params"] = int(x_full.shape[1])
    info["identifiable"] = (r_full > r_red) and (r_full == x_full.shape[1]) and (n_s > x_full.shape[1])
    if x_red.shape[1] > 1 and grp.var() > 0:
        coef, *_ = np.linalg.lstsq(x_red, grp, rcond=None)
        info["group_collinearity_r2"] = float(max(0.0, 1.0 - ((grp - x_red @ coef).var() / grp.var())))
    else:
        info["group_collinearity_r2"] = 0.0
    if not info["identifiable"]:
        info["warning"] = ("Group is not separable from the covariates (or n is too small for "
                           "the design) — the adjusted detection test is not identifiable. "
                           "The contrast is confounded; use the unadjusted view with caution or "
                           "balance the groups.")
        empty = pd.DataFrame(columns=["feature", "det_A", "n_A", "det_B", "n_B",
                                      "rate_A", "rate_B", "logOR", "p", "q"])
        return empty, info

    rows = []
    for pep in det_matrix.index:
        y = det_matrix.loc[pep, samples].to_numpy(dtype=float)
        na = int(y[: len(a_cols)].sum())
        nb = int(y[len(a_cols):].sum())
        if y.sum() == 0 or y.sum() == len(y):
            logor, p = 0.0, 1.0  # no detection variation -> nothing to test
        else:
            try:
                bf, llf, cvf = _firth_logit(x_full, y)
                _, llr, cvr = _firth_logit(x_red, y)
                p = float(chi2.sf(max(2.0 * (llf - llr), 0.0), 1))
                logor = float(bf[1])
                if not (cvf and cvr):
                    info["n_nonconverged"] = info.get("n_nonconverged", 0) + 1
            except Exception:  # noqa: BLE001
                logor, p = np.nan, np.nan
        rows.append({"feature": pep, "det_A": na, "n_A": len(a_cols), "det_B": nb,
                     "n_B": len(b_cols), "rate_A": na / len(a_cols) if a_cols else np.nan,
                     "rate_B": nb / len(b_cols) if b_cols else np.nan, "logOR": logor, "p": p})
    out = pd.DataFrame(rows)
    ok = out["p"].notna()
    out.loc[ok, "q"] = multipletests(out.loc[ok, "p"], method="fdr_bh")[1]
    out = out.sort_values("p").reset_index(drop=True)
    return out, info


def compute_pca(
    expr_log2: pd.DataFrame, sample_ids: list[str], n_components: int = 6
) -> tuple[pd.DataFrame, np.ndarray, int]:
    """PCA of samples on complete-case features (feature-centred log2).

    Returns (scores DataFrame indexed by sample_id with PC1..PCk columns,
    variance-explained ratios, n_features_used).
    """
    M = expr_log2[sample_ids].dropna(axis=0, how="any")
    if M.shape[0] < 2 or M.shape[1] < 2:
        raise ValueError("Not enough complete features / samples for PCA.")
    X = M.to_numpy(dtype=float).T  # samples x features
    X = X - X.mean(axis=0, keepdims=True)  # centre each feature
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    k = int(min(n_components, S.shape[0]))
    scores = U[:, :k] * S[:k]
    var_ratio = (S**2 / np.sum(S**2))[:k]
    return (
        pd.DataFrame(scores, index=sample_ids, columns=[f"PC{i + 1}" for i in range(k)]),
        var_ratio,
        int(M.shape[0]),
    )


def match_genes(
    annot: pd.DataFrame, label_col: str, expr_log2: pd.DataFrame, genes: list[str]
) -> tuple[dict[str, str], list[str]]:
    """Map gene symbols to one representative feature id each (highest mean log2).

    Handles multi-gene leading labels like ``"ANXA2 / ANXA2"`` by splitting on '/'.
    Returns (dict gene->feature_id, list of missing genes). Matching is case-insensitive.
    """
    if label_col not in annot.columns:
        return {}, list(genes)
    gmap: dict[str, list[str]] = {}
    for fid, gname in annot[label_col].astype(str).items():
        for g in gname.split("/"):
            g = g.strip().upper()
            if g:
                gmap.setdefault(g, []).append(fid)
    found, missing = {}, []
    for g in genes:
        fids = list(dict.fromkeys(gmap.get(g.strip().upper(), [])))
        if not fids:
            missing.append(g)
        elif len(fids) == 1:
            found[g] = fids[0]
        else:
            found[g] = expr_log2.loc[fids].mean(axis=1).idxmax()
    return found, missing


def covariate_balance(meta: pd.DataFrame, group_a, group_b, col: str) -> pd.DataFrame | None:
    """Crosstab of a categorical metadata column against group A/B selection."""
    if col not in meta.columns:
        return None
    sel = pd.Series({**{s: "A" for s in group_a}, **{s: "B" for s in group_b}}, name="group")
    sub = meta.loc[sel.index].copy()
    sub["group"] = sel.values
    return pd.crosstab(sub[col].astype(str), sub["group"])


def differential(
    expr_log2: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    covariates: pd.DataFrame | None = None,
    min_per_group: int = 2,
) -> tuple[pd.DataFrame, dict]:
    """limma empirical-Bayes moderated t-test: group B vs group A.

    Parameters
    ----------
    expr_log2 : features x samples (log2); columns are sample_ids.
    group_a, group_b : sample_id lists. B is "treatment" (+logFC = higher in B).
    covariates : optional DataFrame (index=sample_id). Numeric columns enter the
        model centred; categorical columns enter as drop-first dummies. Columns
        that are constant, confounded with group, or contain missing values in
        the selected samples are dropped with a logged message.
    min_per_group : minimum non-NaN samples required per group for a feature.

    Returns (results DataFrame sorted by P.Value, info dict).
    """
    info: dict = {"messages": []}
    cols = list(group_a) + list(group_b)
    if len(set(cols)) != len(cols):
        raise ValueError("A sample is assigned to both groups. Groups must be disjoint.")
    if len(group_a) < min_per_group or len(group_b) < min_per_group:
        raise ValueError(
            f"Each group needs at least {min_per_group} samples "
            f"(got A={len(group_a)}, B={len(group_b)})."
        )

    M = expr_log2[cols].copy()
    grp = np.array([0.0] * len(group_a) + [1.0] * len(group_b))
    design = {"Intercept": np.ones(len(cols)), "groupB": grp}

    if covariates is not None and covariates.shape[1] > 0:
        cov = covariates.reindex(cols)
        for name in cov.columns:
            s = cov[name]
            if s.isna().any():
                info["messages"].append(
                    f"Covariate '{name}' has missing values in selected samples - skipped."
                )
                continue
            if pd.api.types.is_numeric_dtype(s):
                v = s.to_numpy(dtype=float)
                if np.unique(v).size < 2:
                    info["messages"].append(f"Covariate '{name}' is constant - skipped.")
                    continue
                design[name] = v - v.mean()
            else:
                dummies = pd.get_dummies(s.astype(str), prefix=name, drop_first=True)
                for c in dummies.columns:
                    col = dummies[c].to_numpy(dtype=float)
                    if np.unique(col).size < 2:
                        continue
                    if np.allclose(col, grp) or np.allclose(col, 1 - grp):
                        info["messages"].append(
                            f"Covariate level '{c}' is confounded with group - dropped."
                        )
                        continue
                    design[c] = col

    design_df = pd.DataFrame(design, index=cols)
    coef_idx = list(design_df.columns).index("groupB")
    n, p = design_df.shape
    if n - p < 1:
        raise ValueError(
            f"Not enough residual degrees of freedom (n={n}, params={p}). "
            "Use more samples or fewer covariates."
        )
    if np.linalg.matrix_rank(design_df.to_numpy(dtype=float)) < p:
        raise ValueError(
            "Design matrix is rank-deficient (covariates collinear with each other "
            "or with group). Remove a covariate."
        )
    info["covariates_used"] = [c for c in design_df.columns if c not in ("Intercept", "groupB")]

    a_ok = M[group_a].notna().sum(axis=1)
    b_ok = M[group_b].notna().sum(axis=1)
    keep = M.notna().all(axis=1) & (a_ok >= min_per_group) & (b_ok >= min_per_group)
    info["n_features_total"] = int(len(M))
    info["n_features_tested"] = int(keep.sum())
    info["n_features_dropped"] = int((~keep).sum())
    if keep.sum() == 0:
        raise ValueError("No features had complete data across the selected samples.")

    Mk = M.loc[keep]
    fit = L.lmFit(Mk.to_numpy(dtype=float), design_df.to_numpy(dtype=float))
    sigma = np.asarray(fit.sigma, dtype=float)
    df_res = np.asarray(fit.df_residual, dtype=float)
    sv = L.squeezeVar(sigma**2, df_res)
    s2_post = np.asarray(sv["var_post"], dtype=float)
    df_prior = np.ravel(np.asarray(sv["df_prior"], dtype=float))

    coef = np.asarray(fit.coefficients, dtype=float)[:, coef_idx]
    stdev_unscaled = np.asarray(fit.stdev_unscaled, dtype=float)[:, coef_idx]
    t = coef / (stdev_unscaled * np.sqrt(s2_post))
    df_total = df_res + df_prior
    pval = 2.0 * st.t.cdf(-np.abs(t), df=df_total)
    qval = multipletests(pval, method="fdr_bh")[1]

    res = pd.DataFrame(
        {
            "logFC": coef,
            "FC": np.power(2.0, coef),
            "AveExpr": np.asarray(fit.Amean, dtype=float),
            "t": t,
            "P.Value": pval,
            "adj.P.Val": qval,
            "mean_A": Mk[group_a].mean(axis=1).to_numpy(),
            "mean_B": Mk[group_b].mean(axis=1).to_numpy(),
            "n_A": int(len(group_a)),
            "n_B": int(len(group_b)),
        },
        index=Mk.index,
    )
    info["df_prior"] = float(np.median(df_prior))
    info["df_residual"] = float(np.median(df_res))
    return res.sort_values("P.Value"), info


# --- Significance finder ---------------------------------------------------
# Scans many ways to split one categorical column's levels into two groups,
# runs limma for each, and ranks by signal. EXPLORATORY: searching contrasts
# inflates false positives, so the UI calibrates the top hit with a label
# permutation (how much signal a random split of the same sizes gives).

def _enumerate_splits(levels: list[str], scope: str, max_tests: int):
    """Yield (a_levels, b_levels) partitions of `levels` for the given scope."""
    import itertools
    K = len(levels)
    if scope == "One-vs-rest":
        for lv in levels:
            yield [lv], [x for x in levels if x != lv]
    elif scope == "Pairwise":
        for a, b in itertools.combinations(levels, 2):
            yield [a], [b]
    else:  # "Subset pairs (thorough)": disjoint non-empty A,B (levels may be excluded)
        count = 0
        seen = set()
        idx = list(range(K))
        # assign each level to A(0)/B(1)/exclude(2); dedupe A<->B mirrors
        for assign in itertools.product((0, 1, 2), repeat=K):
            a = [levels[i] for i in idx if assign[i] == 0]
            b = [levels[i] for i in idx if assign[i] == 1]
            if not a or not b:
                continue
            key = frozenset((frozenset(a), frozenset(b)))
            if key in seen:
                continue
            seen.add(key)
            yield a, b
            count += 1
            if count >= max_tests:
                return


def significance_scan(expr_log2: pd.DataFrame, meta: pd.DataFrame, col: str, scope: str,
                      covariates: pd.DataFrame | None = None, q_cut: float = 0.05,
                      lfc_cut: float = 1.0, min_n: int = 3, max_tests: int = 400,
                      experimental_only: bool = True) -> tuple[pd.DataFrame, bool]:
    """Run limma for many splits of `col`; return (ranked results, truncated?).

    Results columns: A, B, n_A, n_B, n_sig, min_q, plus hidden _a/_b level lists.
    """
    pool = (meta[meta["sample_type"] == "experimental"]
            if experimental_only and "sample_type" in meta.columns else meta)
    vc = pool[col].dropna().astype(str).value_counts()
    levels = [lv for lv in vc.index if vc[lv] >= min_n]
    colvals = pool[col].astype(str)
    ids_by_level = {lv: pool.index[colvals == lv].tolist() for lv in levels}

    rows, n_done, truncated = [], 0, False
    for a_lvls, b_lvls in _enumerate_splits(levels, scope, max_tests):
        ga = [s for lv in a_lvls for s in ids_by_level[lv]]
        gb = [s for lv in b_lvls for s in ids_by_level[lv]]
        if len(ga) < min_n or len(gb) < min_n:
            continue
        n_done += 1
        if n_done > max_tests:
            truncated = True
            break
        try:
            res, _ = differential(expr_log2, ga, gb, covariates=covariates, min_per_group=min_n)
        except Exception:  # noqa: BLE001
            continue
        nsig = int(((res["adj.P.Val"] < q_cut) & (res["logFC"].abs() >= lfc_cut)).sum())
        rows.append({"A": " + ".join(a_lvls), "B": " + ".join(b_lvls),
                     "n_A": len(ga), "n_B": len(gb), "n_sig": nsig,
                     "min_q": float(res["adj.P.Val"].min()),
                     "_a": a_lvls, "_b": b_lvls})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["n_sig", "min_q"], ascending=[False, True]).reset_index(drop=True)
    return df, truncated


def scan_summary(scan_df: pd.DataFrame | None) -> dict:
    """Headline counts for a `significance_scan` result (sorted by n_sig desc).

    Returns n_contrasts, best_n_sig (significant features in the best contrast),
    n_with_hits (contrasts with >=1 significant feature), and the best contrast's
    A/B labels. Lets the UI state plainly how many features are significant rather
    than only listing per-contrast n_sig.
    """
    if scan_df is None or len(scan_df) == 0:
        return {"n_contrasts": 0, "best_n_sig": 0, "n_with_hits": 0,
                "best_A": None, "best_B": None}
    top = scan_df.iloc[0]
    return {
        "n_contrasts": int(len(scan_df)),
        "best_n_sig": int(scan_df["n_sig"].max()),
        "n_with_hits": int((scan_df["n_sig"] > 0).sum()),
        "best_A": str(top["A"]),
        "best_B": str(top["B"]),
    }


def permute_calibrate(expr_log2: pd.DataFrame, group_a: list[str], group_b: list[str],
                      covariates: pd.DataFrame | None, q_cut: float, lfc_cut: float,
                      n_perm: int = 200, min_n: int = 3) -> list[int]:
    """Null distribution of n_sig under random label shuffles of the same group sizes.

    NB: this calibrates ONE split's labels — it does not account for the search over
    many partitions. For family-wise calibration use scan_permutation_null.
    """
    pool = list(group_a) + list(group_b)
    n_a = len(group_a)
    rng = np.random.default_rng(12345)
    null = []
    for _ in range(n_perm):
        perm = list(rng.permutation(pool))
        pa, pb = perm[:n_a], perm[n_a:]
        try:
            res, _ = differential(expr_log2, pa, pb, covariates=covariates, min_per_group=min_n)
            null.append(int(((res["adj.P.Val"] < q_cut) & (res["logFC"].abs() >= lfc_cut)).sum()))
        except Exception:  # noqa: BLE001
            pass
    return null


def scan_permutation_null(expr_log2: pd.DataFrame, meta: pd.DataFrame, col: str, scope: str,
                          covariates: pd.DataFrame | None = None, q_cut: float = 0.05,
                          lfc_cut: float = 1.0, min_n: int = 3, max_tests: int = 400,
                          n_perm: int = 100, seed: int = 12345) -> list[int]:
    """Family-wise null for a Finder search: permute the column labels, RE-RUN the whole
    scan, and record the **max** n_sig per permutation. The observed top hit's empirical
    p against this null accounts for the selection over all partitions (the within-split
    permutation does not).

    Labels are permuted across the experimental samples while each sample keeps its own
    covariate values (so the null breaks the group↔feature association; with strong
    covariate↔group confounding this is approximate — interpret alongside identifiability).
    """
    pool_mask = (meta["sample_type"] == "experimental") if "sample_type" in meta.columns \
        else pd.Series(True, index=meta.index)
    ids = meta.index[pool_mask].tolist()
    labels = meta.loc[ids, col].astype("object").to_numpy()
    rng = np.random.default_rng(seed)
    meta_perm = meta.copy()
    null_max = []
    for _ in range(n_perm):
        meta_perm.loc[ids, col] = rng.permutation(labels)
        df, _ = significance_scan(expr_log2, meta_perm, col, scope, covariates=covariates,
                                  q_cut=q_cut, lfc_cut=lfc_cut, min_n=min_n, max_tests=max_tests)
        null_max.append(int(df["n_sig"].max()) if len(df) else 0)
    return null_max


# --------------------------------------------------------------------------- #
# Streamlit UI
# --------------------------------------------------------------------------- #

# Sensible defaults (edit freely in the UI). Currently the LatimerCryptic plasma study —
# cryptic peptides need Feature level = peptide.
_DEFAULT_OUTDIR = r"G:\Manuscripts\LatimerCryptic\plasma_analysis\prism_output_dir"
_DEFAULT_CLINICAL = r"G:\Manuscripts\LatimerCryptic\Metadata\Req242_Plasma_MS_merged_metadata.csv"


# --- Native file/folder picker (local runs only) -------------------------------
# The OS dialog runs in a SEPARATE process so tkinter never touches Streamlit's
# script thread (mixing them crashes/hangs). The dialog opens on the machine
# running the server, so this only works when the app runs locally — on a remote
# deploy just type the path into the box.
_PICK_CODE = r"""
import sys, os, tkinter as tk
from tkinter import filedialog
kind, initial = sys.argv[1], sys.argv[2]
root = tk.Tk(); root.withdraw(); root.wm_attributes("-topmost", 1)
if kind == "dir":
    picked = filedialog.askdirectory(initialdir=initial or None,
                                     title="Select PRISM output_dir")
else:
    start = os.path.dirname(initial) if initial else None
    picked = filedialog.askopenfilename(initialdir=start, title="Select a CSV file",
                                        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
root.destroy()
sys.stdout.write(picked or "")
"""


def longitudinal_time_slope(expr_log2: pd.DataFrame, subject: pd.Series, time: pd.Series,
                            min_timepoints: int = 2) -> tuple[pd.DataFrame, dict]:
    """Within-subject time trend per feature: fit ``~ subject + time`` (fixed subject
    intercepts + continuous time) and test the time-slope coefficient.

    `subject` and `time` are per-sample Series indexed by expr_log2 columns. Subjects
    with fewer than `min_timepoints` distinct timepoints carry no within-subject slope
    and are dropped. EXPLORATORY: returns slope (log2 per time unit), raw P and BH adj.P
    — small cohorts rarely clear FDR. Returns (results sorted by P.Value, info dict).
    """
    from statsmodels.stats.multitest import multipletests

    subj = subject.reindex(expr_log2.columns).astype(str)
    tnum = pd.to_numeric(time.reindex(expr_log2.columns), errors="coerce")
    ok = (subj.notna() & ~subj.str.lower().isin(["nan", "na", "none", "", "<na>"])
          & tnum.notna())
    keep_subj = [s for s, g in tnum[ok].groupby(subj[ok]) if g.nunique() >= min_timepoints]
    ok = ok & subj.isin(keep_subj)
    cols = [c for c in expr_log2.columns if bool(ok.get(c, False))]
    if len(keep_subj) < 2 or len(cols) < 4:
        raise ValueError(
            f"Need at least 2 subjects each with >={min_timepoints} distinct timepoints "
            f"(got {len(keep_subj)} usable subject(s), {len(cols)} sample(s)). Check the "
            "subject and timepoint columns.")
    subj = subj[cols]
    tnum = tnum[cols]

    M = expr_log2[cols].dropna(how="any")
    if M.shape[0] < 2:
        raise ValueError("No features are complete across the longitudinal samples.")
    Y = M.to_numpy(dtype=float)
    dummies = pd.get_dummies(pd.Categorical(subj.values), drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(cols)), dummies, tnum.to_numpy(float)])
    time_i = X.shape[1] - 1
    dof = X.shape[0] - X.shape[1]
    if dof < 1:
        raise ValueError(f"Not enough residual df (samples={len(cols)}, params={X.shape[1]}). "
                         "Need more samples relative to subjects.")
    xtx_inv = np.linalg.inv(X.T @ X)
    betas = Y @ (xtx_inv @ X.T).T
    resid = Y - betas @ X.T
    s2 = np.einsum("ij,ij->i", resid, resid) / dof
    se = np.sqrt(s2 * xtx_inv[time_i, time_i])
    slope = betas[:, time_i]
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, slope / se, np.nan)
    pval = 2.0 * st.t.sf(np.abs(tstat), dof)   # module-level st is scipy.stats
    qval = np.full(pval.shape, np.nan)
    good = np.isfinite(pval)
    if good.any():
        qval[good] = multipletests(pval[good], method="fdr_bh")[1]
    res = pd.DataFrame({"slope": slope, "t": tstat, "P.Value": pval, "adj.P.Val": qval,
                        "AveExpr": Y.mean(axis=1)}, index=M.index)
    info = {"n_subjects": len(keep_subj), "n_samples": len(cols),
            "residual_df": int(dof), "n_features": int(M.shape[0]),
            "subjects": keep_subj}
    return res.sort_values("P.Value"), info


def heterogeneity_table(expr_log2: pd.DataFrame, det: pd.DataFrame,
                        individuals: list[str]) -> pd.DataFrame:
    """Per-feature heterogeneity across individual samples, from genuine detection.

    `det` is a binary detection matrix (feature x sample; 1 = real detection, i.e.
    DetectionQValue below threshold — NOT the dense borrowed-baseline matrix). Returns
    per feature: prevalence (fraction of individuals that genuinely detect it), n_detected,
    SD and mean of log2 among detecting individuals, and pooling_dilution (~1/prevalence).
    Low-prevalence features are exactly the ones a pooled sample averages toward baseline.
    """
    idv = [s for s in individuals if s in det.columns and s in expr_log2.columns]
    n = len(idv)
    D = det.reindex(index=expr_log2.index, columns=idv).fillna(0).astype(bool)
    n_det = D.sum(axis=1)
    prev = (n_det / n) if n else n_det.astype(float)
    detected_log2 = expr_log2[idv].where(D)
    with np.errstate(invalid="ignore"):
        sd = detected_log2.std(axis=1)
        mean_det = detected_log2.mean(axis=1)
    return pd.DataFrame({
        "prevalence": prev,
        "n_detected": n_det.astype(int),
        "n_individuals": n,
        "sd_when_detected": sd,
        "mean_log2_detected": mean_det,
        "pooling_dilution": 1.0 / prev.where(prev > 0),
    })


def _build_html_report(title: str, intro: str, sections: list) -> str:
    """Assemble a self-contained HTML report (Plotly.js embedded inline once, so the
    file works offline). `sections` is a list of dicts: {'heading', 'paras' (list[str]),
    'figs' (list[plotly fig]), 'table' (DataFrame or None)}.
    """
    import html as _html

    from plotly.offline import get_plotlyjs

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{_html.escape(title)}</title>",
        "<style>body{font-family:system-ui,Segoe UI,Arial,sans-serif;max-width:940px;margin:2rem "
        "auto;padding:0 1.2rem;color:#1a1a1a;line-height:1.5}h1{font-size:1.7rem}h2{font-size:"
        "1.25rem;border-bottom:1px solid #eee;padding-bottom:.3rem;margin-top:2rem}table{border-"
        "collapse:collapse;font-size:.82rem;margin:.5rem 0}td,th{border:1px solid #ddd;padding:"
        "3px 8px;text-align:right}th{background:#f6f6f6}.note{color:#666;font-size:.82rem}</style>"
        f"<script type='text/javascript'>{get_plotlyjs()}</script></head><body>",
        f"<h1>{_html.escape(title)}</h1>",
    ]
    if intro:
        parts.append(f"<p>{_html.escape(intro)}</p>")
    fig_i = 0
    for sec in sections:
        parts.append(f"<h2>{_html.escape(sec['heading'])}</h2>")
        for para in sec.get("paras", []):
            parts.append(f"<p>{_html.escape(para)}</p>")
        for fig in sec.get("figs", []):
            parts.append(fig.to_html(full_html=False, include_plotlyjs=False,
                                     div_id=f"pdxfig{fig_i}"))
            fig_i += 1
        tbl = sec.get("table")
        if tbl is not None:
            parts.append(tbl.to_html(border=0))
        raw = sec.get("html")
        if raw:
            parts.append(raw)  # trusted, pre-built fragment (e.g. the feature explorer)
    parts.append("<p class='note'>Generated by PRISM Differential Explorer. Charts are interactive "
                 "(hover / zoom). Detection-based metrics use DetectionQValue, not the dense "
                 "matrix.</p></body></html>")
    return "\n".join(parts)


def _feature_explorer_html(expr_log2, group_a, group_b, label_a, label_b, annot,
                           preselect=None, block_id="pdxexp", max_cells=4_000_000):
    """Build a self-contained interactive box-plot explorer as an HTML fragment.

    Embeds the contrast's log2 matrix (group_a + group_b samples only) plus a small
    JS renderer so a reader can search ANY feature by accession/sequence and see its
    box plot client-side — no server, no recompute. Returns (html_fragment, note).

    If the matrix would exceed `max_cells`, keeps the most-variable features so the
    file stays shareable; the cap is disclosed in the returned note (never silent).
    """
    import html as _html
    import json

    import numpy as np
    import pandas as pd

    samples = [s for s in list(group_a) + list(group_b) if s in expr_log2.columns]
    groups = [label_a if s in set(group_a) else label_b for s in samples]
    sub = expr_log2[samples]

    note = ""
    n_samp = max(len(samples), 1)
    cap = max(int(max_cells // n_samp), 200)
    if len(sub) > cap:
        keep = sub.var(axis=1, skipna=True).sort_values(ascending=False).head(cap).index
        # keep the preselected feature even if it isn't top-variance
        if preselect is not None and preselect in sub.index and preselect not in set(keep):
            keep = keep.insert(0, preselect)
        sub = sub.loc[keep]
        note = (f"Feature explorer capped to the {len(sub):,} most-variable features "
                f"(of {len(expr_log2):,}) to keep the file shareable.")

    feats = list(sub.index)
    # searchable label per feature: id + any string annotation columns (truncated)
    lab_cols = [c for c in annot.columns
                if annot[c].dtype == object or str(annot[c].dtype).startswith("string")]
    ann = annot.reindex(feats)
    labels = []
    for f in feats:
        bits, seen = [str(f)], {str(f).lower()}
        for c in lab_cols:
            v = ann.at[f, c] if f in ann.index else None
            s = str(v)[:40] if v is not None else ""
            if s.lower() not in ("nan", "none", "") and s.lower() not in seen:
                bits.append(s)
                seen.add(s.lower())
        labels.append(" | ".join(bits)[:140])

    arr = sub.to_numpy(dtype=float)
    values = [[None if (x != x) else round(float(x), 2) for x in row] for row in arr]

    if preselect is not None and preselect in feats:
        pre_i = feats.index(preselect)
    else:
        pre_i = 0

    data = {
        "samples": samples, "groups": groups,
        "glabels": [label_a, label_b],
        "colors": {label_a: "#1f77b4", label_b: "#d62728"},
        "labels": labels, "values": values, "pre": pre_i,
    }
    payload = json.dumps(data, separators=(",", ":"))
    # Neutralise HTML-significant chars so a feature id/annotation containing "</script>"
    # (or "<!--") can't break out of the inline <script>. Escaping <,>,& to \uXXXX keeps
    # the JSON valid (single backslash) and JS parses them back to the literal chars.
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    bid = block_id
    frag = (
        f"<div class='pdx-explorer' id='{bid}'>"
        f"<p class='note'>Type an accession or peptide sequence, pick a feature, and its "
        f"box plot ({_html.escape(label_b)} vs {_html.escape(label_a)}) renders below. "
        f"Every feature in the contrast is here — search any of them.</p>"
        f"<input id='{bid}-q' type='text' placeholder='search accession / sequence…' "
        f"style='width:100%;max-width:520px;padding:6px 8px;font-size:.9rem;"
        f"border:1px solid #ccc;border-radius:4px'/>"
        f"<select id='{bid}-sel' size='8' style='display:block;width:100%;max-width:520px;"
        f"margin:.4rem 0;font-size:.82rem'></select>"
        f"<div id='{bid}-plot' style='height:440px'></div></div>"
        f"<script>(function(){{"
        f"var D={payload};var BID={json.dumps(bid)};"
        f"var q=document.getElementById(BID+'-q');"
        f"var sel=document.getElementById(BID+'-sel');"
        f"var plot=document.getElementById(BID+'-plot');"
        f"var low=D.labels.map(function(s){{return s.toLowerCase();}});"
        f"function fill(idxs){{sel.innerHTML='';idxs.slice(0,300).forEach(function(i){{"
        f"var o=document.createElement('option');o.value=i;o.textContent=D.labels[i];"
        f"sel.appendChild(o);}});}}"
        f"function render(i){{var vals=D.values[i];var traces=D.glabels.map(function(g){{"
        f"var ys=[];var txt=[];for(var k=0;k<D.samples.length;k++){{"
        f"if(D.groups[k]===g&&vals[k]!==null){{ys.push(vals[k]);txt.push(D.samples[k]);}}}}"
        f"return{{type:'box',name:g,y:ys,text:txt,boxpoints:'all',jitter:0.5,pointpos:0,"
        f"marker:{{color:D.colors[g],size:6}},line:{{color:D.colors[g]}},"
        f"hovertemplate:'%{{text}}<br>log2 %{{y:.2f}}<extra>'+g+'</extra>'}};}});"
        f"Plotly.newPlot(plot,traces,{{title:D.labels[i],height:440,"
        f"yaxis:{{title:'log2 abundance'}},margin:{{l:50,r:10,t:40,b:30}},"
        f"showlegend:false}},{{responsive:true,displaylogo:false}});}}"
        f"q.addEventListener('input',function(){{var t=q.value.toLowerCase().trim();"
        f"var idxs=[];for(var i=0;i<low.length;i++){{if(!t||low[i].indexOf(t)>=0)idxs.push(i);"
        f"if(idxs.length>=300)break;}}fill(idxs);}});"
        f"sel.addEventListener('change',function(){{if(sel.value!=='')render(+sel.value);}});"
        f"fill(D.values.map(function(_,i){{return i;}}));"
        f"sel.value=D.pre;render(D.pre);"
        f"}})();</script>"
    )
    return frag, note


def _viewer_config():
    """Return the viewer-mode config dict, or None for the normal (local) app.

    Viewer mode requires the env var ``PDX_VIEWER=1`` (so a normal local launch is never
    affected). The config is read from, in order: the JSON file named by
    ``PDX_VIEWER_CONFIG`` (or ``.share_config.json`` in the app dir) that the in-app
    "Share this analysis" button writes; otherwise a ``[viewer]`` section in
    ``.streamlit/secrets.toml``. It fixes the data source, hides the sidebar
    load/Browse controls (a visitor can't point the app elsewhere), optionally pins a
    comparison, and optionally password-gates. Keys: password, out_dir, clin_path,
    level, display_name, comparison{grouping_column, group_a_values, group_b_values,
    suggested_covariates, label_a, label_b}, q_cut, lfc_cut.
    """
    import json
    import os
    from pathlib import Path

    import streamlit as st
    if not os.environ.get("PDX_VIEWER"):
        return None
    cfg_path = os.environ.get("PDX_VIEWER_CONFIG") or str(
        Path(__file__).with_name(".share_config.json"))
    try:
        if Path(cfg_path).is_file():
            with open(cfg_path, encoding="utf-8") as fh:
                v = json.load(fh)
            return v if v.get("out_dir") else None
    except Exception:  # noqa: BLE001 - unreadable/partial file -> fall through
        pass
    try:
        v = st.secrets.get("viewer")
    except Exception:  # noqa: BLE001 - no secrets.toml at all
        return None
    if not v or not v.get("out_dir"):
        return None
    return dict(v)


def _viewer_gate(vcfg) -> bool:
    """Render a password prompt in viewer mode; True once authorised (or no password set)."""
    import hmac

    import streamlit as st
    pw = vcfg.get("password")
    if not pw:
        return True  # view-only but open — data source is still locked
    if st.session_state.get("_viewer_ok"):
        return True
    st.markdown("#### 🔒 Shared, read-only PRISM analysis")
    st.caption("Enter the password your colleague gave you to view this analysis.")
    entered = st.text_input("Password", type="password", key="_viewer_pw",
                            label_visibility="collapsed")
    if entered:
        if hmac.compare_digest(str(entered), str(pw)):
            st.session_state["_viewer_ok"] = True
            return True
        st.error("Incorrect password.")
    return False


def _pick_path(kind: str, initial: str) -> str:
    """Open a native folder ('dir') or file ('file') picker. '' if cancelled/failed."""
    import os
    import subprocess
    import sys
    if os.environ.get("PDX_NO_PICKER"):
        return ""  # test seam: skip the native dialog
    try:
        r = subprocess.run([sys.executable, "-c", _PICK_CODE, kind, initial or ""],
                           capture_output=True, text=True, timeout=300)
        return (r.stdout or "").strip()
    except Exception:  # noqa: BLE001
        return ""


# NB: module-level `st` is scipy.stats (line ~43); streamlit is imported locally,
# as in main(), so these callbacks touch st.session_state (streamlit), not scipy.
def _browse_outdir() -> None:
    import os
    import streamlit as st
    p = _pick_path("dir", st.session_state.get("out_dir", ""))
    if p:
        st.session_state["out_dir"] = os.path.normpath(p)


def _browse_clin() -> None:
    import os
    import streamlit as st
    p = _pick_path("file", st.session_state.get("clin_path", ""))
    if p:
        st.session_state["clin_path"] = os.path.normpath(p)


def _share_paths():
    """(repo_dir, .share_config.json, current_share_url.txt) for the shared-view launcher."""
    import os
    repo = os.path.dirname(os.path.abspath(__file__))
    return (repo, os.path.join(repo, ".share_config.json"),
            os.path.join(os.path.expanduser("~"), "Apps", "current_share_url.txt"))


def _start_share(cfg: dict) -> None:
    """Write the shared-view config and spawn the view-only app + Cloudflare tunnel (detached)."""
    import json
    import os
    import subprocess
    repo, cfg_path, url_path = _share_paths()
    with open(cfg_path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2)
    try:
        if os.path.exists(url_path):
            os.remove(url_path)
    except OSError:
        pass
    ps1 = os.path.join(repo, "share_tunnel.ps1")
    # CREATE_NO_WINDOW (0x08000000): runs hidden but reliably; DETACHED_PROCESS made
    # powershell exit before running the script.
    flags = (subprocess.CREATE_NEW_PROCESS_GROUP | 0x08000000) if os.name == "nt" else 0
    subprocess.Popen(
        ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", ps1],
        cwd=repo, creationflags=flags, close_fds=True,
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _read_share_url() -> str:
    """The current tunnel URL written by the launcher, or '' if not ready."""
    import os
    _, _, url_path = _share_paths()
    try:
        with open(url_path, encoding="utf-8-sig") as fh:
            u = fh.read().strip()
        return u if u.startswith("http") else ""
    except OSError:
        return ""


def _url_live(url: str) -> bool:
    """True if the public URL actually serves the app (tunnel edge is routable)."""
    import urllib.request
    try:
        with urllib.request.urlopen(url, timeout=12) as r:
            return r.status == 200
    except Exception:  # noqa: BLE001 - 5xx/DNS/TLS while the edge warms up
        return False


def _share_running() -> bool:
    """True if a shared view-only app appears to be listening on :8502."""
    import socket
    s = socket.socket()
    s.settimeout(0.3)
    try:
        return s.connect_ex(("127.0.0.1", 8502)) == 0
    except OSError:
        return False
    finally:
        s.close()


def _stop_share() -> None:
    """Take the shared link offline: stop only OUR cloudflared (the :8502 tunnel) and the
    view-only app on :8502, and remove the shared-view config so nothing lingers."""
    import os
    import subprocess
    # Kill only cloudflared processes whose command line targets :8502 (not unrelated ones),
    # then whatever is listening on :8502.
    ps = ("Get-CimInstance Win32_Process -Filter \"Name='cloudflared.exe'\" "
          "-ErrorAction SilentlyContinue | Where-Object { $_.CommandLine -match '8502' } | "
          "ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue };"
          "$c=Get-NetTCPConnection -LocalPort 8502 -State Listen -ErrorAction SilentlyContinue;"
          "if($c){$c.OwningProcess|Select-Object -Unique|"
          "ForEach-Object{Stop-Process -Id $_ -Force -ErrorAction SilentlyContinue}}")
    subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True)
    _, cfg_path, _ = _share_paths()
    try:
        if os.path.exists(cfg_path):
            os.remove(cfg_path)  # drop the config (contains the password) once sharing ends
    except OSError:
        pass


# ============================ Functional enrichment ============================
# Over-representation (GO / Reactome / KEGG) of the significant hits via g:Profiler,
# plus disease-association "novelty" via Open Targets. Both are OPT-IN network calls
# that send only gene SYMBOLS (never intensities, sample ids, or clinical metadata).
# Kept at module level (Streamlit-free) so they are unit-testable with a mocked
# transport; ``import requests`` stays lazy, matching the rest of the module.

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def _post_json(url: str, payload: dict, timeout: int = 30, _post=None) -> dict:
    """POST ``payload`` as JSON and return the parsed response.

    ``_post`` is the transport (defaults to ``requests.post``) and is injectable so
    tests can supply canned responses without touching the network.
    """
    if _post is None:
        import requests
        _post = requests.post
    resp = _post(url, json=payload, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _clean_symbols(genes) -> list[str]:
    """De-duplicate and stringify gene symbols, dropping empty/NaN (order-preserving).

    Delimited multi-gene fields (e.g. a protein group whose ``leading_gene_name`` is
    ``"GENEA;GENEB"``) are split on ``; , / |`` and whitespace before cleaning —
    real HGNC symbols never contain those characters, so a delimited value would
    otherwise be sent verbatim to g:Profiler/Open Targets, match nothing, and drop
    the feature from the analysis silently.
    """
    out: list[str] = []
    seen: set[str] = set()
    for g in genes:
        if g is None:
            continue
        for tok in re.split(r"[;,/|\s]+", str(g).strip()):
            s = tok.strip()
            if not s or s.lower() in {"nan", "none"} or s in seen:
                continue
            seen.add(s)
            out.append(s)
    return out


def sig_and_background_genes(res, annot, label_col, q_cut, lfc_cut, direction="both"):
    """Return ``(significant gene symbols, background = all TESTED gene symbols)``.

    Each tested feature id is mapped to a gene via ``annot[label_col]``. Using the
    tested set as the enrichment universe (not the whole genome) is essential in
    proteomics, where the *detectable* proteome is a biased slice — otherwise ORA
    just rediscovers "what is abundant in the sample." ``direction`` in
    {"both","up","down"} filters significant hits by ``logFC`` sign.
    """
    if label_col not in annot.columns:
        return [], []
    gene = annot[label_col]
    background = _clean_symbols(gene.reindex(res.index).tolist())
    mask = (res["adj.P.Val"] < q_cut) & (res["logFC"].abs() >= lfc_cut)
    if direction == "up":
        mask = mask & (res["logFC"] > 0)
    elif direction == "down":
        mask = mask & (res["logFC"] < 0)
    sig = _clean_symbols(gene.reindex(res.index[mask]).tolist())
    return sig, background


def gprofiler_enrichment(genes, background=None,
                         sources=("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"),
                         organism="hsapiens", user_threshold=0.05, _post=None):
    """Over-representation analysis via g:Profiler g:GOSt.

    Returns a tidy DataFrame (one row per g:SCS-significant term), or an empty
    DataFrame when nothing is enriched. When ``background`` is supplied it is used
    as a custom statistical domain (strongly recommended for proteomics).
    """
    query = _clean_symbols(genes)
    if not query:
        return pd.DataFrame()
    payload = {
        "organism": organism,
        "query": query,
        "sources": list(sources),
        "user_threshold": float(user_threshold),
        "significance_threshold_method": "g_SCS",
        "no_evidences": True,
    }
    bg = _clean_symbols(background or [])
    if bg:
        payload["background"] = bg
        payload["domain_scope"] = "custom_annotated"
    data = _post_json(GPROFILER_URL, payload, _post=_post)
    rows = data.get("result", []) if isinstance(data, dict) else []
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    rename = {"source": "source", "native": "term_id", "name": "term_name",
              "p_value": "p_value", "term_size": "term_size",
              "query_size": "query_size", "intersection_size": "intersection_size",
              "effective_domain_size": "domain_size"}
    have = {k: v for k, v in rename.items() if k in df.columns}
    out = df[list(have)].rename(columns=have)
    need = {"intersection_size", "query_size", "term_size", "domain_size"}
    if need <= set(out.columns):
        denom = (out["term_size"] / out["domain_size"]).replace(0, np.nan)
        out["fold_enrichment"] = (out["intersection_size"] / out["query_size"]) / denom
    return out.sort_values("p_value").reset_index(drop=True)


def opentargets_search_disease(name, _post=None):
    """Search Open Targets for a disease by name → list of ``(efo_id, label)``."""
    if not str(name).strip():
        return []
    q = ("query($q:String!){search(queryString:$q,entityNames:[\"disease\"],"
         "page:{index:0,size:10}){hits{id name entity}}}")
    data = _post_json(OPENTARGETS_URL, {"query": q, "variables": {"q": str(name)}}, _post=_post)
    hits = (((data or {}).get("data") or {}).get("search") or {}).get("hits") or []
    return [(h["id"], h.get("name", h["id"])) for h in hits if h.get("entity") == "disease"]


def opentargets_disease_targets(efo_id, size=3000, _post=None):
    """Map ``approvedSymbol`` (upper-cased) → association score for a disease.

    Pulls the top ``size`` associated targets (Open Targets returns them ordered by
    descending overall association score). NOTE: a gene ranked beyond ``size`` is
    absent from the returned map and therefore classified 'candidate-novel'; the
    default (3000) comfortably covers every non-trivially-associated gene for a
    typical disease, but a disease with more associations than ``size`` could
    mislabel a weakly-known gene. ``classify_novelty``'s threshold makes this moot
    for genes below the 'known' cutoff regardless.
    """
    q = ("query($efo:String!,$size:Int!){disease(efoId:$efo){associatedTargets("
         "page:{index:0,size:$size}){rows{score target{approvedSymbol}}}}}")
    data = _post_json(OPENTARGETS_URL,
                      {"query": q, "variables": {"efo": efo_id, "size": int(size)}}, _post=_post)
    rows = ((((data or {}).get("data") or {}).get("disease") or {})
            .get("associatedTargets") or {}).get("rows") or []
    scores: dict[str, float] = {}
    for r in rows:
        sym = ((r.get("target") or {}).get("approvedSymbol") or "").strip().upper()
        if sym:
            scores[sym] = float(r.get("score") or 0.0)
    return scores


def classify_novelty(sig_genes, assoc_scores, known_threshold=0.1):
    """Tag significant genes 'known' vs 'candidate-novel' by Open Targets score.

    A gene is 'known' when its association score is ``>= known_threshold``. Absence
    from ``assoc_scores`` (score ``None``/NaN) is treated as *not currently
    annotated in Open Targets* — a reflection of database coverage and recency, NOT
    proof of biological novelty.
    """
    rows = []
    for g in _clean_symbols(sig_genes):
        s = assoc_scores.get(g.upper())
        known = s is not None and s >= float(known_threshold)
        rows.append({"gene": g, "assoc_score": s,
                     "status": "known" if known else "candidate-novel"})
    df = pd.DataFrame(rows, columns=["gene", "assoc_score", "status"])
    if df.empty:
        return df
    return df.sort_values(["status", "assoc_score"], ascending=[True, False],
                          na_position="last").reset_index(drop=True)


def resolve_launch_prism_dir(env, query=None) -> str:
    """Resolve a pre-loaded PRISM output dir handed in at launch (the Skyline handoff).

    Precedence: the ``PDX_PRISM_DIR`` environment variable, then a ``prism_dir``
    query-string parameter (for a URL/tunnel launch). Returns a normalized path, or
    ``""`` when neither is set. Pure and injectable (``env`` and ``query`` are plain
    mappings) so it is unit-testable without a Streamlit runtime.

    Only consulted in normal (local) mode — a locked viewer/tunnel instance ignores
    it, so this never widens filesystem exposure on a shared deployment.
    """
    import os
    raw = str((env or {}).get("PDX_PRISM_DIR", "") or "").strip()
    if not raw and query:
        raw = str(query.get("prism_dir", "") or "").strip()
    return os.path.normpath(raw) if raw else ""


def main() -> None:
    import os

    import plotly.express as px
    import streamlit as st

    viewer = _viewer_config()
    _vtitle = (viewer or {}).get("display_name") if viewer else None
    st.set_page_config(page_title=_vtitle or "PRISM Differential Explorer", layout="wide")
    st.title(_vtitle or "PRISM Differential Explorer")

    # Fail CLOSED: a shared/tunnelled instance sets PDX_VIEWER=1. If the view config is
    # missing/corrupt (viewer is None), do NOT fall through to the normal browsable app on
    # a public endpoint — hard-stop instead, so a visitor can never reach the file controls.
    if viewer is None and os.environ.get("PDX_VIEWER"):
        st.error("This shared instance is misconfigured (no valid view configuration). "
                 "Refusing to start the browsable app on a public endpoint.")
        st.stop()

    if viewer is not None and not _viewer_gate(viewer):
        st.stop()  # locked until the shared-view password is entered

    if viewer is not None:
        st.caption("Shared, read-only view · explore the analysis below (data is fixed).")
    else:
        st.caption(
            "limma empirical-Bayes moderated t-test on PRISM corrected outputs. "
            "Join clinical metadata, define two groups, add covariates, explore hits."
        )

    cached_load = st.cache_data(load_prism)
    cached_attach = st.cache_data(attach_clinical)

    if viewer is not None:
        # Shared/tunnelled deployment: fixed data source, no path/Browse/Load controls
        # (a visitor must not be able to point the app at other files on the host).
        out_dir = str(viewer["out_dir"])
        clin_path = str(viewer.get("clin_path", "") or "")
        level = str(viewer.get("level", "protein"))
        load_clicked = False
        with st.sidebar:
            st.header("1 - Data")
            st.caption(f"📊 Shared dataset: **{viewer.get('display_name', 'PRISM analysis')}** "
                       f"· {level} level (read-only)")
        if st.session_state.get("prism_key") != (out_dir, level, clin_path):
            try:
                data = cached_load(out_dir, level)
                cinfo = None
                if clin_path.strip():
                    data = dict(data)
                    data["meta"], cinfo = cached_attach(data["meta"], clin_path.strip(),
                                                        viewer.get("clin_key"))
                st.session_state["prism"] = data
                st.session_state["prism_key"] = (out_dir, level, clin_path)
                st.session_state["clin_info"] = cinfo
                st.session_state.pop("analysis", None)
            except Exception as e:  # noqa: BLE001
                st.error(f"Shared dataset failed to load: {e}")
                st.stop()
        # Pin the shared comparison once so collaborators land on the analysis as set up
        # (reuses the proven ai_apply seeding path + a one-shot auto-run flag).
        comp = viewer.get("comparison")
        if comp and comp.get("grouping_column") and not st.session_state.get("_viewer_pinned"):
            st.session_state["_viewer_pinned"] = True
            st.session_state["ai_apply"] = {
                "grouping_column": comp.get("grouping_column"),
                "group_a_values": comp.get("group_a_values", []),
                "group_b_values": comp.get("group_b_values", []),
                "suggested_covariates": comp.get("suggested_covariates", []),
            }
            if viewer.get("q_cut") is not None:
                st.session_state["q_cut"] = float(viewer["q_cut"])
            if viewer.get("lfc_cut") is not None:
                st.session_state["lfc_cut"] = float(viewer["lfc_cut"])
            st.session_state["_viewer_autorun"] = True
    else:
        # Skyline handoff: if launched with a PRISM output dir (env var PDX_PRISM_DIR,
        # or a ?prism_dir=... query param), pre-fill it and auto-load once so the user
        # lands ready to add clinical metadata — no manual path entry needed.
        _launch_dir = resolve_launch_prism_dir(os.environ, dict(st.query_params))
        if _launch_dir and not st.session_state.get("_prism_launch_applied"):
            st.session_state["out_dir"] = _launch_dir
            st.session_state["_prism_launch_applied"] = True
            st.session_state["_prism_autoload"] = True
        with st.sidebar:
            st.header("1 - Data")
            if st.session_state.get("_prism_launch_applied"):
                st.caption("✓ PRISM output loaded from Skyline. Add your **clinical "
                           "metadata CSV** below to enable disease-vs-control analysis.")
            st.session_state.setdefault("out_dir", _DEFAULT_OUTDIR)
            st.session_state.setdefault("clin_path", _DEFAULT_CLINICAL)
            oc1, oc2 = st.columns([5, 1], vertical_alignment="bottom")
            out_dir = oc1.text_input("PRISM output_dir", key="out_dir")
            oc2.button("📁", key="browse_out", help="Browse for a folder", on_click=_browse_outdir)
            level = st.radio("Feature level", ["protein", "peptide"], horizontal=True)
            cc1, cc2 = st.columns([5, 1], vertical_alignment="bottom")
            clin_path = cc1.text_input("Clinical metadata CSV (optional)", key="clin_path")
            cc2.button("📁", key="browse_clin", help="Browse for a CSV file", on_click=_browse_clin)
            if clin_path.strip():
                try:
                    _clin_cols = list(pd.read_csv(clin_path.strip(), nrows=0).columns)
                except Exception:  # noqa: BLE001 - unreadable path; picker just stays on auto
                    _clin_cols = []
                if _clin_cols:
                    st.selectbox(
                        "Clinical ID column", ["(auto-detect)"] + _clin_cols, key="clin_key",
                        help="Which CSV column holds the sample identifier that matches your "
                             "PRISM sample names. Leave on auto unless it picks the wrong column.")
            load_clicked = st.button("Load", type="primary")
            # Skyline handoff: auto-load once without requiring a manual Load click.
            if st.session_state.pop("_prism_autoload", False):
                load_clicked = True

    clin_key = (None if st.session_state.get("clin_key", "(auto-detect)") == "(auto-detect)"
                else st.session_state.get("clin_key"))
    if load_clicked:
        try:
            data = cached_load(out_dir, level)
            cinfo = None
            if clin_path.strip():
                data = dict(data)
                data["meta"], cinfo = cached_attach(data["meta"], clin_path.strip(), clin_key)
            st.session_state["prism"] = data
            st.session_state["prism_key"] = (out_dir, level, clin_path)
            st.session_state["clin_info"] = cinfo
            st.session_state.pop("analysis", None)  # drop results from a prior dataset
        except Exception as e:  # noqa: BLE001
            st.error(f"Load failed: {e}")
            st.stop()

    if "prism" not in st.session_state:
        st.info("Set paths in the sidebar and click **Load**.")
        st.stop()

    data = st.session_state["prism"]
    meta: pd.DataFrame = data["meta"]
    expr_log2: pd.DataFrame = data["expr_log2"]
    annot: pd.DataFrame = data["annot"]
    label_col: str = data["label_col"]
    cinfo = st.session_state.get("clin_info")

    name_to_id = dict(zip(meta["sample"], meta["sample_id"]))
    id_to_name = dict(zip(meta["sample_id"], meta["sample"]))

    n_feat, n_samp = expr_log2.shape
    msg = f"Loaded **{data['level']}**: {n_feat:,} features x {n_samp} samples."
    st.success(msg)
    if st.session_state.get("prism_key") != (out_dir, level, clin_path):
        st.warning(f"Sidebar settings changed (you're viewing **{data['level']}** level) — "
                   "click **Load** to apply them.", icon="⚠️")
    if cinfo is not None:
        how = "manually chosen" if cinfo.get("manual") else "auto-detected"
        if cinfo["key_column"] and cinfo["match_rate"] > 0:
            st.caption(f"Clinical join: matched **{cinfo['match_rate']:.0%}** of samples via "
                       f"column **`{cinfo['key_column']}`** ({how}) — "
                       f"added {cinfo['n_clinical_cols']} metadata columns.")
            if cinfo["match_rate"] < 0.9:
                st.caption("⚠️ Low match rate — if it picked the wrong column, set "
                           "**Clinical ID column** in the sidebar and click **Load** again.")
        else:
            st.warning("Clinical CSV loaded but **no column matched the sample names**. "
                       "Pick the correct **Clinical ID column** in the sidebar (the one holding "
                       "sample/replicate identifiers) and click **Load** again.")

    # metadata columns available for grouping / covariates
    meta_cols = [c for c in meta.columns if c not in _NON_META_COLS and meta[c].notna().any()]
    # Grouping ("Group by") categories = text columns + low-cardinality numeric columns
    # (discrete codes like stage/grade/batch), whose values are cast to strings when
    # grouping. Continuous numerics stay numeric-only (covariate / longitudinal time axis).
    cat_cols = _grouping_columns(meta, meta_cols)

    exp_pool = meta[meta["sample_type"] == "experimental"] if "sample_type" in meta.columns else meta

    def _validate_suggestion(s: dict) -> dict | None:
        """Normalise a suggestion (heuristic or LLM) against the live experimental data.

        Drops out-of-vocabulary values/columns and fills missing fields so the UI
        is safe even if a local model omits some keys.
        """
        col = s.get("grouping_column")
        if col not in cat_cols:
            return None
        present = set(exp_pool[col].dropna().astype(str))
        a = [v for v in s.get("group_a_values", []) if v in present]
        b = [v for v in s.get("group_b_values", []) if v in present and v not in a]
        if not a or not b:
            return None
        cov = [c for c in s.get("suggested_covariates", []) if c in meta_cols and c != col]
        la = str(s.get("label_a") or "Group A")
        lb = str(s.get("label_b") or "Group B")
        return {
            "name": str(s.get("name") or f"{lb} vs {la} ({col})"),
            "grouping_column": col, "group_a_values": a, "group_b_values": b,
            "label_a": la, "label_b": lb, "suggested_covariates": cov,
            "rationale": str(s.get("rationale") or ""),
            "confound_warning": str(s.get("confound_warning") or ""),
        }

    # Apply a pending AI suggestion by seeding widget state BEFORE the widgets render.
    pend = st.session_state.pop("ai_apply", None)
    if pend and cat_cols:
        st.session_state["grp_mode"] = "By metadata column"
        st.session_state["grp_col"] = pend["grouping_column"]
        st.session_state["grp_a_vals"] = pend["group_a_values"]
        st.session_state["grp_b_vals"] = pend["group_b_values"]
        st.session_state["cov_select"] = pend["suggested_covariates"]
        # remember if this contrast came from a Finder search (selection bias)
        st.session_state["search_contrast"] = pend.get("_from_search")

    with st.sidebar:
        st.header("2 - Define groups")
        types = sorted(meta["sample_type"].dropna().unique().tolist())
        keep_types = st.multiselect("Include sample types", types, default=["experimental"])
        pool_ids = meta[meta["sample_type"].isin(keep_types)]["sample_id"].tolist()

        if cat_cols:
            oll_running, oll_models = ollama_status()
            engine = "Heuristic"
            oll_model = None
            if oll_running and oll_models:
                engine = st.radio("Suggestion engine",
                                  ["Heuristic (instant)", "Local LLM (Ollama)"],
                                  horizontal=True, key="sug_engine")
                if engine.startswith("Local"):
                    oll_model = st.selectbox("Ollama model", oll_models, key="oll_model")
            elif oll_running:
                st.caption("Ollama is running but no model is pulled "
                           "(`ollama pull llama3.2:3b`). Using the heuristic.")
            if st.button("✨ Suggest comparisons", help="Reads your metadata and proposes "
                         "control-vs-case contrasts. Fully local — no cloud, no cost."):
                use_llm = engine.startswith("Local") and oll_model
                try:
                    if use_llm:
                        with st.spinner(f"Asking {oll_model} (local, ~1-2 min on CPU)…"):
                            raw = suggest_groups_ollama(meta, oll_model)
                    else:
                        raw = suggest_groups(meta)
                except Exception as e:  # noqa: BLE001
                    st.warning(f"Ollama suggest failed ({e}); using the heuristic.")
                    raw = suggest_groups(meta)
                sug = [v for v in (_validate_suggestion(s) for s in raw) if v]
                st.session_state["ai_suggestions"] = sug
                if not sug:
                    st.info("No usable comparisons found — define groups manually below.")

        modes = ["By metadata column", "Manual / pattern"] if cat_cols else ["Manual / pattern"]
        mode = st.radio("Grouping mode", modes, horizontal=True, key="grp_mode")

        def _derive(vals):
            return " + ".join(vals) if vals else ""

        group_a, group_b = [], []
        default_la = default_lb = ""
        if mode == "By metadata column":
            default_col = "Condition" if "Condition" in cat_cols else cat_cols[0]
            gcol = st.selectbox("Group by", cat_cols, index=cat_cols.index(default_col),
                                key="grp_col")
            gser = meta.loc[pool_ids, gcol]
            vals = [v for v in _as_cat_str(gser.dropna()).unique().tolist() if v != ""]
            try:  # numeric codes sort by value ("2" before "10"), text sorts lexically
                vals = sorted(vals, key=float)
            except (TypeError, ValueError):
                vals = sorted(vals)
            # Drop any seeded/stale selections that aren't valid for this column.
            if "grp_a_vals" in st.session_state:
                st.session_state["grp_a_vals"] = [v for v in st.session_state["grp_a_vals"] if v in vals]
            a_vals = st.multiselect("Group A values", vals, key="grp_a_vals")
            b_opts = [v for v in vals if v not in a_vals]
            if "grp_b_vals" in st.session_state:
                st.session_state["grp_b_vals"] = [v for v in st.session_state["grp_b_vals"] if v in b_opts]
            b_vals = st.multiselect("Group B values", b_opts, key="grp_b_vals")
            colvals = _as_cat_str(gser)
            group_a = [i for i in pool_ids if colvals.get(i) in a_vals]
            group_b = [i for i in pool_ids if colvals.get(i) in b_vals]
            default_la, default_lb = _derive(a_vals), _derive(b_vals)
        else:
            pool_names = [id_to_name[i] for i in pool_ids]
            pat_a = st.text_input("Group A name contains", value="")
            pat_b = st.text_input("Group B name contains", value="")
            da = [s for s in pool_names if pat_a and pat_a.lower() in s.lower()]
            db = [s for s in pool_names if pat_b and pat_b.lower() in s.lower() and s not in da]
            an = st.multiselect("Group A samples", pool_names, default=da)
            bn = st.multiselect("Group B samples", [s for s in pool_names if s not in an], default=db)
            group_a = [name_to_id[s] for s in an]
            group_b = [name_to_id[s] for s in bn]
            default_la, default_lb = pat_a.strip(), pat_b.strip()

        # Editable labels used everywhere in place of "A"/"B". Default to the
        # selected values; shorten as you like (e.g. "Sporadic ADD + ..." -> "ADD").
        label_a = st.text_input("Label for Group A", value=default_la or "Group A").strip() or "Group A"
        label_b = st.text_input("Label for Group B", value=default_lb or "Group B").strip() or "Group B"
        st.caption(f"**{label_a}**: n={len(group_a)}  |  **{label_b}**: n={len(group_b)}")

        st.header("3 - Model covariates")
        cov_choices = [c for c in meta_cols if mode == "Manual / pattern" or c != gcol]
        if "cov_select" in st.session_state:
            st.session_state["cov_select"] = [c for c in st.session_state["cov_select"] if c in cov_choices]
        chosen_cov = st.multiselect(
            "Covariates (limma design)", cov_choices, key="cov_select",
            help="Categorical -> dummies, numeric -> centred. Confounded/constant cols are auto-dropped.",
        )

        st.header("4 - Thresholds")
        sig_basis_label = st.radio(
            "Significance by", ["FDR (BH-adjusted)", "Raw p-value"], key="sig_basis",
            help="FDR controls the false-discovery rate across all tested features "
                 "(rigorous, conservative). Raw p-value is uncorrected (exploratory) — use it "
                 "to reproduce a raw p<0.05 call, but it does not control multiple testing.")
        sig_basis = "raw" if sig_basis_label.startswith("Raw") else "fdr"
        _cut_label = "p-value cutoff" if sig_basis == "raw" else "FDR (adj.P.Val) cutoff"
        q_cut = st.slider(_cut_label, 0.0, 0.25, 0.05, 0.01, key="q_cut")
        lfc_cut = st.slider("|log2 FC| cutoff", 0.0, 3.0, 1.0, 0.1, key="lfc_cut")
        run_clicked = st.button("Run limma", type="primary")
    # Shared view: auto-run the pinned comparison once so it lands on results.
    run_clicked = run_clicked or st.session_state.pop("_viewer_autorun", False)

    # AI-suggested comparisons (rendered in the main area; click to load into the sidebar).
    if st.session_state.get("ai_suggestions"):
        with st.expander(f"✨ {len(st.session_state['ai_suggestions'])} suggested comparisons "
                         "— click Use to load one into the sidebar", expanded=True):
            for i, s in enumerate(st.session_state["ai_suggestions"]):
                c1, c2 = st.columns([5, 1])
                with c1:
                    st.markdown(
                        f"**{s['name']}** — `{s['grouping_column']}`: "
                        f"**{s['label_b']}** ({', '.join(s['group_b_values'])}) vs "
                        f"**{s['label_a']}** ({', '.join(s['group_a_values'])})"
                    )
                    st.caption(s["rationale"]
                               + (f"  ⚠️ {s['confound_warning']}" if s.get("confound_warning") else "")
                               + (f"  · covariates: {', '.join(s['suggested_covariates'])}"
                                  if s.get("suggested_covariates") else ""))
                if c2.button("Use", key=f"use_sug_{i}"):
                    st.session_state["ai_apply"] = s
                    st.rerun()

    # Run the fit only when the button is pressed, then CACHE it (compute runs
    # regardless of the active tab, so results survive tab switches / reruns and
    # the PCA tab works without a fit).
    if run_clicked:
        cov_df = meta[chosen_cov] if chosen_cov else None
        try:
            res_fit, rinfo = differential(expr_log2, group_a, group_b, covariates=cov_df)
        except Exception as e:  # noqa: BLE001
            st.error(f"Analysis failed: {e}")
        else:
            res_fit = res_fit.join(annot[[c for c in [label_col] if c in annot.columns]], how="left")
            # is this exactly the contrast a Finder search selected? (selection bias)
            sc = st.session_state.get("search_contrast")
            from_search = bool(
                sc and mode == "By metadata column" and sc.get("col") == gcol
                and sorted(a_vals) == sc.get("a") and sorted(b_vals) == sc.get("b"))
            st.session_state["analysis"] = {
                "res": res_fit, "rinfo": rinfo,
                "group_a": group_a, "group_b": group_b, "label_col": label_col,
                "label_a": label_a, "label_b": label_b,
                "from_search": from_search,
                "n_searched": sc.get("n_tested") if (from_search and sc) else None,
            }

    def detect_cryptic(term: str) -> tuple[list, dict]:
        """Find cryptic feature ids (annotation match + peptide-level merged_data)."""
        mask = pd.Series(False, index=expr_log2.index)
        if term:
            for c in annot.columns:
                col = annot[c]
                if col.dtype == object or str(col.dtype).startswith("string"):
                    hit = col.astype(str).str.contains(re.escape(term), case=False, na=False)
                    mask = mask | hit.reindex(mask.index, fill_value=False)
            idx_hit = pd.Series(expr_log2.index.astype(str), index=expr_log2.index) \
                .str.contains(re.escape(term), case=False, na=False)
            mask = mask | idx_hit
        cids = list(expr_log2.index[mask])
        cmap: dict = {}
        if data["level"] == "peptide" and term:
            mp = Path(st.session_state.get("prism_key", ("",))[0] or "") / "merged_data.parquet"
            if mp.exists():
                try:
                    cmap = st.cache_data(cryptic_peptide_map)(str(mp), term)
                    cids = sorted(set(cids) | {p for p in expr_log2.index if p in cmap})
                except Exception:  # noqa: BLE001
                    pass
        return cids, cmap

    def render_finder(mat: pd.DataFrame, key: str, note: str = "",
                      unit: str = "features") -> None:
        """Significance-search UI over `mat` (all features, or a cryptic subset).

        `unit` is the noun used in the summary (e.g. "cryptic peptides") so the
        cryptic finder states significant *cryptic peptides*, not generic features.
        """
        if not cat_cols:
            st.info("No categorical metadata columns to scan. Load clinical metadata first.")
            return
        nlev = {c: meta[c].dropna().nunique() for c in cat_cols}
        scan_default = max(cat_cols, key=lambda c: nlev[c])
        f1, f2, f3 = st.columns([2, 2, 1])
        scan_col = f1.selectbox("Scan column", cat_cols,
                                index=cat_cols.index(scan_default), key=f"{key}_col")
        scope = f2.selectbox("Search", ["One-vs-rest", "Pairwise", "Subset pairs (thorough)"],
                             key=f"{key}_scope")
        min_n = f3.slider("Min group n", 2, 10, 3, key=f"{key}_minn")
        n_levels = (meta[meta.get("sample_type", "") == "experimental"][scan_col].dropna().nunique()
                    if "sample_type" in meta.columns else meta[scan_col].dropna().nunique())
        max_tests = 400
        if scope.startswith("Subset"):
            max_tests = st.slider("Max splits to test (subset search grows fast)",
                                  50, 2000, 400, 50, key=f"{key}_max")
        cov_df_scan = meta[chosen_cov] if chosen_cov else None
        st.caption(f"`{scan_col}` has **{n_levels}** levels (≥{min_n} samples). Scanning "
                   f"**{len(mat):,}** {unit}" + (f" ({note})" if note else "")
                   + f" at FDR<{q_cut:g}, |log2FC|≥{lfc_cut:g}"
                   + (f", covariates: {', '.join(chosen_cov)}" if chosen_cov else "") + ".")
        st.warning("Exploratory — searching many contrasts inflates false positives. The top hit "
                   "is permutation-calibrated below; still validate any winner independently.",
                   icon="⚠️")
        if st.button("Run significance scan", type="primary", key=f"{key}_run"):
            with st.spinner(f"Scanning {scope.lower()} splits of {scan_col}…"):
                df, trunc = significance_scan(mat, meta, scan_col, scope, covariates=cov_df_scan,
                                              q_cut=q_cut, lfc_cut=lfc_cut, min_n=min_n,
                                              max_tests=max_tests)
            st.session_state[f"{key}_result"] = {
                "df": df, "trunc": trunc, "col": scan_col, "q": q_cut, "lfc": lfc_cut,
                "cov": list(chosen_cov), "min_n": min_n, "scope": scope, "max_tests": max_tests}
            st.session_state.pop(f"{key}_null", None)  # stale calibration

        sr = st.session_state.get(f"{key}_result")
        if sr is not None and len(sr["df"]):
            df = sr["df"]
            s = scan_summary(df)
            st.write(f"**{s['n_contrasts']} contrasts tested** on `{sr['col']}`"
                     + ("  (truncated at the split cap)" if sr["trunc"] else "")
                     + f"  ·  best contrast (**{s['best_B']}** vs **{s['best_A']}**): "
                     + f"**{s['best_n_sig']}** significant {unit}.")
            st.caption(f"{s['n_with_hits']} of {s['n_contrasts']} contrasts had ≥1 significant "
                       f"{unit} at FDR<{sr['q']:g}, |log2FC|≥{sr['lfc']:g}. Per-contrast counts "
                       "are in the **n_sig** column below.")
            top = df.iloc[0]
            pool = (meta[meta.get("sample_type", "") == "experimental"]
                    if "sample_type" in meta.columns else meta)
            cv = pool[sr["col"]].astype(str)
            ga = [s for lv in top["_a"] for s in pool.index[cv == lv]]
            gb = [s for lv in top["_b"] for s in pool.index[cv == lv]]
            cov_df = meta[sr["cov"]] if sr["cov"] else None
            obs = int(top["n_sig"])
            cc1, cc2 = st.columns([2, 1])
            familywise = cc1.checkbox(
                "Family-wise calibration (re-runs the whole search under permutations — "
                "the statistically correct null; slower)", value=False, key=f"{key}_fw")
            n_perm = cc2.slider("Permutations", 100, 1000, 200, 100, key=f"{key}_nperm")
            # Permutation is heavy (re-runs the model many times) — compute ONLY on click,
            # then cache, so it doesn't re-run on every Streamlit interaction.
            if st.button("Calibrate top hit (permutation test)", key=f"{key}_calib"):
                if familywise:
                    with st.spinner(f"Family-wise calibration: re-running the search ×{n_perm}…"):
                        null = scan_permutation_null(mat, meta, sr["col"], sr["scope"], cov_df,
                                                     sr["q"], sr["lfc"], sr["min_n"],
                                                     sr.get("max_tests", 400), n_perm=n_perm)
                    null_label = "max n_sig over the search, permuted"
                else:
                    with st.spinner(f"Calibrating the top split ×{n_perm}…"):
                        null = permute_calibrate(mat, ga, gb, cov_df, sr["q"], sr["lfc"],
                                                 n_perm=n_perm, min_n=sr["min_n"])
                    null_label = "n_sig for this split, permuted labels"
                st.session_state[f"{key}_null"] = {"null": null, "familywise": familywise,
                                                   "obs": obs, "label": null_label}
            cal = st.session_state.get(f"{key}_null")
            if not (cal and cal["null"]):
                st.caption("Click **Calibrate** to permutation-test whether the top hit is beyond "
                           "chance (turn on family-wise for the search-corrected version).")
            else:
                null, familywise, obs = cal["null"], cal["familywise"], cal["obs"]
                null_label = cal["label"]
                null_arr = np.array(null)
                emp_p = (int((null_arr >= obs).sum()) + 1) / (len(null_arr) + 1)
                c1, c2, c3 = st.columns(3)
                c1.metric(f"Best contrast: significant {unit}", obs)
                c2.metric("Null n_sig", f"{np.median(null_arr):.0f} "
                          f"(95th: {np.percentile(null_arr, 95):.0f})")
                c3.metric("Empirical p", f"{emp_p:.3f}"
                          + ("" if familywise else "  (within-split)"))
                st.caption(f"Null = {null_label}; min achievable p = {1 / (len(null_arr) + 1):.3f}."
                           + ("" if familywise else "  This does **not** correct for the search over "
                              "partitions — turn on family-wise calibration for that."))
                if emp_p > 0.1:
                    st.error(f"The top hit ({obs}) is **not beyond chance** (p={emp_p:.2f}; null "
                             f"median {np.median(null_arr):.0f}). Likely a search/multiplicity "
                             "artifact.")
                elif familywise:
                    st.success(f"The top hit ({obs}) survives **family-wise** permutation "
                               f"(p={emp_p:.3f}) — a real, selection-corrected signal.")
                else:
                    st.info(f"The top hit ({obs}) beats permuted labels for this split "
                            f"(p={emp_p:.3f}), but this isn't search-corrected — confirm with "
                            "family-wise calibration above.")
            show = df[["A", "B", "n_A", "n_B", "n_sig", "min_q"]].head(30)
            st.dataframe(show.style.format({"min_q": "{:.2e}"}), width="stretch", height=360)
            st.download_button("Download full scan (CSV)",
                               df.drop(columns=["_a", "_b"]).to_csv(index=False).encode(),
                               file_name=f"significance_scan_{sr['col']}.csv", key=f"{key}_dl")
            st.caption("Load a contrast into the **Differential** tab to inspect it:")
            for i in range(min(5, len(df))):
                r = df.iloc[i]
                if st.button(f"Use #{i + 1}:  {r['B']}  vs  {r['A']}  (n_sig={int(r['n_sig'])})",
                             key=f"{key}_use_{i}"):
                    st.session_state["ai_apply"] = {
                        "name": f"{r['B']} vs {r['A']}", "grouping_column": sr["col"],
                        "group_a_values": list(r["_a"]), "group_b_values": list(r["_b"]),
                        "label_a": "A", "label_b": "B", "suggested_covariates": sr["cov"],
                        "rationale": "", "confound_warning": "",
                        "_from_search": {"col": sr["col"], "a": sorted(r["_a"]),
                                         "b": sorted(r["_b"]), "n_tested": int(len(df))}}
                    st.rerun()
        elif sr is not None:
            st.info("No contrasts met the minimum group size — lower **Min group n**.")

    (tab_diff, tab_pca, tab_ev, tab_find, tab_cfind, tab_cryptic,
     tab_enrich, tab_long, tab_het, tab_report) = st.tabs(
        ["Differential", "PCA", "EV markers", "Finder", "Cryptic finder", "Cryptic",
         "Enrichment", "Longitudinal", "Heterogeneity", "Report"])

    # ============================ Differential tab ============================
    with tab_diff:
        if group_a and group_b:
            with st.expander("Group balance vs metadata (check confounding)", expanded=True):
                check_cols = [c for c in (chosen_cov + ["batch"]) if c in cat_cols]
                shown = False
                for c in dict.fromkeys(check_cols):
                    ct = covariate_balance(meta, group_a, group_b, c)
                    if ct is not None and ct.shape[0] > 1:
                        ct = ct.rename(columns={"A": label_a, "B": label_b})
                        st.write(f"**{c}**")
                        st.dataframe(ct)
                        if ((ct == 0).any(axis=1) & (ct.sum(axis=1) > 0)).any():
                            st.warning(
                                f"A level of '{c}' falls entirely in one group - partially "
                                f"confounded. Add '{c}' as a covariate and interpret with care."
                            )
                        shown = True
                if not shown:
                    st.caption("No categorical covariate selected to cross-tabulate yet.")

        st.info(
            "The corrected matrix is dense, so this tests **integrated-window intensity** "
            "(detection conflated with abundance). For true on/off detection, recover it from "
            "transition-level `merged_data`.",
            icon="ℹ️",
        )

        if "analysis" not in st.session_state:
            st.info("Set up the groups, then click **Run limma**.")
        else:
            A = st.session_state["analysis"]
            res = A["res"].copy()
            rinfo = A["rinfo"]
            res_group_a, res_group_b = A["group_a"], A["group_b"]
            res_label_col = A["label_col"]
            la, lb = A.get("label_a", "Group A"), A.get("label_b", "Group B")
            mean_a_col, mean_b_col = f"mean[{la}]", f"mean[{lb}]"
            res = res.rename(columns={"mean_A": mean_a_col, "mean_B": mean_b_col})

            st.caption(f"Contrast: **{lb}** vs **{la}**  (log2FC > 0 = higher in {lb})")
            if data["level"] == "peptide":
                st.caption("ℹ️ Peptide-level FDR treats peptides as independent — peptides from "
                           "one protein are correlated, so these q-values are anti-conservative "
                           "for *protein*-level claims. Use protein level as the primary readout "
                           "(peptide/cryptic as targeted follow-up).")
            if A.get("from_search"):
                st.warning(
                    f"⚠️ **This contrast was selected from a {A.get('n_searched', 'multi')}-way "
                    "Finder search** — it was chosen *because* it maximized the hit count. The "
                    "q-values below are **not adjusted for that selection** and will be "
                    "anti-conservative. Treat as hypothesis generation; confirm on independent "
                    "data or a pre-specified contrast. (The Finder's permutation p is the honest "
                    "family-wise signal.)", icon="⚠️")
            confound_msgs = [m for m in rinfo["messages"] if "confounded with group" in m]
            other_msgs = [m for m in rinfo["messages"] if "confounded with group" not in m]
            if confound_msgs:
                st.warning("⚠️ A requested covariate could not be separated from the group and "
                           "was **dropped** — the contrast is partly confounded and the adjusted "
                           "effect is not identifiable along that axis:\n\n- "
                           + "\n- ".join(confound_msgs))
            for msg in other_msgs:
                st.caption("• " + msg)
            if rinfo.get("covariates_used"):
                st.caption("Model covariates used: " + ", ".join(rinfo["covariates_used"]))

            sig_series, p_thresh = _significance(res, sig_basis, q_cut, lfc_cut)
            res["significant"] = sig_series
            _basis_word = "raw p" if sig_basis == "raw" else "FDR"

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Features tested", f"{rinfo['n_features_tested']:,}")
            c2.metric("Dropped (incomplete)", f"{rinfo['n_features_dropped']:,}")
            c3.metric(f"Significant ({_basis_word}<{q_cut:g})", f"{int(res['significant'].sum()):,}")
            c4.metric("Median df (prior+resid)", f"{rinfo['df_prior']:.1f} + {rinfo['df_residual']:.0f}")

            label = res[res_label_col] if res_label_col in res.columns else res.index.to_series()
            res = res.assign(
                _label=label.astype(str),
                _neglog10p=-np.log10(res["P.Value"].clip(lower=1e-300)),
                _dir=np.where(res["significant"],
                              np.where(res["logFC"] > 0, f"up in {lb}", f"up in {la}"), "ns"),
            )
            fig = px.scatter(
                res, x="logFC", y="_neglog10p", color="_dir",
                color_discrete_map={f"up in {lb}": "#d62728", f"up in {la}": "#1f77b4", "ns": "#b0b0b0"},
                hover_name="_label",
                hover_data={"adj.P.Val": ":.3g", "logFC": ":.2f", "_neglog10p": False, "_dir": False},
                labels={"logFC": f"log2 fold change ({lb} / {la})", "_neglog10p": "-log10 P"}, height=520,
            )
            _add_volcano_guides(fig, p_thresh, lfc_cut, sig_basis, q_cut)
            if p_thresh is None:
                st.caption(f"No features pass FDR {q_cut:g} — no cutoff line drawn. "
                           "(Try 'Raw p-value' in the sidebar for an exploratory view.)")
            fig.update_layout(legend_title_text="", margin=dict(l=10, r=10, t=30, b=10))
            st.plotly_chart(fig, width="stretch")

            # C4/M1: per-group detection % (peptide level) — abundance on the dense matrix
            # conflates detection with abundance; flag features that are mostly borrowed-baseline.
            det_extra = []
            mp_d = Path(st.session_state.get("prism_key", ("",))[0] or "") / "merged_data.parquet"
            if data["level"] == "peptide" and mp_d.exists():
                if st.checkbox("Annotate per-group detection % (flags borrowed-baseline features)",
                               value=False, key="diff_detann",
                               help="For each feature, the fraction of samples genuinely detected "
                                    "(DetectionQValue<0.01). Abundance is unreliable where a group "
                                    "is mostly undetected — the value is integrated baseline noise."):
                    try:
                        dm = st.cache_data(cryptic_detection_matrix)(str(mp_d), None)
                        ac = [s for s in res_group_a if s in dm.columns]
                        bc = [s for s in res_group_b if s in dm.columns]
                        idx = res.index.intersection(dm.index)
                        res.loc[idx, "det%_A"] = dm.loc[idx, ac].mean(axis=1) * 100
                        res.loc[idx, "det%_B"] = dm.loc[idx, bc].mean(axis=1) * 100
                        res["abundance_unreliable"] = (res[["det%_A", "det%_B"]].min(axis=1) < 50)
                        det_extra = ["det%_A", "det%_B", "abundance_unreliable"]
                        n_unrel = int((res["significant"] & res["abundance_unreliable"]).sum())
                        if n_unrel:
                            st.warning(f"⚠️ {n_unrel} of your 'significant' features are <50% detected "
                                       "in a group — their abundance logFC reflects borrowed baseline, "
                                       "not real signal. Prefer the **detection** test for those.")
                    except Exception as e:  # noqa: BLE001
                        st.caption(f"Detection annotation unavailable: {e}")

            show_cols = [c for c in [res_label_col, "logFC", "FC", "AveExpr", "t", "P.Value",
                                     "adj.P.Val", mean_a_col, mean_b_col] + det_extra
                         + ["significant"] if c in res.columns]
            only_sig = st.checkbox("Show significant only", value=True)
            table = res[res["significant"]] if only_sig else res
            st.dataframe(
                table[show_cols].style.format(
                    {"logFC": "{:.2f}", "FC": "{:.2f}", "AveExpr": "{:.2f}", "t": "{:.2f}",
                     "P.Value": "{:.2e}", "adj.P.Val": "{:.2e}", mean_a_col: "{:.2f}",
                     mean_b_col: "{:.2f}", "det%_A": "{:.0f}", "det%_B": "{:.0f}"}
                ),
                width="stretch", height=320,
            )
            d1, d2 = st.columns(2)
            _tag = "_SEARCH-SELECTED_not-selection-adjusted" if A.get("from_search") else ""
            d1.download_button("Download all results (CSV)", res[show_cols].to_csv().encode(),
                               file_name=f"prism_diff_results{_tag}.csv")
            d2.download_button("Download significant (CSV)", res[res["significant"]][show_cols].to_csv().encode(),
                               file_name=f"prism_diff_significant{_tag}.csv")

            st.subheader("Feature detail")
            options = table.index.tolist()
            if options:
                pick = st.selectbox(
                    "Feature", options,
                    format_func=lambda i: f"{res.loc[i, res_label_col]}  ({i})" if res_label_col in res.columns else str(i),
                )
                vals = expr_log2.loc[pick, res_group_a + res_group_b]
                long = pd.DataFrame({
                    "log2 abundance": vals.values,
                    "group": [la] * len(res_group_a) + [lb] * len(res_group_b),
                    "sample": [id_to_name.get(s, s) for s in (res_group_a + res_group_b)],
                })
                bx = px.box(long, x="group", y="log2 abundance", points="all", color="group",
                            category_orders={"group": [la, lb]},
                            hover_data=["sample"], color_discrete_map={la: "#1f77b4", lb: "#d62728"}, height=380)
                r = res.loc[pick]
                bx.update_layout(showlegend=False, margin=dict(l=10, r=10, t=40, b=10),
                                 title=f"{r.get(res_label_col, pick)}  |  log2FC={r['logFC']:.2f}  q={r['adj.P.Val']:.2e}")
                st.plotly_chart(bx, width="stretch")

    # ================================ PCA tab ================================
    with tab_pca:
        st.caption("PCA over the selected samples (complete-case features, feature-centred log2). "
                   "Independent of the limma run.")
        pca_types = st.multiselect("Samples to include (by type)", types, default=types, key="pca_types")
        pca_ids = meta[meta["sample_type"].isin(pca_types)]["sample_id"].tolist()

        # Which features to run PCA on. "All complete" is unsupervised over the whole
        # proteome (global variance dominates -> disease often not on PC1/2). "Top
        # variable" focuses on the most dynamic features. "Significant hits" projects
        # onto the limma-significant set -> groups WILL separate, but that is circular
        # (features were chosen *because* they differ by group) so it is illustrative,
        # not evidence of separation.
        feat_choices = ["All complete", "Top variable", "Significant hits"]
        if data["level"] == "peptide":
            feat_choices.append("Cryptic peptides")
        feat_mode = st.radio("Features", feat_choices, horizontal=True, key="pca_feat")
        pca_expr = expr_log2
        if feat_mode == "Top variable":
            topn = st.slider("N most variable features", 50, 3000, 500, 50, key="pca_topn")
            sub = expr_log2[pca_ids].dropna(how="any")
            pca_expr = expr_log2.loc[sub.var(axis=1).nlargest(topn).index]
        elif feat_mode == "Cryptic peptides":
            cterm_pca = st.text_input("Accession contains", value="cryptic", key="pca_crypt_term").strip()
            cids_pca, _ = detect_cryptic(cterm_pca)
            pca_expr = expr_log2.loc[[p for p in cids_pca if p in expr_log2.index]]
            st.caption(f"PCA on **{pca_expr.shape[0]}** cryptic peptides. Note: cryptic peptides "
                       "are low-abundance/on-off, so the dense Area matrix here is largely "
                       "borrowed baseline — PC structure may reflect detection/noise, not biology.")
        elif feat_mode == "Significant hits":
            if "analysis" in st.session_state:
                rf = st.session_state["analysis"]["res"]
                sig_idx = rf.index[(rf["adj.P.Val"] < q_cut) & (rf["logFC"].abs() >= lfc_cut)]
                pca_expr = expr_log2.loc[expr_log2.index.intersection(sig_idx)]
                st.caption(
                    f"Projecting onto {len(pca_expr):,} significant features from the current run. "
                    "Circular by construction (features selected because they differ by group) - "
                    "use to visualise the hits, not to claim separation."
                )
            else:
                st.warning("Run limma first to use significant hits - showing all features instead.")

        if len(pca_ids) < 3:
            st.warning("Select at least 3 samples for PCA.")
        elif pca_expr.shape[0] < 2:
            st.warning("Not enough features in the selected set for PCA.")
        else:
            try:
                scores, var_ratio, n_used = st.cache_data(compute_pca)(pca_expr, pca_ids)
            except Exception as e:  # noqa: BLE001
                st.warning(str(e))
            else:
                # warn when the complete-case feature set is a small/biased slice
                total_feat = int(pca_expr.shape[0])
                if total_feat and n_used / total_feat < 0.5:
                    st.warning(f"PCA uses only **{n_used:,} of {total_feat:,}** features "
                               f"({n_used / total_feat:.0%}) measured in *every* included sample. "
                               "Missingness is non-random (low-abundance, batch-linked), so this "
                               "slice is biased toward abundant/ubiquitous features — PC structure "
                               "may reflect that, not biology.", icon="⚠️")
                pcnames = list(scores.columns)
                cc0, cc1, cc2 = st.columns(3)
                # default to a cohort/batch column when present — PC1 often tracks it
                _color_opts = ["Selected groups", "sample_type"] + meta_cols
                _cohort = next((c for c in meta_cols
                                if re.search(r"study|cohort|site|batch|plate|adrc", c, re.I)), None)
                color_by = cc0.selectbox("Color by", _color_opts,
                                         index=_color_opts.index(_cohort) if _cohort else 0,
                                         key="pca_color")
                pcx = cc1.selectbox("X axis", pcnames, index=0, key="pcx")
                pcy = cc2.selectbox("Y axis", pcnames, index=min(1, len(pcnames) - 1), key="pcy")
                highlight = st.checkbox("Outline selected A/B samples", value=True, key="pca_hl")

                dfp = scores.join(meta, how="left")
                dfp["_name"] = [id_to_name.get(s, s) for s in dfp.index]
                hov = {c: True for c in ["Condition", "Study_Name", "Sex", "batch"] if c in dfp.columns}

                if color_by == "Selected groups":
                    sa, sb = set(group_a), set(group_b)
                    dfp["_color"] = [label_a if s in sa else label_b if s in sb else "other"
                                     for s in dfp.index]
                    fig = px.scatter(
                        dfp, x=pcx, y=pcy, color="_color", hover_name="_name", hover_data=hov,
                        category_orders={"_color": [label_a, label_b, "other"]},
                        color_discrete_map={label_a: "#1f77b4", label_b: "#d62728",
                                            "other": "#d9d9d9"}, height=560,
                    )
                else:
                    s = dfp[color_by]
                    dfp["_colorby"] = s if pd.api.types.is_numeric_dtype(s) else s.astype(str)
                    fig = px.scatter(dfp, x=pcx, y=pcy, color="_colorby", hover_name="_name",
                                     hover_data=hov, height=560,
                                     labels={"_colorby": color_by})
                    fig.update_layout(legend_title_text=color_by)

                if highlight:
                    sel = [s for s in dfp.index if s in set(group_a) | set(group_b)]
                    if sel:
                        sd = dfp.loc[sel]
                        fig.add_scatter(
                            x=sd[pcx], y=sd[pcy], mode="markers", name="selected",
                            marker=dict(size=13, color="rgba(0,0,0,0)", line=dict(width=2, color="black")),
                            hoverinfo="skip",
                        )

                vx = var_ratio[pcnames.index(pcx)] * 100
                vy = var_ratio[pcnames.index(pcy)] * 100
                fig.update_xaxes(title=f"{pcx} ({vx:.1f}%)")
                fig.update_yaxes(title=f"{pcy} ({vy:.1f}%)")
                fig.update_layout(margin=dict(l=10, r=10, t=30, b=10))
                st.plotly_chart(fig, width="stretch")
                st.caption(
                    f"{len(pca_ids)} samples · {n_used:,} complete features · variance explained: "
                    + ", ".join(f"{n}={v * 100:.1f}%" for n, v in zip(pcnames, var_ratio))
                )
                st.session_state["pca_result"] = {"fig": fig,
                    "summary": f"{feat_mode} feature set · {len(pca_ids)} samples · "
                               f"{n_used:,} complete features. Variance explained: "
                               + ", ".join(f"{n}={v * 100:.1f}%"
                                           for n, v in zip(pcnames, var_ratio)) + "."}

    # ============================== EV markers tab ==============================
    with tab_ev:
        if data["level"] != "protein":
            st.warning("EV markers are matched by gene name - switch to **protein** level.")
        else:
            st.caption("Canonical extracellular-vesicle markers (MISEV-style), matched by gene "
                       "symbol. Positive markers should be high in EV preps; contamination low.")
            cats = list(EV_MARKERS)
            csel = st.multiselect("Marker categories", cats, default=cats, key="ev_cats")
            incl_neg = st.checkbox("Include contamination markers (expect low)", value=True, key="ev_neg")
            custom = st.text_input("Add custom genes (comma-separated)", key="ev_custom")

            pos_genes = [g for c in csel for g in EV_MARKERS[c]]
            neg_genes = EV_NEGATIVE if incl_neg else []
            custom_genes = [g.strip() for g in custom.split(",") if g.strip()]
            req = list(dict.fromkeys(pos_genes + neg_genes + custom_genes))
            gene_cat = {g: c for c in csel for g in EV_MARKERS[c]}
            for g in neg_genes:
                gene_cat[g] = "Contamination"
            for g in custom_genes:
                gene_cat.setdefault(g, "Custom")

            found, missing = match_genes(annot, label_col, expr_log2, req)
            if not found:
                st.warning("None of the selected markers were found in this dataset.")
            else:
                st.caption(f"Found **{len(found)}/{len(req)}** markers."
                           + (f"  Not detected: {', '.join(missing)}" if missing else ""))

                ev_types = st.multiselect("Samples to include (by type)", types, default=types, key="ev_types")
                ev_ids = meta[meta["sample_type"].isin(ev_types)]["sample_id"].tolist()
                gcols = ["sample_type"] + cat_cols
                gcol = st.selectbox("Group / order columns by", gcols,
                                    index=gcols.index("Condition") if "Condition" in gcols else 0,
                                    key="ev_gcol")
                view = st.radio("View", ["Group means", "Per sample"], horizontal=True, key="ev_view")

                if len(ev_ids) < 1:
                    st.warning("No samples selected.")
                else:
                    # marker x sample log2 matrix, rows ordered by category then gene
                    cat_order = list(EV_MARKERS) + ["Contamination", "Custom"]
                    ordered = sorted([g for g in req if g in found],
                                     key=lambda g: (cat_order.index(gene_cat.get(g, "Custom")), g))
                    mat = expr_log2.loc[[found[g] for g in ordered], ev_ids]
                    mat.index = ordered

                    gseries = meta.loc[ev_ids, gcol].astype(str)
                    if view == "Group means":
                        M = mat.T.groupby(gseries.values).mean().T
                    else:
                        order = gseries.sort_values().index.tolist()
                        M = mat[order]
                        M.columns = [id_to_name.get(s, s) for s in order]

                    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
                    hm = px.imshow(Z, color_continuous_scale="RdBu_r", zmin=-2, zmax=2, aspect="auto",
                                   labels=dict(color="z (per marker)"),
                                   height=max(360, 17 * len(ordered)))
                    hm.update_layout(margin=dict(l=10, r=10, t=40, b=10),
                                     title=f"EV markers (row z-scored log2) — {view.lower()} by {gcol}")
                    st.plotly_chart(hm, width="stretch")
                    st.session_state["ev_result"] = {"fig": hm,
                        "summary": f"{len(ordered)} canonical EV markers, row z-scored log2 "
                                   f"({view.lower()} by {gcol})."}
                    if view == "Per sample":
                        st.caption(f"Columns ordered by {gcol}.")

                    st.subheader("Marker detail")
                    pickg = st.selectbox("Marker", ordered, key="ev_pick")
                    yv = expr_log2.loc[found[pickg], ev_ids]
                    long = pd.DataFrame({
                        "log2 abundance": yv.values,
                        gcol: meta.loc[ev_ids, gcol].astype(str).values,
                        "sample": [id_to_name.get(s, s) for s in ev_ids],
                    })
                    bx = px.box(long, x=gcol, y="log2 abundance", points="all", color=gcol,
                                hover_data=["sample"], height=380)
                    bx.update_layout(showlegend=False, margin=dict(l=10, r=10, t=40, b=10),
                                     title=f"{pickg}  ({gene_cat.get(pickg, '')})")
                    st.plotly_chart(bx, width="stretch")

    # =============================== Finder tab ===============================
    with tab_find:
        st.caption("Scan many ways to split one diagnosis/condition column into two groups "
                   "and rank them by signal — to spot promising contrasts (e.g. co-pathology "
                   "groupings) without testing each by hand.")
        render_finder(expr_log2, "scan")

    # =========================== Cryptic finder tab ===========================
    with tab_cfind:
        st.caption("Same search, but limma is run on **cryptic peptides only** — so the FDR "
                   "correction is over just those features (far more powerful for low-abundance "
                   "cryptic signal than correcting over the whole proteome).")
        cterm = st.text_input("Accession contains", value="cryptic", key="cfind_term").strip()
        cf_ids, _ = detect_cryptic(cterm)
        if data["level"] != "peptide":
            st.warning("Load at **peptide** level — cryptic peptides are folded into canonical "
                       "proteins during the protein rollup.")
        elif not cf_ids:
            st.warning(f"No cryptic features matched '{cterm}'.")
        else:
            st.success(f"Restricting the scan to **{len(cf_ids)}** cryptic peptides.")
            render_finder(expr_log2.loc[cf_ids], "cscan",
                          note=f"{len(cf_ids)} cryptic peptides", unit="cryptic peptides")

    # =============================== Cryptic tab ==============================
    with tab_cryptic:
        st.caption("Cryptic / non-canonical features. Cryptic peptides are flagged at the "
                   "**peptide** level (their `Protein` field carries the term, e.g. "
                   "`CRYPTIC_UNIQUE|...`); the protein rollup folds them into canonical "
                   "proteins, so load at **peptide** level to see them.")
        term = st.text_input("Accession contains", value="cryptic", key="crypt_term").strip()

        crypt_map: dict = {}

        # (a) annotation-based match (works when the term is in a protein accession column)
        mask = pd.Series(False, index=expr_log2.index)
        if term:
            for c in annot.columns:
                col = annot[c]
                if col.dtype == object or str(col.dtype).startswith("string"):
                    hit = col.astype(str).str.contains(re.escape(term), case=False, na=False)
                    mask = mask | hit.reindex(mask.index, fill_value=False)
            idx_hit = pd.Series(expr_log2.index.astype(str), index=expr_log2.index) \
                .str.contains(re.escape(term), case=False, na=False)
            mask = mask | idx_hit
        crypt_ids = list(expr_log2.index[mask])

        # (b) peptide-level match via the merged_data Protein field
        if data["level"] == "peptide" and term:
            mp = Path(st.session_state.get("prism_key", ("",))[0] or "") / "merged_data.parquet"
            if mp.exists():
                try:
                    crypt_map = st.cache_data(cryptic_peptide_map)(str(mp), term)
                    crypt_ids = sorted(set(crypt_ids) | {p for p in expr_log2.index if p in crypt_map})
                except Exception as e:  # noqa: BLE001
                    st.warning(f"Couldn't read merged_data.parquet for cryptic flags: {e}")
            else:
                st.info("`merged_data.parquet` not found in the output_dir — needed to flag "
                        "cryptic peptides at peptide level.")

        if not term:
            st.info("Enter an accession term (e.g. `cryptic`).")
        elif not crypt_ids:
            if data["level"] == "protein":
                st.warning(f"No protein matched '{term}'. Cryptic peptides get folded into "
                           "canonical proteins during the protein rollup — **reload at peptide "
                           "level** to see them.")
            else:
                st.warning(f"No peptides matched '{term}'.")
        else:
            st.success(f"Found **{len(crypt_ids)}** cryptic features matching '{term}'.")
            labelcol = label_col if label_col in annot.columns else None
            # gene map (UniProt) — show real gene names instead of TrEMBL accessions
            _gm_path = next((p for p in [
                Path(st.session_state.get("prism_key", ("",))[0] or "") / "cryptic_gene_map.csv",
                Path(__file__).parent / "cryptic_gene_map.csv"] if p.exists()), None)
            gene_map_cr = st.cache_data(load_cryptic_gene_map)(str(_gm_path)) if _gm_path else {}

            def _crypt_label(p):
                g = gene_map_cr.get(p, {}).get("gene")
                if g:
                    return f"{g} · {p}"
                return f"{cryptic_short_label(crypt_map[p])} · {p}" if p in crypt_map else str(p)
            disp = {p: _crypt_label(p) for p in crypt_ids}
            if gene_map_cr:
                _ngene = sum(1 for p in crypt_ids if gene_map_cr.get(p, {}).get("gene"))
                st.caption(f"Gene names resolved for {_ngene}/{len(crypt_ids)} cryptic peptides "
                           "(via UniProt map).")

            # cross-reference the current contrast, if a fit exists
            analysis = st.session_state.get("analysis")
            if analysis is not None:
                resA = analysis["res"]
                la_, lb_ = analysis.get("label_a", "A"), analysis.get("label_b", "B")
                cr = resA.reindex(crypt_ids).dropna(subset=["logFC", "P.Value"])
                # Choose the multiplicity scope for cryptic calls. Genome-wide FDR corrects
                # these features against the whole proteome (very conservative); cryptic-
                # family FDR corrects within the cryptic set only (the pre-specified target).
                if sig_basis == "raw":
                    crypt_basis = "raw"
                    st.caption("Significance: **raw p-value** — exploratory, no multiple-testing "
                               "correction. Switch to FDR in the sidebar for corrected calls.")
                else:
                    _family = st.checkbox(
                        "Correct within the cryptic set only (family-wise FDR)", value=True,
                        key="crypt_family",
                        help="BH-correct across just the cryptic features (the pre-specified "
                             "hypothesis family) rather than the whole proteome — the right "
                             "multiplicity scope when cryptic peptides are the study target. "
                             "Unchecked = genome-wide FDR, which is far more conservative.")
                    crypt_basis = "family" if _family else "fdr"
                sig, p_thresh_cr = _significance(cr, crypt_basis, q_cut, lfc_cut)
                cr = cr.assign(significant=sig)
                if crypt_basis == "raw":
                    call_col = "p.value"
                    cr = cr.assign(**{call_col: cr["P.Value"]})
                elif crypt_basis == "family":
                    call_col = "cryptic.FDR"
                    cr = cr.assign(**{call_col: _bh(cr["P.Value"].to_numpy())})
                else:
                    call_col = "adj.P.Val"
                _bw = {"raw": "raw p", "family": "cryptic-FDR", "fdr": "FDR"}[crypt_basis]
                c1, c2 = st.columns(2)
                c1.metric(f"Cryptic features in the {lb_}-vs-{la_} fit", len(cr))
                c2.metric(f"...significant ({_bw} < {q_cut:g})", int(cr["significant"].sum()))
                if len(cr):
                    lbl = cr.assign(
                        feature=[disp.get(i, str(i)) for i in cr.index],
                        cryptic_protein=[crypt_map.get(i, "") for i in cr.index],
                        _nlp=-np.log10(cr["P.Value"].clip(lower=1e-300)),
                        _dir=np.where(cr["significant"],
                                      np.where(cr["logFC"] > 0, f"up in {lb_}", f"up in {la_}"), "ns"))
                    fig = px.scatter(lbl, x="logFC", y="_nlp", color="_dir", hover_name="feature",
                                     color_discrete_map={f"up in {lb_}": "#d62728",
                                                         f"up in {la_}": "#1f77b4", "ns": "#b0b0b0"},
                                     labels={"logFC": f"log2 FC ({lb_}/{la_})", "_nlp": "-log10 P"},
                                     title=f"Cryptic features — {lb_} vs {la_} ({_bw})", height=420)
                    _add_volcano_guides(fig, p_thresh_cr, lfc_cut, crypt_basis, q_cut)
                    if p_thresh_cr is None:
                        st.caption(f"No cryptic features pass {_bw} < {q_cut:g} at these thresholds.")
                    fig.update_layout(legend_title_text="", margin=dict(l=10, r=10, t=40, b=10))
                    st.plotly_chart(fig, width="stretch")
                    cols = ["feature", "logFC", call_col, "significant", "cryptic_protein"]
                    st.dataframe(lbl.sort_values("P.Value")[cols].style.format(
                        {"logFC": "{:.2f}", call_col: "{:.2e}"}), width="stretch", height=240)
                    st.session_state["cryptic_result"] = {
                        "fig": fig,
                        "summary": (f"{int(cr['significant'].sum())} of {len(cr)} cryptic features "
                                    f"significant in {lb_} vs {la_} at {_bw} < {q_cut:g}, "
                                    f"|log2FC| >= {lfc_cut:g}."),
                        "table": lbl.sort_values("P.Value")[
                            ["feature", "logFC", call_col, "significant"]].round(4).head(25)}
                else:
                    st.session_state.pop("cryptic_result", None)
            else:
                st.session_state.pop("cryptic_result", None)
                st.caption("Run a contrast in the **Differential** tab to see which cryptic "
                           "features are significant.")

            # ---- detection (present/absent) test — the right analysis for cryptic peptides ----
            st.subheader("Detection (present / absent)")
            st.caption("Cryptic peptides are on/off — the dense Area matrix integrates a borrowed "
                       "baseline where they're truly absent, so **abundance is misleading**. This "
                       "tests genuine *detection* (DetectionQValue < 0.01) per peptide between "
                       "groups A and B.")
            mp_det = Path(st.session_state.get("prism_key", ("",))[0] or "") / "merged_data.parquet"
            if data["level"] != "peptide" or not mp_det.exists():
                st.info("Detection needs peptide-level data + `merged_data.parquet`.")
            elif not (group_a and group_b):
                st.info("Define **Group A** and **Group B** in the sidebar to test detection "
                        "(e.g. set PrimaryDx → A=CO, B=FTD).")
            else:
                dtest = st.radio(
                    "Test", ["Firth logistic (covariate-adjusted)", "Fisher exact (unadjusted)"],
                    horizontal=True, key="crypt_dettest",
                    help="Firth logistic adjusts for the sidebar covariates (e.g. cohort, age) and "
                         "handles separation at small n; Fisher is the simpler unadjusted 2x2.")
                detmat = st.cache_data(cryptic_detection_matrix)(str(mp_det), term)
                detmat = detmat.reindex([p for p in crypt_ids if p in detmat.index])
                dinfo, dres = {}, pd.DataFrame()
                try:
                    # cache so the test doesn't re-run on every Streamlit interaction
                    if dtest.startswith("Firth"):
                        cov_for_det = meta[chosen_cov] if chosen_cov else None
                        dres, dinfo = st.cache_data(detection_test_glm)(
                            detmat, group_a, group_b, cov_for_det)
                    else:
                        dres = st.cache_data(detection_test)(detmat, group_a, group_b)
                except Exception as e:  # noqa: BLE001
                    st.warning(f"Detection test failed: {e}")

                if dinfo.get("warning"):
                    st.error(f"⚠️ {dinfo['warning']}  "
                             f"(group is {dinfo.get('group_collinearity_r2', 0):.0%} explained by "
                             "the covariates).")
                if dinfo.get("covariates_used"):
                    st.caption("Adjusted for: " + ", ".join(dinfo["covariates_used"])
                               + (f"  ·  dropped: {', '.join(dinfo['dropped'])}" if dinfo.get("dropped") else ""))
                if dinfo.get("n_nonconverged"):
                    st.caption(f"⚠️ {dinfo['n_nonconverged']} peptide fit(s) did not fully converge "
                               "(estimates retained; interpret those rows with care).")
                if len(dres):
                    dres = dres.assign(feature_label=[disp.get(f, f) for f in dres["feature"]])
                    m1, m2 = st.columns(2)
                    m1.metric("Differentially detected (q<0.05)", int((dres["q"] < 0.05).sum()))
                    m2.metric("...at raw p<0.05 (candidates)", int((dres["p"] < 0.05).sum()))
                    aset = [s for s in group_a if s in detmat.columns]
                    bset = [s for s in group_b if s in detmat.columns]
                    st.caption(f"Detection rate — {label_a}: {detmat[aset].mean().mean():.0%}  ·  "
                               f"{label_b}: {detmat[bset].mean().mean():.0%} (across cryptic peptides).")
                    eff = ["logOR"] if "logOR" in dres.columns else []
                    cols_d = (["feature_label", "det_A", "n_A", "rate_A", "det_B", "n_B", "rate_B"]
                              + eff + ["p", "q"])
                    st.dataframe(
                        dres[cols_d].rename(columns={"rate_A": f"rate[{label_a}]",
                                                     "rate_B": f"rate[{label_b}]"})
                        .style.format({f"rate[{label_a}]": "{:.0%}", f"rate[{label_b}]": "{:.0%}",
                                       "logOR": "{:.2f}", "p": "{:.2e}", "q": "{:.2e}"}),
                        width="stretch", height=300)
                    st.download_button("Download detection table (CSV)", dres.to_csv(index=False).encode(),
                                       file_name="cryptic_detection.csv", key="crypt_det_dl")

                    # per-cohort detection rates — pooled rates hide cohort confounding
                    cohort_choices = [c for c in cat_cols if c in meta.columns]
                    if cohort_choices:
                        strat = st.selectbox("Show detection rate by (check for confounding)",
                                             ["(none)"] + cohort_choices, key="crypt_strat")
                        if strat != "(none)":
                            sub = meta.loc[aset + bset].copy()
                            sub["group"] = [label_a] * len(aset) + [label_b] * len(bset)
                            sub["_rate"] = [detmat[s].mean() for s in aset + bset]
                            piv = sub.pivot_table(index=strat, columns="group", values="_rate",
                                                  aggfunc="mean")
                            cnt = pd.crosstab(sub[strat].astype(str), sub["group"])
                            st.write("Mean cryptic detection rate by " + strat + " (n in parentheses):")
                            st.dataframe(piv.style.format("{:.0%}"), width="stretch")
                            st.caption("Counts per cell: " + cnt.to_dict().__str__())

            # abundance heatmap of cryptic features across samples
            st.subheader("Abundance")
            ctypes = st.multiselect("Samples to include (by type)", types, default=types, key="crypt_types")
            cids = meta[meta["sample_type"].isin(ctypes)]["sample_id"].tolist()
            gcols = ["sample_type"] + cat_cols
            cgcol = st.selectbox("Group / order columns by", gcols,
                                 index=gcols.index("Condition") if "Condition" in gcols else 0,
                                 key="crypt_gcol")
            cview = st.radio("View", ["Group means", "Per sample"], horizontal=True, key="crypt_view")

            # cap rows for readability; keep the most variable cryptic features
            mat = expr_log2.loc[crypt_ids, cids]
            if len(mat) > 50:
                mat = mat.loc[mat.var(axis=1).nlargest(50).index]
                st.caption(f"Showing the 50 most variable of {len(crypt_ids)} cryptic features.")
            mat = mat.rename(index={p: disp.get(p, str(p)) for p in mat.index})

            gseries = meta.loc[cids, cgcol].astype(str)
            if cview == "Group means":
                M = mat.T.groupby(gseries.values).mean().T
            else:
                order = gseries.sort_values().index.tolist()
                M = mat[order]
                M.columns = [id_to_name.get(s, s) for s in order]
            Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
            hm = px.imshow(Z, color_continuous_scale="RdBu_r", zmin=-2, zmax=2, aspect="auto",
                           labels=dict(color="z (per feature)"), height=max(320, 16 * len(M)))
            hm.update_layout(margin=dict(l=10, r=10, t=40, b=10),
                             title=f"Cryptic features (row z-scored log2) — {cview.lower()} by {cgcol}")
            st.plotly_chart(hm, width="stretch")

            # per-feature box / strip plot of the cryptic peptide across groups A vs B
            st.subheader("Feature detail")
            # restrict the picker to cryptic peptides significant in the current contrast
            _an = st.session_state.get("analysis")
            sig_crypt = []
            if _an is not None:
                _r = _an["res"]
                _sig = set(_r.index[(_r["adj.P.Val"] < q_cut) & (_r["logFC"].abs() >= lfc_cut)])
                sig_crypt = [p for p in crypt_ids if p in _sig]
            only_sig_c = st.checkbox("Significant peptides only", value=True, key="crypt_pick_sig")
            if only_sig_c and sig_crypt:
                pick_pool = sig_crypt
            else:
                pick_pool = crypt_ids
                if only_sig_c:
                    st.caption("No cryptic peptides significant at the current FDR/|log2FC| "
                               "thresholds (run a contrast in **Differential** first, or try the "
                               "Detection test above) — showing all.")
            pickc = st.selectbox(f"Cryptic feature ({len(pick_pool)})", pick_pool,
                                 format_func=lambda p: disp.get(p, str(p)), key="crypt_pick")
            if not (group_a and group_b):
                st.info("Define **Group A** and **Group B** in the sidebar to see the "
                        "A-vs-B box plot for this cryptic peptide.")
            else:
                vals = expr_log2.loc[pickc, group_a + group_b]
                long = pd.DataFrame({
                    "log2 abundance": vals.values,
                    "group": [label_a] * len(group_a) + [label_b] * len(group_b),
                    "sample": [id_to_name.get(s, s) for s in (group_a + group_b)],
                })
                bx = px.box(long, x="group", y="log2 abundance", points="all", color="group",
                            category_orders={"group": [label_a, label_b]}, hover_data=["sample"],
                            color_discrete_map={label_a: "#1f77b4", label_b: "#d62728"}, height=400)
                ttl = disp.get(pickc, str(pickc))
                if pickc in crypt_map:
                    ttl += f"  ({crypt_map[pickc].split('|source:')[0]})"
                bx.update_layout(showlegend=False, margin=dict(l=10, r=10, t=40, b=10), title=ttl)
                st.plotly_chart(bx, width="stretch")

            # ---- cryptic vs canonical (stoichiometry) ----
            st.subheader("Cryptic vs canonical (stoichiometry)")
            out_dir_cr = st.session_state.get("prism_key", ("",))[0]
            if not gene_map_cr:
                st.info("Needs `cryptic_gene_map.csv` (UniProt accession→gene map) in the "
                        "output_dir or app folder to link cryptic peptides to their genes.")
            elif not (group_a and group_b):
                st.info("Define Group A and B in the sidebar to compare cryptic vs canonical.")
            else:
                # cryptic peptides that have a resolved gene AND a canonical protein in the data
                prot = st.cache_data(load_prism)(out_dir_cr, "protein")
                pannot = prot["annot"]
                gene_to_pg = {}
                for pg, g in pannot["leading_gene_name"].astype(str).items():
                    for gg in g.split("/"):
                        gene_to_pg.setdefault(gg.strip().upper(), pg)
                linkable = [p for p in crypt_ids
                            if gene_map_cr.get(p, {}).get("gene", "").upper() in gene_to_pg]
                # optional: restrict to cryptic peptides significant in the current contrast
                _ans = st.session_state.get("analysis")
                sig_link = []
                if _ans is not None:
                    _rs = _ans["res"]
                    _sigs = set(_rs.index[(_rs["adj.P.Val"] < q_cut) & (_rs["logFC"].abs() >= lfc_cut)])
                    sig_link = [p for p in linkable if p in _sigs]
                only_sig_s = st.checkbox("Significant cryptic peptides only", value=False,
                                         key="crypt_stoich_sig")
                pool_s = sig_link if (only_sig_s and sig_link) else linkable
                if only_sig_s and not sig_link:
                    st.caption("No significant cryptic peptides (with a canonical match) at current "
                               "thresholds — showing all.")
                if not pool_s:
                    st.caption("None of these cryptic peptides have a canonical protein quantified "
                               "in the protein-level output.")
                else:
                    from scipy.stats import ttest_ind
                    from statsmodels.stats.multitest import multipletests

                    cols = group_a + group_b

                    def _ttest_p(v):
                        a, b = v[group_a].dropna(), v[group_b].dropna()
                        if len(a) < 2 or len(b) < 2:
                            return float("nan")
                        return float(ttest_ind(a, b, equal_var=False).pvalue)

                    def _star(q):
                        if q != q:
                            return "n/a"
                        return ("***" if q < 0.001 else "**" if q < 0.01
                                else "*" if q < 0.05 else "ns")

                    def _bh_map(pmap):
                        ks = [k for k, val in pmap.items() if val == val]
                        if not ks:
                            return {}
                        return dict(zip(ks, multipletests([pmap[k] for k in ks], method="fdr_bh")[1]))

                    # Test A-vs-B across the WHOLE browsable family (all linkable cryptic
                    # peptides + their unique canonical proteins), then BH-correct — so the
                    # marker reflects the multiplicity of scrolling through the set. ~350 cheap
                    # t-tests; computed inline (caching big matrices would cost more than this).
                    crypt_p_all = {p: _ttest_p(expr_log2.loc[p, cols].astype(float)) for p in linkable}
                    canon_pgs = {gene_to_pg[gene_map_cr[p]["gene"].upper()] for p in linkable}
                    canon_p_all = {pg2: _ttest_p(prot["expr_log2"].loc[pg2, cols].astype(float))
                                   for pg2 in canon_pgs}
                    crypt_q_all, canon_q_all = _bh_map(crypt_p_all), _bh_map(canon_p_all)

                    picks = st.selectbox(f"Cryptic peptide ({len(pool_s)} shown)",
                                         pool_s, format_func=lambda p: disp.get(p, p), key="crypt_stoich")
                    gene = gene_map_cr[picks]["gene"]
                    pg = gene_to_pg[gene.upper()]
                    grp = [label_a] * len(group_a) + [label_b] * len(group_b)
                    crypt_v = expr_log2.loc[picks, cols].astype(float)
                    canon_v = prot["expr_log2"].loc[pg, cols].astype(float)
                    cp, cq = crypt_p_all.get(picks, float("nan")), crypt_q_all.get(picks, float("nan"))
                    kp, kq = canon_p_all.get(pg, float("nan")), canon_q_all.get(pg, float("nan"))

                    long = pd.DataFrame({
                        "log2 abundance": list(crypt_v.values) + list(canon_v.values),
                        "group": grp + grp,
                        "track": ["cryptic peptide"] * len(cols) + [f"canonical {gene}"] * len(cols),
                    })
                    fig = px.box(long, x="track", y="log2 abundance", color="group", points="all",
                                 category_orders={"group": [label_a, label_b],
                                                  "track": ["cryptic peptide", f"canonical {gene}"]},
                                 color_discrete_map={label_a: "#1f77b4", label_b: "#d62728"}, height=440)
                    # A-vs-B marker per track: star reflects the BH q-value (corrected over the set)
                    fig.add_annotation(x="cryptic peptide", y=float(crypt_v.max()), yshift=22,
                                       showarrow=False, font=dict(size=13),
                                       text=(f"{_star(cq)}  p={cp:.2g}, q={cq:.2g}" if cp == cp else "n/a"))
                    fig.add_annotation(x=f"canonical {gene}", y=float(canon_v.max()), yshift=22,
                                       showarrow=False, font=dict(size=13),
                                       text=(f"{_star(kq)}  p={kp:.2g}, q={kq:.2g}" if kp == kp else "n/a"))
                    fig.update_layout(margin=dict(l=10, r=10, t=50, b=10),
                                      title=f"{gene}: cryptic peptide vs canonical protein  "
                                            f"({label_b} vs {label_a}, Welch t-test)",
                                      boxmode="group")
                    st.plotly_chart(fig, width="stretch")
                    st.caption(f"Stars = **BH-FDR q** (corrected across the {len(linkable)} cryptic "
                               f"peptides / {len(set(canon_p_all))} canonical proteins in this view); "
                               "p is the raw Welch t-test. Unadjusted for cohort — confirm hits in the "
                               "Detection/Differential tabs. Canonical q here is over the cryptic-linked "
                               "proteins only, not the whole proteome.")
                    # relative cryptic level = cryptic peptide − canonical protein (log2), by group
                    ratio = pd.DataFrame({"relative cryptic (log2 cryptic − canonical)":
                                          (crypt_v.values - canon_v.values), "group": grp})
                    rb = px.box(ratio, x="group", y="relative cryptic (log2 cryptic − canonical)",
                                points="all", color="group",
                                category_orders={"group": [label_a, label_b]},
                                color_discrete_map={label_a: "#1f77b4", label_b: "#d62728"}, height=340)
                    rb.update_layout(showlegend=False, margin=dict(l=10, r=10, t=40, b=10),
                                     title=f"{gene}: cryptic-to-canonical ratio (higher = more cryptic-spliced)")
                    st.plotly_chart(rb, width="stretch")
                    cs = gene_map_cr[picks]["cryptic_string"]
                    region = next((x for x in cs.split("|") if x.startswith("tryptic:")), "")
                    st.caption(f"Cryptic peptide `{picks}` maps to **{gene}** "
                               f"({gene_map_cr[picks]['accession']}), canonical protein `{pg}` "
                               f"({pannot.loc[pg, 'n_peptides']} peptides). Cryptic-exon region: {region}. "
                               "⚠️ Different scales (one peptide vs protein rollup) — read the "
                               "**relative change between groups**, not the absolute ratio.")

            dl = annot.loc[crypt_ids].join(expr_log2.loc[crypt_ids])
            if crypt_map:
                dl.insert(0, "cryptic_protein", [crypt_map.get(i, "") for i in dl.index])
            st.download_button("Download cryptic feature table (CSV)", dl.to_csv().encode(),
                               file_name="cryptic_features.csv")

    # =========================== Longitudinal tab ============================
    with tab_long:
        st.caption("Within-subject change over time. Pick the **subject** and a numeric "
                   "**timepoint** column; fits `~ subject + time` per feature (each subject is "
                   "its own baseline) and reports the time slope. Exploratory — small cohorts "
                   "rarely clear FDR, so read the slopes and trajectories, not the verdict.")
        # a timepoint column needs >=2 distinct numeric values (tolerates missing
        # timepoints, which would sink a fixed populated-fraction threshold).
        time_choices = [c for c in meta_cols
                        if pd.to_numeric(meta[c], errors="coerce").nunique() >= 2]
        if not meta_cols or not time_choices:
            st.info("Need a subject column and a numeric timepoint column in the metadata "
                    "(load a clinical CSV that carries them).")
        else:
            lc1, lc2 = st.columns(2)
            _sdef = next((c for c in meta_cols
                          if re.search(r"subj|nvt|entity|patient|donor|_id$|^id$", c, re.I)),
                         meta_cols[0])
            subj_col = lc1.selectbox("Subject column", meta_cols,
                                     index=meta_cols.index(_sdef), key="long_subj")
            _tdef = next((c for c in time_choices
                          if re.search(r"time|month|day|week|visit|seq|age", c, re.I)),
                         time_choices[0])
            time_col = lc2.selectbox("Timepoint column (numeric)", time_choices,
                                     index=time_choices.index(_tdef), key="long_time")
            grp_col = st.selectbox("Colour trajectories by (optional)", ["(subject)"] + cat_cols,
                                   key="long_grp")
            crypt_only, long_ids = False, None
            if data["level"] == "peptide":
                crypt_only = st.checkbox("Restrict to cryptic peptides", value=False,
                                         key="long_crypt")
                if crypt_only:
                    lterm = st.text_input("Accession contains", value="cryptic",
                                          key="long_cterm").strip()
                    long_ids, _ = detect_cryptic(lterm)
                    st.caption(f"{len(long_ids)} cryptic peptides."
                               if long_ids else f"No cryptic peptides matched '{lterm}'.")

            if st.button("Run longitudinal scan", type="primary", key="long_run"):
                mat = expr_log2.loc[long_ids] if (crypt_only and long_ids) else expr_log2
                try:
                    lres, linfo = longitudinal_time_slope(
                        mat, meta.reindex(expr_log2.columns)[subj_col],
                        meta.reindex(expr_log2.columns)[time_col])
                    st.session_state["long_result"] = {
                        "res": lres, "info": linfo, "time_col": time_col,
                        "subj_col": subj_col, "grp_col": grp_col}
                except Exception as e:  # noqa: BLE001
                    st.error(f"Longitudinal scan failed: {e}")
                    st.session_state.pop("long_result", None)

            lr = st.session_state.get("long_result")
            if lr is not None:
                info, lres, tcol = lr["info"], lr["res"], lr["time_col"]
                st.write(f"**{info['n_features']:,}** features · **{info['n_subjects']}** subjects "
                         f"(≥2 timepoints) · **{info['n_samples']}** samples · residual df "
                         f"**{info['residual_df']}**.")
                nsig = int((lres["adj.P.Val"] < q_cut).sum())
                if nsig == 0:
                    st.warning("No feature clears FDR — expected at this sample size. Rank by "
                               "slope / raw P to prioritise and eyeball the trajectories; don't "
                               "read significance into it.", icon="⚠️")
                else:
                    st.success(f"{nsig} feature(s) at adj.P<{q_cut:g} — treat as leads, not proof.")
                disp = lres.rename(columns={"slope": f"slope_per_{tcol}"})
                st.dataframe(disp.head(300).style.format(
                    {f"slope_per_{tcol}": "{:.3f}", "t": "{:.2f}", "P.Value": "{:.2e}",
                     "adj.P.Val": "{:.2e}", "AveExpr": "{:.1f}"}), width="stretch", height=280)
                st.download_button("Download slope table (CSV)", lres.to_csv().encode(),
                                   file_name="longitudinal_slopes.csv", key="long_dl")

                st.subheader("Per-subject trajectory")
                feat = st.selectbox("Feature", lres.index.tolist(), key="long_feat")
                tser = pd.to_numeric(meta.reindex(expr_log2.columns)[lr["time_col"]],
                                     errors="coerce")
                sser = meta.reindex(expr_log2.columns)[lr["subj_col"]].astype(str)
                pcols = [c for c in expr_log2.columns
                         if pd.notna(tser.get(c)) and sser.get(c, "nan").lower() not in
                         ("nan", "na", "none", "")]
                tdf = pd.DataFrame({
                    "time": tser[pcols].to_numpy(float),
                    "log2 abundance": expr_log2.loc[feat, pcols].to_numpy(float),
                    "subject": sser[pcols].to_numpy()})
                if lr["grp_col"] != "(subject)":
                    tdf["group"] = meta.reindex(pcols)[lr["grp_col"]].astype(str).to_numpy()
                    color = "group"
                else:
                    color = "subject"
                tdf = tdf.dropna(subset=["time", "log2 abundance"]).sort_values("time")
                if len(tdf):
                    fig = px.line(tdf, x="time", y="log2 abundance", color=color,
                                  line_group="subject", markers=True, hover_name="subject",
                                  labels={"time": tcol}, height=430,
                                  title=f"{feat} — per-subject trajectory over {tcol}")
                    fig.update_layout(margin=dict(l=10, r=10, t=40, b=10))
                    st.plotly_chart(fig, width="stretch")

    # =========================== Heterogeneity tab ===========================
    with tab_het:
        st.caption("Disease heterogeneity. Pooling collapses a heterogeneous cohort to one "
                   "average, so a feature carried by a subset of patients gets diluted toward "
                   "baseline and 'disappears'. This uses **genuine detection** (merged_data "
                   "DetectionQValue), not the dense borrowed-baseline matrix, to show how many "
                   "individuals actually carry each feature.")
        mp_h = Path(st.session_state.get("prism_key", ("",))[0] or "") / "merged_data.parquet"
        if data["level"] != "peptide" or not mp_h.exists():
            st.info("Heterogeneity needs **peptide** level + `merged_data.parquet` (real detection).")
        elif not cat_cols:
            st.info("Need a categorical group column in the metadata.")
        else:
            hc1, hc2 = st.columns(2)
            _gdef = next((c for c in cat_cols if re.search(r"dx|diagn|group|condition|key", c, re.I)),
                         cat_cols[0])
            grp_h = hc1.selectbox("Group / diagnosis column", cat_cols,
                                  index=cat_cols.index(_gdef), key="het_grp")
            qthr = hc2.slider("Detection q-value", 0.001, 0.05, 0.01, 0.001, key="het_q")
            gvals = meta.reindex(expr_log2.columns)[grp_h].astype(str)
            uniq = sorted(v for v in gvals.dropna().unique() if v.lower() != "nan")
            _pguess = [v for v in uniq if re.search(r"pool|ref|control|healthy", v, re.I)]
            pool_vals = st.multiselect("Values that are POOLED / reference (excluded from "
                                       "individuals)", uniq, default=_pguess, key="het_pools")
            indiv_vals = [v for v in uniq if v not in pool_vals]
            focus = st.selectbox("Focus group (its individuals are analysed)", indiv_vals or uniq,
                                 key="het_focus")

            if st.button("Compute heterogeneity", type="primary", key="het_run"):
                try:
                    det = st.cache_data(cryptic_detection_matrix)(str(mp_h), None, float(qthr))
                    indiv = [s for s in expr_log2.columns if str(gvals.get(s)) == focus]
                    pools = [s for s in expr_log2.columns
                             if str(gvals.get(s)) in set(pool_vals)]
                    htab = heterogeneity_table(expr_log2, det, indiv)
                    st.session_state["het_result"] = {"tab": htab, "det": det, "indiv": indiv,
                                                      "pools": pools, "focus": focus}
                except Exception as e:  # noqa: BLE001
                    st.error(f"Heterogeneity failed: {e}")
                    st.session_state.pop("het_result", None)

            hr = st.session_state.get("het_result")
            if hr:
                htab, det, indiv = hr["tab"], hr["det"], hr["indiv"]
                pools, focus = hr["pools"], hr["focus"]
                seen = htab[htab["n_detected"] > 0]
                subset = int((seen["prevalence"] < 0.5).sum())
                st.write(f"**{focus}**: {len(indiv)} individuals · {len(pools)} pool/ref samples · "
                         f"**{len(seen):,}** features detected in ≥1 individual.")
                if len(indiv) >= 3 and len(seen):
                    st.warning(f"**{subset:,} features are detected in <50% of {focus} individuals** "
                               "— subset markers a pooled sample dilutes toward baseline. That is the "
                               "heterogeneity, not absence.", icon="⚠️")
                if len(seen):
                    b = pd.cut(seen["prevalence"], [0, .25, .5, .75, .99, 1.001], right=False,
                               labels=["1-25%", "25-50%", "50-75%", "75-99%", "100%"])
                    vc = b.value_counts().sort_index()
                    figh = px.bar(x=vc.index.astype(str), y=vc.to_numpy(), height=300,
                                  labels={"x": f"detected in this fraction of {focus} individuals",
                                          "y": "features"},
                                  title="Detection prevalence across individuals")
                    figh.update_layout(margin=dict(l=10, r=10, t=40, b=10))
                    st.plotly_chart(figh, width="stretch")

                # ---- cohort composition (#4): show the group is not monolithic ----
                with st.expander(f"Cohort composition — how heterogeneous are the {focus} "
                                 f"individuals (n={len(indiv)})?"):
                    imeta = meta.reindex(indiv)
                    shown_any = False
                    for c in meta_cols:
                        if c == grp_h:
                            continue
                        col = imeta[c].dropna()
                        num = pd.to_numeric(col, errors="coerce")
                        if not col.empty and num.notna().mean() > 0.7 and num.nunique() > 2:
                            st.caption(f"**{c}**: {num.min():.0f}–{num.max():.0f} "
                                       f"(median {num.median():.0f})")
                            shown_any = True
                        elif 2 <= col.nunique() <= 10:
                            vc = col.astype(str).value_counts().to_dict()
                            st.caption(f"**{c}**: " + ", ".join(f"{k} ×{v}" for k, v in vc.items()))
                            shown_any = True
                    if not shown_any:
                        st.caption("No additional varying metadata to summarise.")

                st.subheader("Per-feature heterogeneity (rare → common)")
                n_never = int((htab["n_detected"] == 0).sum())
                st.caption(f"Detected in ≥1 {focus} individual, rarest first — the subset markers "
                           f"pooling underrepresents. ({n_never:,} features detected in no {focus} "
                           "individual are hidden here but kept in the CSV.)")
                ranked = (htab[htab["n_detected"] > 0]
                          .sort_values(["prevalence", "sd_when_detected"], ascending=[True, False]))
                st.dataframe(ranked.head(300).style.format(
                    {"prevalence": "{:.0%}", "sd_when_detected": "{:.2f}",
                     "mean_log2_detected": "{:.1f}", "pooling_dilution": "{:.1f}×"}),
                    width="stretch", height=250)
                st.download_button("Download heterogeneity table (CSV)", htab.to_csv().encode(),
                                   file_name="heterogeneity.csv", key="het_dl")

                st.subheader("Look up a protein / peptide")
                lk1, lk2 = st.columns([2, 1])
                term = lk1.text_input("Accession, gene, or peptide contains", value="",
                                      key="het_lookup").strip()
                _strat_opts = [c for c in cat_cols if 2 <= meta[c].nunique() <= 8]
                strat_col = lk2.selectbox("Stratify individuals by (optional)",
                                          ["(none)"] + _strat_opts, key="het_strat")
                if term:
                    lo = term.lower()
                    hitset = {p for p in htab.index if lo in str(p).lower()}    # peptide sequence
                    for c in annot.columns:                                     # accession / protein / gene
                        col = annot[c]
                        if col.dtype == object or str(col.dtype).startswith("string"):
                            m = col.astype(str).str.contains(re.escape(term), case=False, na=False)
                            hitset |= set(col.index[m])
                    hits = (htab.loc[[p for p in htab.index if p in hitset]]
                            .sort_values("prevalence", ascending=False).index.tolist())
                    if not hits:
                        st.caption(f"No feature matched '{term}' (peptide sequence or annotation).")
                    else:
                        pick = st.selectbox(f"{len(hits)} match(es) — most-detected first",
                                            hits[:300], key="het_pick")
                        row = htab.loc[pick]
                        mc1, mc2, mc3 = st.columns(3)
                        mc1.metric("Detected in", f"{int(row['n_detected'])}/{int(row['n_individuals'])} "
                                                  f"({row['prevalence']:.0%})")
                        mc2.metric("Est. pooling dilution",
                                   f"{row['pooling_dilution']:.1f}×"
                                   if pd.notna(row["pooling_dilution"]) else "—")
                        mc3.metric("SD when detected", f"{row['sd_when_detected']:.2f}"
                                   if pd.notna(row["sd_when_detected"]) else "—")
                        strat = (None if strat_col == "(none)"
                                 else meta.reindex(expr_log2.columns)[strat_col].astype(str))
                        rows = []
                        for base, ss in (("individual", indiv), ("pool/ref", pools)):
                            for s in ss:
                                if s not in expr_log2.columns:
                                    continue
                                d = (int(det.loc[pick, s]) if (pick in det.index and s in det.columns)
                                     else 0)
                                cls = base
                                if base == "individual" and strat is not None:
                                    sv = strat.get(s)
                                    cls = str(sv) if (sv and str(sv).lower() != "nan") else "individual: NA"
                                rows.append({"log2 abundance": float(expr_log2.loc[pick, s]),
                                             "class": cls,
                                             "detection": "detected" if d else "baseline (not detected)",
                                             "sample": id_to_name.get(s, s)})
                        bdf = pd.DataFrame(rows).dropna(subset=["log2 abundance"])
                        if len(bdf):
                            figb = px.strip(bdf, x="class", y="log2 abundance", color="detection",
                                            hover_name="sample", stripmode="overlay", height=420,
                                            color_discrete_map={"detected": "#1f77b4",
                                                                "baseline (not detected)": "#c8c8c8"},
                                            title=f"{pick} — individuals vs pools"
                                                  + (f", split by {strat_col}"
                                                     if strat_col != "(none)" else ""))
                            figb.update_layout(margin=dict(l=10, r=10, t=40, b=10))
                            st.plotly_chart(figb, width="stretch")
                            st.caption("Grey = integrated baseline, not a genuine peak. A pool sits near "
                                       "the cohort average, so a feature carried by a few individuals "
                                       "reads as low/absent in it — heterogeneity, not absence.")

    # ============================ Enrichment tab =============================
    with tab_enrich:
        st.subheader("Functional enrichment & disease novelty")
        A_enr = st.session_state.get("analysis")
        if A_enr is None:
            st.info("Run **limma** on the Differential tab first, then come back.")
        elif data.get("level") != "protein" or label_col != _PROTEIN_LABEL:
            st.info("Enrichment needs gene symbols, available only at **protein** level "
                    "(annotation column `leading_gene_name`). Reload at protein level.")
        else:
            eres = A_enr["res"]
            la_e = A_enr.get("label_a", "A")
            lb_e = A_enr.get("label_b", "B")
            direction = st.radio(
                "Which hits to analyse", ["both", f"up in {lb_e}", f"up in {la_e}"],
                horizontal=True, key="enr_dir",
                help="Enriching up- and down-regulated genes separately is usually "
                     "more interpretable than pooling both directions.")
            dcode = {"both": "both", f"up in {lb_e}": "up", f"up in {la_e}": "down"}[direction]
            sig_g, bg_g = sig_and_background_genes(
                eres, annot, label_col, q_cut, lfc_cut, dcode)
            m1, m2 = st.columns(2)
            m1.metric("Significant genes", f"{len(sig_g):,}")
            m2.metric("Background (tested genes)", f"{len(bg_g):,}")
            if not sig_g:
                st.warning("No significant genes at the current thresholds — "
                           "adjust the FDR / log2FC cutoffs on the left.")
            else:
                # -------------------- GO / pathway over-representation --------------
                with st.expander("GO / pathway enrichment (g:Profiler)", expanded=True):
                    srcs = st.multiselect(
                        "Annotation sources",
                        ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"],
                        default=["GO:BP", "REAC"], key="enr_src")
                    st.caption("Sends **gene symbols only** to g:Profiler (no intensities or "
                               "metadata). Your tested genes are the statistical background, so "
                               "terms reflect enrichment above the detectable EV proteome — not "
                               "just what is abundant.")
                    if st.button("Run enrichment", key="enr_go_btn", type="primary"):
                        try:
                            with st.spinner("Querying g:Profiler…"):
                                _enr = st.cache_data(show_spinner=False)(gprofiler_enrichment)
                                st.session_state["enr_go"] = _enr(
                                    tuple(sig_g), tuple(bg_g), tuple(srcs))
                        except Exception as exc:  # noqa: BLE001 - surface any network/API error
                            st.session_state["enr_go"] = None
                            st.error(f"g:Profiler request failed: {exc}")
                    go = st.session_state.get("enr_go")
                    if go is not None and not go.empty:
                        top = go.head(25).copy()
                        top["-log10 p"] = -np.log10(top["p_value"].clip(lower=1e-300))
                        fig = px.bar(
                            top.sort_values("-log10 p"), x="-log10 p", y="term_name",
                            color="source", orientation="h",
                            height=min(720, 70 + 22 * len(top)),
                            hover_data={"term_id": True, "p_value": ":.2e",
                                        "intersection_size": True, "term_size": True})
                        fig.update_layout(yaxis_title="", legend_title="")
                        st.plotly_chart(fig, width="stretch")
                        cols = [c for c in ["source", "term_id", "term_name", "p_value",
                                            "fold_enrichment", "intersection_size", "term_size"]
                                if c in go.columns]
                        st.dataframe(
                            go[cols].style.format(
                                {"p_value": "{:.2e}", "fold_enrichment": "{:.2f}"}),
                            width="stretch", height=320)
                        st.download_button(
                            "Download enrichment (CSV)", go.to_csv(index=False),
                            file_name="enrichment_gprofiler.csv", mime="text/csv")
                    elif go is not None:
                        st.info("No terms passed g:Profiler's significance threshold.")

                # -------------------- Disease novelty (Open Targets) ----------------
                with st.expander("Known vs candidate-novel for a disease (Open Targets)"):
                    st.caption("Flags each significant gene by its Open Targets association "
                               "score to a disease you name (top ~3000 associations pulled). "
                               "**'Candidate-novel' means 'not currently annotated in Open "
                               "Targets' — a reflection of database coverage, not proof of "
                               "biological novelty.**")
                    disq = st.text_input("Disease name", key="enr_dis",
                                         placeholder="e.g. Parkinson disease")
                    if st.button("Search disease", key="enr_dis_btn") and disq.strip():
                        try:
                            with st.spinner("Searching Open Targets…"):
                                st.session_state["enr_efo"] = st.cache_data(
                                    show_spinner=False)(opentargets_search_disease)(disq)
                        except Exception as exc:  # noqa: BLE001
                            st.session_state["enr_efo"] = []
                            st.error(f"Open Targets search failed: {exc}")
                    efos = st.session_state.get("enr_efo") or []
                    if efos:
                        pick = st.selectbox(
                            "Matched disease", efos,
                            format_func=lambda t: f"{t[1]}  ({t[0]})", key="enr_efo_pick")
                        thr = st.slider(
                            "Known-association score cutoff", 0.0, 1.0, 0.1, 0.05,
                            key="enr_thr",
                            help="Open Targets overall association score; genes at or above "
                                 "this count as 'known' for the disease.")
                        if st.button("Classify hits", key="enr_cls_btn", type="primary"):
                            try:
                                with st.spinner("Fetching associations…"):
                                    scores = st.cache_data(show_spinner=False)(
                                        opentargets_disease_targets)(pick[0])
                                st.session_state["enr_nov"] = classify_novelty(sig_g, scores, thr)
                            except Exception as exc:  # noqa: BLE001
                                st.session_state["enr_nov"] = None
                                st.error(f"Open Targets association fetch failed: {exc}")
                        nov = st.session_state.get("enr_nov")
                        if nov is not None and not nov.empty:
                            n_known = int((nov["status"] == "known").sum())
                            n_novel = int((nov["status"] == "candidate-novel").sum())
                            k1, k2 = st.columns(2)
                            k1.metric("Known for disease", f"{n_known:,}")
                            k2.metric("Candidate-novel", f"{n_novel:,}")
                            st.dataframe(
                                nov.style.format({"assoc_score": "{:.3f}"}),
                                width="stretch", height=360)
                            st.download_button(
                                "Download classification (CSV)", nov.to_csv(index=False),
                                file_name="novelty_opentargets.csv", mime="text/csv")

    # ============================== Report tab ===============================
    with tab_report:
        ss = st.session_state

        # ----- Share this analysis as a live link (view-only tunnel) -----------
        if viewer is None:
            with st.expander("🚀 Share this analysis as a live link (view-only)",
                             expanded=("share_url" in ss)):
                st.caption("Spin up a password-protected, read-only copy of this app at a public "
                           "URL. Your collaborator lands on the comparison you have loaded and can "
                           "explore it (volcano, box plots, PCA) but can't change the data or see "
                           "your files. Your PC must stay on and this app must keep running.")
                # Recover control of a share still running from a previous app session
                # (e.g. the app was closed and reopened without clicking Stop).
                if "share_url" not in ss and _share_running():
                    _u = _read_share_url()
                    if _u:
                        ss["share_url"] = _u
                        ss.setdefault("share_live", False)
                        st.info("Reconnected to a share that was already running. "
                                "Use Re-check or Stop sharing below.")
                if "share_url" in ss:
                    if ss.get("share_live"):
                        st.success("🟢 Sharing is live — send this link and password to your "
                                   "collaborator:")
                    else:
                        st.warning("⏳ Link created, but the public tunnel is still warming up "
                                   "(Cloudflare free tier can take a minute). Send it once the "
                                   "check below goes green, or Re-check.")
                    st.code(ss["share_url"], language=None)
                    st.code(f"password: {ss.get('share_pw_active', '')}", language=None)
                    st.caption("The link works only while your PC and this app stay on. "
                               "The URL changes if you stop and start again.")
                    bc1, bc2 = st.columns(2)
                    if not ss.get("share_live") and bc1.button("🔄 Re-check", key="share_recheck"):
                        ss["share_live"] = _url_live(ss["share_url"])
                        st.rerun()
                    if bc2.button("⏹ Stop sharing", key="share_stop"):
                        _stop_share()
                        ss.pop("share_url", None)
                        ss.pop("share_live", None)
                        st.rerun()
                else:
                    have = "analysis" in ss
                    gmode = ss.get("grp_mode")
                    gav = ss.get("grp_a_vals") or []
                    gbv = ss.get("grp_b_vals") or []
                    default_name = (f"{ss['analysis'].get('label_b', 'B')} vs "
                                    f"{ss['analysis'].get('label_a', 'A')}") if have else ""
                    share_name = st.text_input("Title collaborators will see",
                                               value=default_name, key="share_name")
                    share_pw = st.text_input("Password", value=ss.get("share_pw_active", "MacCoss2026"),
                                             key="share_pw")
                    if not have:
                        st.info("Run a comparison first (Differential tab), then share it here.")
                    elif gmode != "By metadata column" or not (gav and gbv):
                        st.info("Sharing needs a **metadata-column** comparison — in the sidebar use "
                                "*Group by* a column with values in Group A and Group B, then Run limma.")
                    else:
                        a = ss["analysis"]
                        st.caption(f"Will share **{share_name or 'PRISM analysis'}** — "
                                   f"{a.get('label_b')} vs {a.get('label_a')} · {level} · "
                                   f"FDR<{q_cut:g}, |log2FC|≥{lfc_cut:g}")
                        if st.button("🚀 Share this analysis", type="primary", key="share_go"):
                            cfg = {
                                "password": share_pw.strip() or "MacCoss2026",
                                "display_name": share_name.strip() or "PRISM analysis",
                                "level": level, "out_dir": out_dir, "clin_path": clin_path,
                                "clin_key": clin_key,
                                "comparison": {
                                    "grouping_column": ss.get("grp_col"),
                                    "group_a_values": list(gav), "group_b_values": list(gbv),
                                    "suggested_covariates": list(ss.get("cov_select") or []),
                                    "label_a": a.get("label_a"), "label_b": a.get("label_b"),
                                },
                                "q_cut": float(q_cut), "lfc_cut": float(lfc_cut),
                            }
                            ss["share_pw_active"] = cfg["password"]
                            import time
                            with st.spinner("Starting the shared app and opening a public "
                                            "tunnel (~20–60s)…"):
                                try:
                                    _start_share(cfg)
                                except Exception as _e:  # noqa: BLE001
                                    st.error(f"Couldn't start sharing: {_e}")
                                    st.stop()
                                url = ""
                                for _ in range(30):
                                    time.sleep(2)
                                    url = _read_share_url()
                                    if url:
                                        break
                                live = False
                                if url:
                                    for _ in range(10):  # let the Cloudflare edge become routable
                                        if _url_live(url):
                                            live = True
                                            break
                                        time.sleep(3)
                            if url:
                                ss["share_url"] = url
                                ss["share_live"] = live
                                st.rerun()
                            else:
                                st.error("Couldn't get the tunnel URL in time. Check that "
                                         "cloudflared is installed at "
                                         "`%USERPROFILE%\\Apps\\cloudflared.exe` and try again.")
            st.divider()

        st.caption("Build one self-contained HTML report to share with a collaborator — no data, "
                   "no install, opens in any browser with the charts interactive. Tick the "
                   "analyses to include (only those you've run this session are available).")
        rtitle = st.text_input("Report title", value="PRISM analysis report", key="rep_title")
        rintro = st.text_area("Intro / notes (optional)", value="", key="rep_intro", height=80)

        available = {
            "Differential": "analysis" in ss,
            "Interactive box plots": "analysis" in ss,
            "PCA": "pca_result" in ss,
            "EV markers": "ev_result" in ss and data["level"] == "protein",
            "Cryptic": "cryptic_result" in ss,
            "Longitudinal": "long_result" in ss,
            "Heterogeneity": "het_result" in ss,
        }
        chosen = []
        for name, ok in available.items():
            label = name if ok else f"{name}  — run it first"
            if st.checkbox(label, value=ok, disabled=not ok, key=f"rep_inc_{name}") and ok:
                chosen.append(name)

        if not chosen:
            st.info("Run at least one analysis (Differential, Longitudinal, or Heterogeneity), "
                    "then tick it here to include it.")
        else:
            sections = []
            if "Differential" in chosen:
                a = ss["analysis"]
                res = a["res"]
                la_, lb_ = a.get("label_a", "A"), a.get("label_b", "B")
                sig, p_thresh_r = _significance(res, sig_basis, q_cut, lfc_cut)
                _bw = "raw p" if sig_basis == "raw" else "FDR"
                vdf = res.assign(_nlp=-np.log10(res["P.Value"].clip(lower=1e-300)),
                                 _sig=np.where(sig, "significant", "ns"))
                vfig = px.scatter(vdf, x="logFC", y="_nlp", color="_sig", height=430,
                                  color_discrete_map={"significant": "#d62728", "ns": "#b0b0b0"},
                                  labels={"logFC": f"log2 FC ({lb_}/{la_})", "_nlp": "-log10 P"},
                                  title=f"Differential: {lb_} vs {la_}")
                _add_volcano_guides(vfig, p_thresh_r, lfc_cut, sig_basis, q_cut)
                top = res.sort_values("P.Value").head(25)[["logFC", "P.Value", "adj.P.Val"]].round(3)
                sections.append({"heading": f"Differential — {lb_} vs {la_}",
                                 "paras": [f"{int(sig.sum())} significant of {len(res):,} features at "
                                           f"{_bw}<{q_cut:g}, |log2FC| >= {lfc_cut:g}."],
                                 "figs": [vfig], "table": top})
            if "Interactive box plots" in chosen:
                a = ss["analysis"]
                la_, lb_ = a.get("label_a", "A"), a.get("label_b", "B")
                pre = a["res"].sort_values("P.Value").index[0] if len(a["res"]) else None
                frag, cap_note = _feature_explorer_html(
                    expr_log2, a["group_a"], a["group_b"], la_, lb_, annot,
                    preselect=pre, block_id="pdxexp")
                paras = ["Search any feature by accession or peptide sequence and its box plot "
                         f"({lb_} vs {la_}) renders in your collaborator's browser — no data file, "
                         "no install. Opens on the top hit."]
                if cap_note:
                    paras.append(cap_note)
                sections.append({"heading": "Explore any feature — interactive box plots",
                                 "paras": paras, "html": frag})
            if "PCA" in chosen:
                p = ss["pca_result"]
                sections.append({"heading": "PCA", "paras": [p["summary"]], "figs": [p["fig"]]})
            if "EV markers" in chosen:
                e = ss["ev_result"]
                sections.append({"heading": "EV markers", "paras": [e["summary"]], "figs": [e["fig"]]})
            if "Cryptic" in chosen:
                cy = ss["cryptic_result"]
                sections.append({"heading": "Cryptic peptides — current contrast",
                                 "paras": [cy["summary"]],
                                 "figs": [cy["fig"]], "table": cy.get("table")})
            if "Longitudinal" in chosen:
                lo = ss["long_result"]
                lres, info, tcol = lo["res"], lo["info"], lo["time_col"]
                lv = lres.assign(_nlp=-np.log10(lres["P.Value"].clip(lower=1e-300)))
                lfig = px.scatter(lv, x="slope", y="_nlp", height=430,
                                  labels={"slope": f"slope per {tcol}", "_nlp": "-log10 P"},
                                  title="Within-subject time trend")
                lt = lres.head(25)[["slope", "P.Value", "adj.P.Val"]].round(4)
                sections.append({"heading": "Longitudinal — within-subject time trend",
                                 "paras": [f"{info['n_subjects']} subjects, {info['n_samples']} samples, "
                                           f"residual df {info['residual_df']}. Exploratory — small "
                                           "cohorts rarely clear FDR."],
                                 "figs": [lfig], "table": lt})
            if "Heterogeneity" in chosen:
                h = ss["het_result"]
                htab2, focus2 = h["tab"], h["focus"]
                seen2 = htab2[htab2["n_detected"] > 0]
                sub2 = int((seen2["prevalence"] < 0.5).sum())
                bb = pd.cut(seen2["prevalence"], [0, .25, .5, .75, .99, 1.001], right=False,
                            labels=["1-25%", "25-50%", "50-75%", "75-99%", "100%"]) \
                    .value_counts().sort_index()
                hfig = px.bar(x=bb.index.astype(str), y=bb.to_numpy(), height=380,
                              labels={"x": f"detected in this fraction of {focus2} individuals",
                                      "y": "features"}, title=f"Detection prevalence — {focus2}")
                ht = (seen2.sort_values("prevalence")
                      .head(25)[["prevalence", "n_detected", "pooling_dilution"]].round(2))
                sections.append({"heading": f"Heterogeneity — {focus2}",
                                 "paras": [f"{len(h['indiv'])} individuals, {len(h['pools'])} pool/ref "
                                           f"samples. {sub2:,} features detected in <50% of {focus2} "
                                           "individuals — subset markers a pooled sample dilutes toward "
                                           "baseline (heterogeneity, not absence)."],
                                 "figs": [hfig], "table": ht})
            report_html = _build_html_report(rtitle, rintro, sections)
            st.download_button("⬇ Download HTML report", report_html.encode("utf-8"),
                               file_name="prism_report.html", mime="text/html", type="primary",
                               key="rep_dl")
            st.caption(f"{len(sections)} section(s) · self-contained, interactive charts embedded.")


if __name__ == "__main__":
    main()
