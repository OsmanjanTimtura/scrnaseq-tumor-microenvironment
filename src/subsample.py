"""LuCA atlas subsampling.

The full Salcher 2022 LuCA atlas (``luca_core.h5ad``, ~13 GB) carries two huge
layers — ``/layers/count`` (raw counts) and ``/layers/counts_length_scaled`` —
that, together with ``/X`` (atlas-normalized values), don't fit in 16 GB of
RAM if loaded naively.

This module reads the source file with **h5py** directly so we can:

1. Inspect ``/obs`` cheaply and choose cell indices stratified by study × disease
   *before* materializing any matrix.
2. Slice raw counts out of ``/layers/count`` for the selected indices — those
   are what scVI expects in notebook 02.
3. Optionally keep the atlas's pre-normalized ``/X`` values in a layer so we
   can compare our own normalization to LuCA's published one.

The default subsample of 100,000 cells across 19 studies × 5 diseases produces
a ~770 MB h5ad that fits comfortably on a 16 GB workstation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import anndata as ad
import h5py
import numpy as np
import pandas as pd


# Columns kept from /obs by default. Includes everything notebooks 01-04 use
# (study/disease for stratification, pre-computed QC, clinical metadata, and
# the published cell-type labels we benchmark against in notebook 03).
DEFAULT_OBS_COLUMNS: tuple[str, ...] = (
    "study", "dataset", "disease", "sample", "donor_id",
    "cell_type", "cell_type_major", "ann_coarse", "ann_fine",
    "ever_smoker", "sex", "age", "tumor_stage", "uicc_stage",
    "origin", "origin_fine", "platform", "assay",
    "n_genes_by_counts", "total_counts", "pct_counts_mito",
    "total_counts_mito", "doublet_status",
    "EGFR_mutation", "KRAS_mutation", "TP53_mutation",
    "ALK_mutation", "BRAF_mutation",
)


def _read_h5_col(grp: h5py.Group, name: str):
    """Read an h5ad obs/var column; handles categoricals and byte strings.

    AnnData stores categoricals as a sub-group with ``codes`` and ``categories``
    datasets, and string columns as either fixed-width bytes (dtype kind "S")
    or variable-length objects.
    """
    item = grp[name]
    if isinstance(item, h5py.Group):  # categorical
        codes = item["codes"][:]
        cats = item["categories"][:]
        if cats.dtype.kind in ("S", "O"):
            cats = np.array(
                [c.decode() if isinstance(c, bytes) else c for c in cats]
            )
        return pd.Categorical.from_codes(codes, cats)
    arr = item[:]
    if arr.dtype.kind == "S":
        arr = np.array([x.decode() for x in arr])
    return arr


def _import_sparse_dataset():
    """Import sparse_dataset from whichever anndata module exposes it.

    anndata renamed the path from anndata.experimental to anndata.io between
    minor releases; we accept either.
    """
    try:
        from anndata.io import sparse_dataset  # newer
        return sparse_dataset
    except ImportError:
        from anndata.experimental import sparse_dataset  # older
        return sparse_dataset


def subsample_luca(
    source_h5ad: Path | str,
    n_target: int = 100_000,
    seed: int = 42,
    obs_columns: Optional[Sequence[str]] = None,
    stratify_by: Sequence[str] = ("study", "disease"),
    include_atlas_norm: bool = True,
) -> ad.AnnData:
    """Stratified subsample of the LuCA full atlas, with raw counts in ``.X``.

    Parameters
    ----------
    source_h5ad
        Path to ``luca_core.h5ad`` from the CELLxGENE portal (see
        ``data/README.md`` for download instructions). ~13 GB.
    n_target
        Approximate cell count to keep. Actual count is per-stratum capped
        (``n_target // n_strata`` cells per stratum), so the final count is
        usually a bit below ``n_target``.
    seed
        RNG seed. Default 42 produces the same 92,452-cell selection on
        repeated runs.
    obs_columns
        Which ``/obs`` columns to retain. Defaults to :data:`DEFAULT_OBS_COLUMNS`.
        Missing columns are silently skipped.
    stratify_by
        ``/obs`` column names that define the strata. Default
        ``("study", "disease")`` balances the 19 source studies across the 5
        disease labels (LUAD, LUSC, NSCLC, normal, etc.).
    include_atlas_norm
        If True, also slice ``/X`` (the atlas's pre-normalized values) into
        ``layers["normalized_atlas"]``. Adds ~30-90s and ~400 MB to the result.

    Returns
    -------
    AnnData
        ``.X`` is raw counts (CSR sparse, float32). ``.obs`` carries the
        selected columns; ``.var`` carries all gene metadata from the source.

    Raises
    ------
    FileNotFoundError
        If ``source_h5ad`` doesn't exist.
    KeyError
        If ``/layers/count`` is missing from the source file (would indicate
        the file is not the canonical LuCA atlas).
    """
    source_h5ad = Path(source_h5ad)
    if not source_h5ad.exists():
        raise FileNotFoundError(
            f"Source not found: {source_h5ad}\n"
            f"See data/README.md for download instructions."
        )

    rng = np.random.default_rng(seed)
    cols = list(obs_columns) if obs_columns is not None else list(DEFAULT_OBS_COLUMNS)
    sparse_dataset = _import_sparse_dataset()

    with h5py.File(source_h5ad, "r") as f:
        if "layers" not in f or "count" not in f["layers"]:
            raise KeyError(
                "Source h5ad does not contain /layers/count. "
                "Confirm the file is the LuCA atlas (luca_core.h5ad) and not "
                "a different export."
            )

        # ---- /obs (small, full) ----
        obs_grp = f["obs"]
        idx_key = obs_grp.attrs.get("_index", "_index")
        obs_data = {}
        for c in cols:
            if c in obs_grp:
                try:
                    obs_data[c] = _read_h5_col(obs_grp, c)
                except Exception:
                    pass  # tolerate odd column encodings
        obs_index = _read_h5_col(obs_grp, idx_key)
        obs_full = pd.DataFrame(
            obs_data, index=pd.Index(obs_index, name="cell_id")
        )

        # ---- stratified sample ----
        missing = [c for c in stratify_by if c not in obs_full.columns]
        if missing:
            raise KeyError(f"stratify_by columns not in obs: {missing}")
        strat = obs_full[stratify_by[0]].astype(str)
        for c in stratify_by[1:]:
            strat = strat + "|" + obs_full[c].astype(str)
        counts = strat.value_counts()
        per_stratum = max(1, n_target // len(counts))
        all_pos = np.arange(len(obs_full))
        picks = [
            rng.choice(
                all_pos[(strat == s).values],
                size=min(per_stratum, n),
                replace=False,
            )
            for s, n in counts.items()
        ]
        sel = np.sort(np.concatenate(picks))

        # ---- raw counts → .X ----
        X_raw = sparse_dataset(f["layers/count"])[sel, :]
        if hasattr(X_raw, "to_memory"):
            X_raw = X_raw.to_memory()

        # ---- atlas-normalized → layer (optional) ----
        X_norm = None
        if include_atlas_norm and "X" in f:
            X_norm = sparse_dataset(f["X"])[sel, :]
            if hasattr(X_norm, "to_memory"):
                X_norm = X_norm.to_memory()

        # ---- /var (small, full) ----
        var_grp = f["var"]
        var_idx_key = var_grp.attrs.get("_index", "_index")
        var_data = {}
        for c in var_grp.keys():
            if c == var_idx_key:
                continue
            try:
                var_data[c] = _read_h5_col(var_grp, c)
            except Exception:
                pass
        var_index = _read_h5_col(var_grp, var_idx_key)
        var_df = pd.DataFrame(
            var_data, index=pd.Index(var_index, name="gene_id")
        )

    obs_sub = obs_full.iloc[sel].copy()
    adata = ad.AnnData(X=X_raw, obs=obs_sub, var=var_df)
    if X_norm is not None:
        adata.layers["normalized_atlas"] = X_norm
    return adata


def main() -> None:  # pragma: no cover - CLI convenience wrapper
    """``python -m src.subsample`` — produce the default 100K subsample.

    Reads from ``data/raw/luca_core.h5ad`` (or path in env LUCA_CORE_PATH) and
    writes to ``data/luca_subsample_100k.h5ad``.
    """
    import os
    import time

    repo_root = Path(__file__).resolve().parent.parent
    source = Path(
        os.environ.get("LUCA_CORE_PATH", repo_root / "data" / "raw" / "luca_core.h5ad")
    )
    out_path = repo_root / "data" / "luca_subsample_100k.h5ad"

    print(f"Source: {source}")
    print(f"Output: {out_path}")
    t0 = time.time()
    adata = subsample_luca(source, n_target=100_000, seed=42)
    print(adata)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(
        f"\nWrote {out_path}  "
        f"({out_path.stat().st_size / 1e9:.2f} GB, {time.time() - t0:.0f}s)"
    )


if __name__ == "__main__":  # pragma: no cover
    main()
