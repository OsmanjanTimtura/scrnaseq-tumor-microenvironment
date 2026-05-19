"""Quality-control helpers for scRNA-seq data."""
from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def compute_qc_metrics(adata: AnnData, mt_prefix: str = "MT-") -> AnnData:
    """Compute standard per-cell QC metrics in-place.

    Adds to adata.obs:
        - n_genes_by_counts
        - total_counts
        - pct_counts_mt
    """
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def filter_cells(
    adata: AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mt: float = 20.0,
) -> Tuple[AnnData, pd.Series]:
    """Filter low-quality cells. Returns the filtered AnnData and a Series
    indicating which cells were kept (boolean, indexed by original cell IDs).

    Default thresholds are reasonable starting points for 10X Chromium data;
    inspect the QC violin plots first before committing.
    """
    keep = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    )
    return adata[keep, :].copy(), keep


def qc_summary(adata: AnnData) -> pd.DataFrame:
    """Return a per-sample QC summary table.

    Assumes adata.obs has a 'sample' or 'patient' column; falls back to whole-data summary if not.
    """
    sample_col = next(
        (c for c in ("sample", "patient", "donor_id", "Sample") if c in adata.obs.columns),
        None,
    )
    if sample_col is None:
        rows = [{
            "sample": "ALL",
            "n_cells": adata.n_obs,
            "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
            "median_counts": float(np.median(adata.obs["total_counts"])),
            "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"])),
        }]
    else:
        rows = []
        for s, sub in adata.obs.groupby(sample_col):
            rows.append({
                "sample": s,
                "n_cells": len(sub),
                "median_genes": float(np.median(sub["n_genes_by_counts"])),
                "median_counts": float(np.median(sub["total_counts"])),
                "median_pct_mt": float(np.median(sub["pct_counts_mt"])),
            })
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)
