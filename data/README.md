# Data download instructions — LuCA lung cancer atlas

This project works on the **LuCA core atlas of lung cancer cell states**
(Salcher et al. 2022, *Cancer Cell*; DOI [10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008)).
The published atlas integrates 19 NSCLC scRNA-seq studies into a single
harmonized AnnData (892,296 cells × 17,764 genes, ~13 GB on disk).

Raw data is **not** committed to this repo (GitHub free tier caps at 100 MB
per file; scientific atlases this size don't belong in a code repo). The
analyses below run on a stratified 100K-cell subsample produced locally from
the published h5ad.

---

## TL;DR — two-step setup

```bash
# 1. Download luca_core.h5ad from the CELLxGENE portal (see Option A below).
#    Save anywhere with ≥15 GB free; default expected at data/raw/luca_core.h5ad.

# 2. Run the subsampler. Produces data/luca_subsample_100k.h5ad (~770 MB).
python -m src.subsample
```

Notebook 01 picks up from the subsample. Total wall time: ~3 minutes after the
download finishes.

---

## A. Get the full atlas

### A.1: Direct download from CELLxGENE (recommended)

Cross-platform, no extra dependencies:

1. Visit https://cellxgene.cziscience.com/
2. Search for **"LuCA"** or **"core atlas of lung cancer cell states"**
3. Open the Salcher 2022 dataset → **Download** → choose **h5ad**
4. Save as `luca_core.h5ad`. Default location: `data/raw/luca_core.h5ad`. If
   you keep it elsewhere (e.g. on a separate drive), set the environment
   variable `LUCA_CORE_PATH` before running `python -m src.subsample`.

File size: ~13 GB.

### A.2 (alternative): `cellxgene-census` Python package — Linux / macOS only

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        obs_value_filter=(
            "tissue_general == 'lung' "
            "and disease in ['lung adenocarcinoma', "
            "'squamous cell lung carcinoma']"
        ),
    )
```

**Windows caveat:** `cellxgene-census` depends on `tiledbsoma`, which as of
May 2026 does **not** publish a Windows wheel. The package fails to install on
Windows from PyPI. Use Option A.1 above. (See the upstream issue tracker for
status updates.)

### A.3 (alternative): direct download via wget

```bash
mkdir -p data/raw
wget -O data/raw/luca_core.h5ad \
    "https://datasets.cellxgene.cziscience.com/<dataset_uuid>.h5ad"
```

Substitute `<dataset_uuid>` with the LuCA dataset UUID shown on the CELLxGENE
portal page. Saves a round-trip through the browser.

---

## B. Build the 100K subsample

The full atlas's `.X` carries the published normalized values; the **raw
counts live in `/layers/count`**. Notebook 02 (scVI) requires raw counts.
`src/subsample.py` slices `/layers/count` directly via h5py — no need to load
the 13 GB file into memory.

### B.1: via CLI

```bash
# uses data/raw/luca_core.h5ad by default
python -m src.subsample
```

Or with a custom source path:

```bash
# Windows PowerShell
$env:LUCA_CORE_PATH = "D:\singleCell_rnaSeq\dataset\luca_core.h5ad"
python -m src.subsample

# bash / zsh
LUCA_CORE_PATH=~/data/luca_core.h5ad python -m src.subsample
```

### B.2: programmatically

```python
from src.subsample import subsample_luca

adata = subsample_luca(
    source_h5ad="data/raw/luca_core.h5ad",
    n_target=100_000,
    seed=42,
)
adata.write_h5ad("data/luca_subsample_100k.h5ad", compression="gzip")
```

The subsampler is stratified by `study × disease` so each of the 19 source
cohorts and 5 disease labels (LUAD, LUSC, NSCLC unspecified, normal, and
LUAD-adjacent normal) is represented proportionally. With `seed=42` and
`n_target=100_000`, the result is 92,452 cells (`100_000 // 30 strata = 3,333`
cells per stratum × 30 strata, less the strata that hold fewer than 3,333
cells).

Wall time: ~2–3 minutes on SSD, ~2 GB peak RAM.

---

## Resulting files

| Path                                                  | Size      | Source                                              | Notes                                                                  |
| ----------------------------------------------------- | --------- | --------------------------------------------------- | ---------------------------------------------------------------------- |
| `<LUCA_CORE_PATH>` (default `data/raw/luca_core.h5ad`) | ~13 GB    | CELLxGENE portal                                    | Full published atlas. Keep outside repo.                               |
| `data/luca_subsample_100k.h5ad`                       | ~770 MB   | `python -m src.subsample`                           | 92K cells, raw counts in `.X`, atlas-normalized values as a layer.     |
| `data/processed/luca_filtered.h5ad`                   | ~1.1 GB   | `notebooks/01_data_download_and_qc.ipynb`           | Filtered + log-normalized; raw counts preserved in `layers["counts"]`. |
| `models/scvi/`                                        | ~50 MB    | `notebooks/02_scvi_integration.ipynb`               | scVI checkpoint. Gitignored.                                           |

---

## Expected metadata

```python
print(adata.obs.columns.tolist())
# ['study', 'dataset', 'disease', 'sample', 'donor_id', 'cell_type',
#  'cell_type_major', 'ann_coarse', 'ann_fine', 'ever_smoker', 'sex', 'age',
#  'tumor_stage', 'uicc_stage', 'origin', 'origin_fine', 'platform', 'assay',
#  'n_genes_by_counts', 'total_counts', 'pct_counts_mito',
#  'total_counts_mito', 'doublet_status', 'EGFR_mutation', 'KRAS_mutation',
#  'TP53_mutation', 'ALK_mutation', 'BRAF_mutation']
```

Key fields used downstream:

- `study` — source publication (19 cohorts; primary batch variable for scVI)
- `disease` — LUAD, LUSC, NSCLC (unspecified), normal, lung-adjacent normal
- `cell_type_major` — LuCA's published coarse labels; benchmark in notebook 03
- `ever_smoker` — used in the stretch never-smoker vs smoker analysis
- `n_genes_by_counts`, `total_counts`, `pct_counts_mito` — **pre-computed QC
  metrics**; notebook 01 filters on these directly rather than recomputing
- `doublet_status` — pre-computed; non-singlets are dropped in notebook 01
- `is_highly_variable` (in `var`) — atlas-flagged HVGs

---

## Disk + RAM planning

| Artifact                  | Disk     | Working RAM (load + process)         |
| ------------------------- | -------- | ------------------------------------ |
| Full atlas (`luca_core`)  | ~13 GB   | 50+ GB (loaded naively — don't)      |
| Subsampling step          | —        | ~2 GB peak                           |
| 100K subsample            | ~770 MB  | ~2 GB                                |
| QC-filtered processed     | ~1.1 GB  | ~3 GB                                |
| scVI training (CPU)       | —        | ~4–6 GB                              |

16 GB system RAM is enough for the full pipeline once the subsample is built.

---

## Why we don't use cellxgene-census on Windows

The original project plan reached for `cellxgene-census` for clean query-based
download. Its underlying TileDB binding (`tiledbsoma`) ships pre-built wheels
for Linux and macOS but not for Windows; `pip install cellxgene-census` on
Windows triggers a C++ source build that fails on the standard MSVC toolchain.
Rather than patch the dependency, we download the published h5ad directly and
subsample locally — same scientific result, less environment friction.

If you're on Linux or macOS and want the Census-API path, `cellxgene-census`
installs cleanly there and Option A.2 above will work.
