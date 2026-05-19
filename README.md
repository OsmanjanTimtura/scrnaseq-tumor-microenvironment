# scRNA-seq Tumor Microenvironment Analysis — Lung (LuCA)

🚧 **In progress** — repo bootstrapped May 14, 2026; target completion 6–8 weeks.

End-to-end single-cell RNA-seq pipeline applied to the **LuCA core atlas of lung cancer cell states** (Salcher et al. 2022, *Cancer Cell*; 892K cells, 19 integrated NSCLC studies). Uses the 2026-standard pharma comp-bio stack — **scanpy + scVI** plus direct HDF5 reads via `h5py` for memory-efficient subsampling — for QC, embedding evaluation, cell-type annotation, and downstream compositional analysis of LUAD vs LUSC.

**Author:** Osmanjan Timtura, Ph.D. · [LinkedIn](https://www.linkedin.com/in/osmanjan-timtura-b08316245/) · [GitHub](https://github.com/OsmanjanTimtura)

---

## Project plan

See [PROJECT_PLAN.md](PROJECT_PLAN.md) for the full 6–8 week scope, dataset choice, biological question, phases, risks, and what's deliberately out of scope.

**Short version:** Salcher et al. 2022 lung cancer atlas (LuCA) → stratified subsample (~100–200K cells) → QC + evaluation of the published scVI embedding → cell-type annotation against published LuCA labels → per-patient cell-state composition → LightGBM classifier predicting LUAD vs LUSC histology with leave-one-study-out cross-validation. **Stretch:** within LUAD only, never-smoker vs smoker compositional comparison (Sherlock-Lung adjacent).

## Status

| Phase | Status |
|---|---|
| 1. Setup + LuCA subsample (h5py streaming, no Census-API dep) | ✅ done |
| 2. QC + filter + normalize (notebook 01) | ✅ done |
| 3. scVI integration (notebook 02) | 🚧 in progress |
| 4. Cell-type annotation vs LuCA benchmark (notebook 03) | ⏳ pending |
| 5. LUAD vs LUSC compositional classifier (notebook 04) | ⏳ pending |
| 6. Writing + polish | ⏳ pending |
| 7. (Stretch) never-smoker vs smoker analysis | ⏳ pending |

## Getting started

```bash
git clone https://github.com/OsmanjanTimtura/scrnaseq-tumor-microenvironment.git
cd scrnaseq-tumor-microenvironment

# Recommended: conda env so h5py / scanpy stay on a known numpy ABI
conda create -n scrnaseq -c conda-forge python=3.11 scanpy anndata leidenalg \
    python-igraph jupyter ipykernel h5py pandas numpy matplotlib seaborn \
    scikit-learn pytest -y
conda activate scrnaseq
pip install scvi-tools  # not on conda-forge yet

# One-time data setup (see data/README.md for the full guide):
#   1. Download luca_core.h5ad (~13 GB) from the CELLxGENE portal.
#   2. Build the 100K subsample (~770 MB, ~3 min):
python -m src.subsample

jupyter notebook notebooks/01_data_download_and_qc.ipynb
```

## Tools

`Python 3.11` · `scanpy` · `anndata` · `scvi-tools` · `h5py` · `leiden` · `lightgbm` · `shap` · `scikit-learn` · `pandas` · `numpy` · `matplotlib` · `seaborn` · `Jupyter` · `pytest`

## What this repo demonstrates

- Streaming subsample of a multi-study scRNA-seq atlas straight out of HDF5
  (sparse-CSR slicing via `h5py` + `anndata.io.sparse_dataset`) — no need to
  materialize the 13 GB file in RAM
- Stratified sampling preserving per-study and per-disease balance (19 cohorts
  × 5 disease labels)
- QC filtering against an atlas that already carried pre-computed metrics
  (`n_genes_by_counts`, `pct_counts_mito`, `doublet_status`) — comparing our
  thresholds against the atlas's own
- scVI batch correction with `study` as the batch key
- Marker-gene-based cell-type annotation, compared against LuCA's published
  labels for benchmarking
- Cell-state composition analysis at the patient level
- LightGBM + SHAP classifier predicting NSCLC histology from cell-state
  composition with **leave-one-study-out cross-validation** (the right CV for
  multi-study data)
- Honest limitations section

## Provenance

Phase 4 of my self-directed computational-biology transition plan. Third pinned GitHub repo after `codon-discovery-pca-kmeans` and `cell-confluency-segmentation`. The dataset choice (LuCA) is deliberately adjacent to the **Sherlock-Lung Study at NCI / FNLCR** — the next-tier-up research environment I'd plausibly target as a future career move.

## License

[MIT](LICENSE) — reuse freely.
