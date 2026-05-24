# Project Plan — scRNA-seq Tumor Microenvironment Analysis (Lung)

**Author:** Osmanjan Timtura

## Scope

A self-contained portfolio piece demonstrating end-to-end single-cell RNA-seq analysis using the 2026-standard pharma comp-bio stack (scanpy + scVI). **Not a publishable paper.** The goal is one clean GitHub repo that proves competence in the workflow recruiters at pharma comp-bio teams screen for.

Specifically the repo should demonstrate:

- Loading a large multi-study public scRNA-seq atlas (the 1.2M-cell LuCA lung-cancer atlas)
- Stratified subsampling that preserves per-study and per-subtype balance
- QC filtering on the subsample
- Embedding-quality evaluation of the published scVI integration (vs. building one from scratch on a smaller dataset)
- Cell-type annotation using marker-gene signatures, compared against LuCA's published labels
- Per-patient compositional analysis
- ML classifier predicting tumor histology (LUAD vs LUSC) from cell-state composition
- Reproducible, version-pinned environment with honest limitations

## Dataset

**Primary:** Salcher et al. 2022, *Cancer Cell* — the **"core atlas of lung cancer cell states" (LuCA)**. ~1.2M cells integrating 29 published non-small-cell lung cancer (NSCLC) scRNA-seq studies into a single harmonized atlas. Covers LUAD (adenocarcinoma), LUSC (squamous-cell carcinoma), and normal adjacent lung tissue. Pre-integrated scVI embedding included.

**Why this one:**
- Cross-cohort integration is already done — the value-add of this portfolio piece becomes "I can evaluate and use a large public atlas," not "I can train scVI from scratch" (still demonstrate scVI competence in the appendix)
- Multi-study scale is what real pharma comp-bio teams work at; smaller atlases under-represent the messiness
- Adjacent to Sherlock-Lung at NCI — the next-tier-up role I'd plausibly target in 12–18 months
- Public on CELLxGENE Census; downloadable via the `cellxgene-census` Python package without registration
- Both major NSCLC histologies (LUAD, LUSC) are represented in sufficient numbers for a supervised-classification question
- Some sub-cohorts include smoking-history metadata (stretch analysis)

**Subsampling strategy:** the full atlas requires ~50 GB RAM and is impractical on a laptop. Default pipeline loads a **100,000–200,000 cell stratified subsample** preserving per-study and per-histology balance. A `--full` flag in the data-loader honors the full atlas for users with appropriate compute.

**Fallback datasets if LuCA proves too large or its CELLxGENE Census endpoint changes:**
- Kim et al. 2020 *Nat Commun* lung adenocarcinoma scRNA-seq (single study, ~200K cells)
- Lambrechts et al. 2018 *Nat Med* lung cancer atlas (older but smaller, ~50K cells)
- Wu et al. 2021 *Nature Genetics* breast cancer atlas (the original v0 plan; switch if lung doesn't work)

## Biological question

**Primary:** Across the LuCA cohort, do cell-state composition signatures distinguish LUAD from LUSC, and which cell-state features are most predictive of histology?

This is small-scope but legitimate. It demonstrates: cell-type-composition analysis at scale, multi-study handling, supervised ML on compositional features. Cross-validation must be **leave-one-study-out** rather than leave-one-patient-out because cohort effects in multi-study data are the dominant nuisance signal.

**Stretch:** Within LUAD only, do never-smokers vs smokers have distinct immune-cell-composition signatures? This is the Sherlock-Lung-adjacent question. Requires that the smoking metadata in the LuCA subset is sufficient — flagged as stretch goal because it depends on per-cohort metadata completeness.

## Phases — 6 to 8 weeks

| Phase | Weeks | Deliverable |
|---|---|---|
| 1. Setup + data download | 1 | Env installed (scanpy + scvi-tools + cellxgene-census), subsample h5ad on disk |
| 2. QC + embedding evaluation | 2 | Filtered AnnData, evaluation of the published scVI embedding (kBET-style mixing metrics) |
| 3. Cell-type annotation | 3 | Major cell types annotated; comparison against LuCA's published labels |
| 4. Compositional analysis (LUAD vs LUSC) | 4–5 | Per-patient cell-state proportions + classifier with leave-one-study-out CV |
| 5. Writing + polish | 6 | README, figures, code cleanup, pytest tests |
| 6. (Optional) Smoking-status stretch | 7–8 | Within-LUAD never-smoker vs smoker compositional comparison |

## Repo structure

```
scrnaseq-tumor-microenvironment/
├── README.md
├── PROJECT_PLAN.md                    (this file)
├── requirements.txt
├── LICENSE
├── .gitignore
├── data/
│   ├── README.md                      (CELLxGENE Census download instructions)
│   └── .gitkeep
├── notebooks/
│   ├── 01_data_download_and_qc.ipynb  (uses cellxgene-census; subsamples to 100-200K)
│   ├── 02_scvi_integration.ipynb      (evaluate published embedding + optional retrain)
│   ├── 03_cell_type_annotation.ipynb
│   └── 04_compositional_analysis.ipynb
├── src/
│   ├── __init__.py
│   ├── data_io.py
│   ├── qc.py
│   ├── integration.py
│   └── annotation.py
├── figures/
└── tests/
    └── test_qc.py
```

## Compute requirements

- **Subsampled (100–200K cells, default): CPU laptop is fine.** scVI re-training on the subsample takes ~30–60 min on CPU, faster on GPU.
- **Full atlas (1.2M cells): GPU required.** ~50 GB RAM, plus a GPU with 12+ GB VRAM for scVI training. Cloud rental: Colab Pro ($10/month) or a single-day AWS p3 instance (~$3).
- Total disk for raw + processed data: ~10–15 GB.
- Skip full-atlas re-training if you don't have a GPU; the published scVI embedding is good enough to use directly.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| LuCA endpoint changes on CELLxGENE Census | `cellxgene-census` is version-pinned in requirements.txt; fallback datasets listed above |
| 1.2M cells overwhelms laptop | Subsample to 100–200K with stratified sampling; documented in `src/data_io.py` |
| Published scVI embedding has artifacts | QC the embedding (Phase 2) and report what we find honestly |
| LUAD vs LUSC classification too easy (trivially separable) | If AUROC > 0.95, switch focus to LUAD subtyping or within-LUAD smoker-status |
| Smoking metadata sparse in LuCA | Stretch goal only — drop if Phase 1 reveals coverage too thin |
| Schedule slips beyond 8 weeks | Drop the smoking stretch first; core 6-week scope still produces a complete repo |

## What this repo does NOT attempt

- Not a publishable paper. (The SLE plan is separate.)
- Not a Sherlock-Lung-style mutational analysis — this is transcriptome only.
- No CAR-T / therapy-response framing — stays diagnostic / mechanistic, outside uBriGene PIIA scope.
- No multi-omics integration; single-modality (transcriptome).

## Provenance

This is Phase 4 of the broader Computational Biology Transition Plan
(`CompBio_Transition_Stepwise_Plan.docx`). Third pinned GitHub repo after
`codon-discovery-pca-kmeans` and `cell-confluency-segmentation`. Demonstrates
the specific tool stack (scanpy + scVI + CELLxGENE Census) that pharma
comp-bio teams screen for, on a dataset adjacent to NCI Sherlock-Lung
research — the kind of substrate a future Computational Scientist II
posting at FNLCR / NCI would value.
