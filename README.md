# Companion code — *A compact C/EBP–AP-1–TWIST1 switch drives the plastic, perivascular mesenchymal cell state of IDH-wildtype glioblastoma*

Analysis code accompanying the manuscript by **Chunying Luo¹, Qingyun Zhang², Yaoyi Xie³, Chuanyu Li⁴˒⁵, Huangde Fu¹˒\*** (manuscript in preparation).

The study is **purely computational (dry-lab)**: it re-analyses public single-cell, spatial-surrogate and bulk
clinical cohorts of glioma. No new experimental data were generated. All reported associations are
exploratory/associational and are stated as such in the manuscript.

---

## ⚠️ Status: partial release — read this first

This repository contains **the code that exists and runs today**, not the complete analysis of the paper.
The following analyses were run in a separate workspace that is **not** included here (see
[What is *not* included](#what-is-not-included)):

- de-novo **pySCENIC** regulon inference per cohort (Figure 3b, 3d–f; Supplementary Figure S2e),
- the **IvyGAP** anatomic-region niche analysis (Figure 5a–c),
- **NicheNet** ligand–target modelling and receptor expression (Figure 5d),
- **CGGA-325/693** MES-score Cox models and the pooled meta-analysis (Figure 5f; Supplementary Figure S2a, c),
- per-cohort trajectory/PAGA panels (Figure 4b, c),
- the figure-composition scripts (`build_figures.py`) and `nichenet_ligand_activity.R` referenced in the
  manuscript's Data-and-code statement.

Consequently, **Figure 1 and Figure 2 are reproducible from this repository**, together with the
prior-based (DoRothEA) part of Figure 3 and the TCGA part of the prognostic analysis; the remaining panels
are **not** independently reproducible from this release. We state this explicitly rather than imply full
reproducibility. **Do not describe this repository as the complete analysis code of the paper.**

---

## What the included code does

| # | Module | Scripts | Manuscript output |
|---|--------|---------|-------------------|
| 1 | GSE103224 discovery pipeline: QC → Harmony integration → inferCNV → Neftel four-state scoring → diffusion pseudotime → DoRothEA TF activity | `00_…`, `01_…`–`05_…` | Figure 1; Figure 3a, c; core MES-transition program |
| 2 | Cross-cohort malignant-state composition (GSE131928, GSE182109, GSE70630) and H1 (primary vs recurrence) | `mc_cohort_state.py`, `mc_gse182109_h1.py` | Figure 2a–f |
| 3 | Cross-cohort TF-activity consensus from prior-based regulons (DoRothEA/decoupleR) | `mc_h2_regulon.py`, `mc_regulon_detail.py` | Figure 3a, c |
| 4 | TCGA-GBM/LGG bulk data acquisition, MES-transition score, Kaplan–Meier and censored Cox models | `download_tcga.py`, `fetch_tcga_*.py`, `parse_idh.py`, `tcga_*.py` | Figure 5f (TCGA panel); Supplementary Figure S2b |
| 5 | Accession verification | `verify_accessions.py` | `tables/datasets.tsv` |

Seeds: **42** (Python) throughout. The statistical unit is the **patient** (patient-level means/proportions),
except for per-cell regulon/trajectory correlations, which are labelled as such.

## Repository layout

```
.
├── README.md               this file
├── DATA_MANIFEST.md        every dataset: accession, source, expected local path, size, consumer script
├── .gitignore
├── env/
│   ├── requirements.txt    pip pins — versions as used (September 2026)
│   └── environment.yml     conda route (Python 3.12)
├── scripts/                24 scripts (pipeline, multi-cohort, TCGA, accessions)
└── tables/
    ├── datasets.tsv           accession-level metadata (verified with verify_accessions.py)
    └── neftel_meta_modules.tsv Neftel 2019 meta-modules (MES1/MES2/AC/OPC/NPC1/NPC2 + G1/S, G2/M)
```

`tables/` holds **inputs/metadata**, not result tables: the pipeline reads
`tables/neftel_meta_modules.tsv` for all state scoring, and the manuscript points readers to
`tables/datasets.tsv` for accessions. **Result tables and intermediate objects are not in this repository.**

## Data

All datasets are public. Accessions, exact source URLs, expected local file layouts and sizes are in
**[DATA_MANIFEST.md](DATA_MANIFEST.md)**. Summary:

| Accession | Cohort | Role in the paper |
|-----------|--------|-------------------|
| GSE103224 (Yuan et al. 2018) | 8 IDH-wt GBM patients, per-sample UMI matrices | discovery cohort (Figures 1) |
| GSE131928 (Neftel et al. 2019) | Smart-seq2 + 10x IDH-wt GBM | cross-cohort composition/TF consensus (Figure 2, 3a) |
| GSE182109 (Abdelfattah et al. 2022) | 44 samples, newly diagnosed vs recurrent | H1 recurrence test (Figure 2d, e) |
| GSE70630 (Venteicher et al. 2017) | IDH-mut astrocytoma | IDH-mutant biological control (Figure 2) |
| TCGA-GBM + TCGA-LGG (GDC) | bulk RNA-seq (STAR counts), clinical, IDH | prognostic analysis (Figure 5f) |

Data are **not** redistributed here; download them from the sources listed in the manifest and respect
GEO/GDC/CGGA/IvyGAP terms of use.

## Environment

Python 3.12 with scanpy, harmonypy, infercnvpy, decoupleR, lifelines, statsmodels
(see `env/requirements.txt` for the versions **as used**). The authors ran this code in a Python
virtualenv; the file is a version record, not a frozen lockfile. Conda users can start from
`env/environment.yml`.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r env/requirements.txt
```

## Quick start

The scripts locate the repository root and the data directory from their own location, and both can be
overridden:

| Variable | Meaning | Default |
|----------|---------|---------|
| `GLIOMA_PROJ` | repository root (where `results/`, `logs/`, `tables/` are written) | parent of `scripts/` |
| `GLIOMA_DATA` | directory holding the downloaded public data (see manifest) | a `Datasets/` directory next to the repo, else `<repo>/data` |

```bash
export GLIOMA_DATA=/path/to/Datasets          # required unless Datasets/ sits next to the repo

# 1) discovery cohort (GSE103224): steps 0-5 in order
bash scripts/run_all.sh

# 2) cross-cohort composition + H1 (needs GSE131928/GSE182109/GSE70630 + step 1-3 outputs)
python3 scripts/mc_gse182109_h1.py
python3 scripts/mc_cohort_state.py

# 3) cross-cohort TF-activity consensus (needs step 5 outputs)
python3 scripts/mc_h2_regulon.py && python3 scripts/mc_regulon_detail.py

# 4) TCGA bulk prognostic analysis (GDC download; parallelisable)
python3 scripts/download_tcga.py            # STAR-count expression
python3 scripts/fetch_tcga_clinical_xml.py  # per-case BCR clinical XML -> censored OS
python3 scripts/fetch_tcga_idh.py           # masked MAFs -> IDH status
python3 scripts/tcga_mes_score.py
python3 scripts/tcga_survival_final.py
```

Outputs go to `results/<step>/`, logs to `logs/`.

## Script reference

| Script | Purpose | Key inputs (under `GLIOMA_DATA`) | Outputs |
|--------|---------|----------------------------------|---------|
| `00_build_gse103224_raw.py` | concatenate the 8 GSE103224 samples into one raw AnnData (idempotent; `--force` to rebuild) | `GSE103224/GSM*_PJ0*.filtered.matrix.txt.gz` | `processed/gse103224_raw.h5ad` |
| `01_qc.py` | protein-coding/lncRNA filtering, mito/UMI/gene thresholds, scrublet doublets | `processed/gse103224_raw.h5ad` | `processed/gse103224_qc.h5ad`, `results/step01_qc/` |
| `02_integration_infercnv.py` | Harmony integration (batch = sample), Leiden clustering, marker-based malignant annotation, inferCNV support | `processed/gse103224_qc.h5ad`, `raw/gencode/gencode.v44.basic.annotation.gtf.gz` | `processed/gse103224_integrated.h5ad`, `results/step02_integration/` |
| `03_state_scoring_shift.py` | Neftel meta-module scoring, four-state assignment, per-sample composition, plasticity entropy, MES1/MES2 split | `processed/gse103224_qc.h5ad`, `…_integrated.h5ad`, `tables/neftel_meta_modules.tsv` | `processed/gse103224_states.h5ad`, `results/step03_state/` |
| `04_trajectory.py` | diffusion pseudotime rooted at the NPC/OPC pole; MES-transition program | `processed/gse103224_states.h5ad` | `processed/gse103224_traj.h5ad`, `results/step04_traj/` (incl. `mes_transition_program_core.tsv`, 131 genes) |
| `05_regulon.py` | DoRothEA prior TF activity (decoupleR MLM) vs MES score | `processed/gse103224_states.h5ad` | `processed/gse103224_regulon.h5ad`, `results/step05_regulon/` |
| `lib.py` | shared paths/IO/loading helpers | — | — |
| `genesets.py` | Neftel module loading + four-state mapping | `tables/neftel_meta_modules.tsv` | — |
| `verify_accessions.py` | verify accessions against NCBI E-utilities | — | `tables/datasets.tsv` |
| `mc_gse182109_h1.py` | H1: MES-like proportion newly diagnosed vs recurrent (patient level, Mann–Whitney) | `cohorts/gse182109/*` | `results/multicohort/gse182109_*` |
| `mc_cohort_state.py` | per-patient four-state composition across cohorts; IDH-wt vs IDH-mut test; per-cohort mean MES figure | `GSE131928/GSM382867{2,3}_*.tsv.gz`, `GSE70630_OG_processed_data_v2.txt.gz` | `results/multicohort/composition_all_cohorts.tsv`, `idh_wt_vs_mut.tsv`, … |
| `mc_h2_regulon.py` | prior-based TF-activity vs MES score in each cohort, BH-FDR | same + `results/step05_regulon/` | `results/multicohort/h2_tf_*.tsv` |
| `mc_regulon_detail.py` | candidate-TF detail plots (Neftel vs Venteicher) | `results/multicohort/h2_tf_*.tsv` | `results/multicohort/regulon_detail.*` |
| `download_tcga.py` / `redownload_tcga_missing.py` | GDC expression download (manifest + STAR counts) | GDC API | `raw/tcga/*.tsv`, `raw/tcga/manifest.tsv` |
| `fetch_tcga_clinical.py` / `download_tcga_clinical.py` | case-level clinical table via the GDC case API | GDC API | `raw/tcga/clinical/` |
| `fetch_tcga_clinical_xml.py` | per-case BCR clinical XML → censored survival table | GDC API | `raw/tcga/clinical_xml/clinical_xml_survival.tsv` |
| `fetch_tcga_idh.py` / `parse_idh.py` | masked-MAF download and IDH1/2 status parsing | GDC API | `raw/tcga/mutation/idh_status.tsv` |
| `tcga_mes_score.py` | MES-transition score (z-scored mean of MES1+MES2 genes, log1p-TPM) + clinical join | `raw/tcga/*.tsv`, clinical, IDH | `results/step07_clinical/tcga_mes_clinical.tsv` |
| `tcga_survival.py` / `tcga_survival_final.py` | IDH-stratified Kaplan–Meier, log-rank, multivariable censored Cox (age, IDH, grade) | `results/step07_clinical/tcga_mes_clinical.tsv`, `raw/tcga/clinical_xml/clinical_xml_survival.tsv` | `results/step07_clinical/tcga_survival.tsv`, `tcga_cox.tsv` |
| `tcga_km_fig.py` | standalone KM figure | `results/step07_clinical/tcga_survival.tsv` | `tcga_km.png` |

## What is *not* included

1. **de-novo regulon inference** — pySCENIC/GRNBoost2 + cisTarget motif pruning + AUCell per cohort.
   Required for Figure 3b, 3d–f and Supplementary Figure S2e.
2. **IvyGAP niche analysis** — per-region MES/hypoxia scoring (GSE107559), Kruskal–Wallis + mixed model.
   Required for Figure 5a–c.
3. **NicheNet** ligand–target modelling and receptor expression on malignant cells; the manuscript cites
   `nichenet_ligand_activity.R`. Required for Figure 5d.
4. **CGGA validation** — MES-score Cox in CGGA-325/CGGA-693 and the fixed/random-effects meta-analysis.
   Required for Figure 5f (CGGA panel) and Supplementary Figure S2a, c.
5. **Per-cohort trajectory/PAGA panels** (Figure 4b, c) and the **figure-composition code**
   (`build_figures.py`) cited in the manuscript's Data-and-code statement.
6. **The Visium spatial extension** — the Ravi et al. 2022 spatial objects were not obtainable at the time
   of analysis (the manuscript notes the spatial extension as pending).
7. **Result tables, intermediate objects and the manuscript itself** — this is a code-only release.

Until items 1–5 are added, treat this repository as the *partial* code companion of the paper.

## Caveats

- Everything here is **in silico, associational**; regulon inference is correlative and does not establish
  binding or causality.
- Cross-cohort single-cell comparisons use cohort-appropriate normalisation; IDH-mutant series are
  represented by Venteicher (GSE70630) among the released scripts.
- The TCGA intermediate tables used for the analysis on the authors' machine were regenerated by the
  scripts above; if you re-run the TCGA arm, re-download from GDC rather than assuming cached files.
- Python version and package versions are given as used; newer scanpy/lifelines releases may change
  defaults.

## Citation

If you use this code, please cite the manuscript (in preparation) and the primary data sources listed in
[DATA_MANIFEST.md](DATA_MANIFEST.md).

## Licence

MIT — see [LICENSE](LICENSE).

## Contact

Chuanyu Li (854309682@qq.com) · Huangde Fu (qt000172@sr.gxmu.edu.cn)
