# Data manifest

Every dataset used by the released scripts: what it is, where it comes from, where the scripts expect it,
how big it is, and which script consumes it. **No data are redistributed in this repository** — download
from the sources below and observe each resource's terms of use.

Set the data root once:

```bash
export GLIOMA_DATA=/path/to/Datasets      # see "Expected layout" below
```

## Expected layout

```
$GLIOMA_DATA/
├── GSE103224/                          # GSE103224 per-sample UMI matrices (genes x cells)
│   └── GSM2758471_PJ016.filtered.matrix.txt.gz …  (8 samples, PJ016–PJ048)
├── GSE131928/                          # Neftel 2019 processed TPM (Smart-seq2 + 10x)
│   ├── GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz
│   └── GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz
├── cohorts/
│   └── gse182109/                      # GSE182109 Cell Ranger per-sample matrices (44 samples)
│       └── *_matrix.mtx.gz (+ barcodes/features)
├── GSE70630_OG_processed_data_v2.txt.gz       # Venteicher 2017, log gene x cell
├── GSE89567_IDH_A_processed_data.txt.gz       # Tirosh 2016 IDH-mut astrocytoma (not used by the released scripts)
├── IvyGAP/                             # GSE107559 anatomic-region FPKM (not used by the released scripts)
├── CGGA/                               # CGGA-325 / CGGA-693 zips (not used by the released scripts)
├── raw/
│   ├── MANIFEST.tsv                    # provenance notes for the raw downloads
│   ├── gencode/gencode.v44.basic.annotation.gtf.gz
│   └── tcga/                           # GDC downloads: STAR counts, clinical, mutation (see below)
└── processed/                          # intermediate AnnData objects (regenerable; ~13 GB)
```

`processed/` is a **working directory**: `00_build_gse103224_raw.py` → `gse103224_raw.h5ad`,
`01_qc.py` → `gse103224_qc.h5ad`, … `05_regulon.py` → `gse103224_regulon.h5ad`. It can be deleted and
rebuilt from the raw downloads.

## Single-cell and spatial-surrogate cohorts

| Accession | Study | Role | Download | Expected path | Size on disk | Consumed by |
|-----------|-------|------|----------|---------------|--------------|-------------|
| **GSE103224** | Yuan et al. 2018, *Genome Med* (PMID 30041684) | discovery cohort, 8 IDH-wt GBM patients (PJ016–PJ048) | GEO FTP: `GSE103224_RAW.tar` (also per-sample `GSM2758471…` matrices) | `GSE103224/` (extracted `*.filtered.matrix.txt.gz`) | 5.5 GB (archive 52 MB) | `00_build_gse103224_raw.py` → `01–05` |
| **GSE131928** | Neftel et al. 2019, *Cell* (PMID 31327527) | cross-cohort composition + prior-based TF consensus (IDH-wt GBM) | GEO: `GSE131928_RAW.tar`, plus processed TPM `GSM3828672`/`GSM3828673` | `GSE131928/GSM382867{2,3}_*.tsv.gz` | 639 MB (archive 638 MB) | `mc_cohort_state.py`, `mc_h2_regulon.py` |
| **GSE182109** | Abdelfattah et al. 2022, *Nat Commun* (PMID 35140215) | H1: newly diagnosed vs recurrent GBM (44 samples) | GEO: `GSE182109_RAW.tar` (Cell Ranger outputs) | `cohorts/gse182109/` | 2.3 GB | `mc_gse182109_h1.py`, `mc_h2_regulon.py` |
| **GSE70630** | Venteicher et al. 2017, *Science* (PMID 28360267) | IDH-mutant astrocytoma control | GEO: `GSE70630_OG_processed_data_v2.txt.gz` | `GSE70630_OG_processed_data_v2.txt.gz` | 74 MB | `mc_cohort_state.py`, `mc_h2_regulon.py` |
| GSE89567 | Tirosh et al. 2016, *Nature* (PMID 27806376) | IDH-mutant astrocytoma (IDH-wt vs IDH-mut comparison) | GEO: `GSE89567_IDH_A_processed_data.txt.gz` | `GSE89567_IDH_A_processed_data.txt.gz` | 119 MB | *no released script* (see README → not included) |
| GSE107559 (IvyGAP) | Puchalski et al. 2018, *Science* (PMID 29748285) | anatomic-region spatial surrogate (LE/IT/CT/MVP/HBV/PNZ/PAN) | GlioVis/IvyGAP portal → GEO GSE107559 | `IvyGAP/` | 38 MB | *no released script* (see README → not included) |
| CGGA mRNAseq-325 / mRNAseq-693 | Zhao et al. 2021, *Neuro Oncol* (PMID 33662628) | external validation of the MES-transition score | http://www.cgga.org.cn (registration required) | `CGGA/` (zips) | 85 MB | *no released script* (see README → not included) |

## Reference files

| File | Purpose | Source | Expected path | Consumed by |
|------|---------|--------|---------------|-------------|
| GENCODE v44 basic annotation (GTF) | gene types/coordinates for expression filtering and inferCNV | https://www.gencodegenes.org/human/release_44.html | `raw/gencode/gencode.v44.basic.annotation.gtf.gz` | `02_integration_infercnv.py` |
| Neftel 2019 meta-modules | MES1/MES2/AC/OPC/NPC1/NPC2 + G1/S, G2/M gene lists — **shipped in this repo** | reproduction of the paper's supplementary tables (PMID 31327527); provenance recorded in `raw/MANIFEST.tsv` | `tables/neftel_meta_modules.tsv` | all state scoring (`03_…`, `mc_*`) |

## TCGA-GBM / TCGA-LGG (GDC)

Downloaded through the GDC REST API by the scripts themselves (no registration needed).

| Layer | Script | Expected path | Notes |
|-------|--------|---------------|-------|
| STAR-count expression (926 files: 391 GBM + 534 LGG, plus manifest) | `download_tcga.py`, `redownload_tcga_missing.py` | `raw/tcga/*.rna_seq.augmented_star_gene_counts.tsv`, `raw/tcga/manifest.tsv` | ~3.7 GB; the download is resumable (re-run `redownload_tcga_missing.py`) |
| Case-level clinical table | `fetch_tcga_clinical.py`, `download_tcga_clinical.py` | `raw/tcga/clinical/` | age, grade, histology, vital status |
| Per-case BCR clinical XML → **censored** survival | `fetch_tcga_clinical_xml.py` | `raw/tcga/clinical_xml/clinical_xml_survival.tsv` | required for the Cox models; `days_to_last_followup` for censored cases (the `demographic` endpoint has none) |
| IDH1/IDH2 status from masked somatic MAFs | `fetch_tcga_idh.py` → `parse_idh.py` | `raw/tcga/mutation/idh_status.tsv` | MAF parts are gzip; `parse_idh.py` deletes the parts after parsing |

**Note.** On the authors' machine the intermediate TCGA tables above were later regenerated by re-running
these scripts; the cached clinical/XML/IDH tables are deliberately not part of this repository, so the TCGA
arm is re-created from GDC by the scripts in the order listed. The released result tables
(`results/step07_clinical/`) were produced from GDC data downloaded in September 2026 — re-running against
today's GDC snapshot may shift individual numbers slightly (cohort n, follow-up updates).

## Get the data (example)

```bash
# GEO cohorts (adjust paths to your data root)
cd "$GLIOMA_DATA"
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103224/suppl/GSE103224_RAW.tar
tar -xf GSE103224_RAW.tar -C GSE103224/
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_RAW.tar
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3828nnn/GSM3828672/suppl/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3828nnn/GSM3828673/suppl/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz"
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182109/suppl/GSE182109_RAW.tar
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE706nnn/GSE70630/suppl/GSE70630_OG_processed_data_v2.txt.gz
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89567/suppl/GSE89567_IDH_A_processed_data.txt.gz

# GENCODE v44 (basic annotation)
wget -c -P raw/gencode https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz

# TCGA via GDC — use the scripts (they page the API and write the manifest)
python3 scripts/download_tcga.py
python3 scripts/fetch_tcga_clinical_xml.py
python3 scripts/fetch_tcga_idh.py
```

GSE103224 ships per-sample matrices in several formats; the pipeline reads the `*.filtered.matrix.txt.gz`
files (column 0 = Ensembl gene id, column 1 = HGNC symbol, then UMI counts). GSE182109 requires the Cell
Ranger `*_matrix.mtx.gz` triplets to be extracted into `cohorts/gse182109/`.

## Disk space

Roughly **25 GB** for the full working set (of which `processed/` ≈ 13 GB is regenerable), plus the
optional datasets that the released scripts do not consume.
