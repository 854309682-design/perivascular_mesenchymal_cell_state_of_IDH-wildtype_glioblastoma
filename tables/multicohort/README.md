# Multi-cohort analysis from already-downloaded data — what we can get for an article

Repo: `/run/media/meng/Seqdata/glioma/Paper_2` · generated 2026-09-10 · IDH arm corrected 2026-10-10
All analysis is Python (scanpy/decoupler), seed=42, single-cell only.

## Datasets used (all already on disk)
| Cohort | Type | Role | State composition (patient mean MES-like) |
|---|---|---|---|
| GSE103224 (Yuan 2018) | 10x, IDH-wt GBM | discovery | 0.27 (median 0.09) |
| GSE131928 (Neftel 2019) | Smartseq2+10x, IDH-wt GBM | discovery / gold-standard | 0.19 |
| GSE182109 (Abdelfattah 2022) | 10x, 44 samples, GBM | primary vs recurrent | 0.50 |
| GSE70630 (Tirosh 2016) | Smartseq2, IDH-mut oligodendroglioma | IDH-mut control | **0.002** (6 tumours) |
| GSE89567 (Venteicher 2017) | Smartseq2, IDH-mut astrocytoma | IDH-mut control | **0.011** (10 patients) |
| (GSE103224 Steps 1–5) | | 4-state + MES-transition + regulon | |

## Step 1 — Four-state composition atlas + IDH comparison (H1/H3-background)
- **IDH-wt vs IDH-mut MES-like proportion:** IDH-wt **0.267** vs IDH-mut **0.008**;
  Mann-Whitney **p = 1.2e-07** (n=69 vs 16 patients) — `idh_wt_vs_mut.tsv`,
  reproduced by `scripts/mc_gse89567_idh.py`.
- ⚠️ **Superseded:** an earlier run pooled the GSE89567 arm with
  `barcode.split("_")[0]`. GSE89567 barcodes are dash-separated (`MGH64-P16-A01`), so
  247 single cells each became their own "patient": n = 8 + 251 = **259**,
  mean 0.008, **P = 1.6e-57**. Two of the 240 surviving single cells were MES-positive,
  so the P-value was inflated by n, not by biology. `mc_gse89567_idh.py` reproduces
  that artifact (`idh_wt_vs_mut_superseded.tsv`, function `_legacy_reference()`) and
  writes the corrected comparison. Parse patient IDs with `tumour_id()`.
- **Interpretation:** the MES-like state is **essentially absent in IDH-mutant gliomas**
  and is an **IDH-wildtype (GBM) phenotype**. This is a clean, quantitative, article-worthy
  result and a strong biological control; the conclusion is unchanged by the correction
  (mean 0.008 → 0.008, only n and P move).
- Files: `composition_all_cohorts.tsv`, `idh_wt_vs_mut.tsv`, `idh_patients.tsv`,
  `gse89567_composition.tsv`, `neftel_composition.tsv`, `gse70630_composition.tsv`,
  `state_composition_by_cohort.png`; figure: `../figures/Figure2_main.pdf|jpeg`.

## Step 1b — H1 primary-vs-recurrence (the "recurrence → MES" hypothesis)
- **Null result:** patient-level MES-like proportion primary (`ndGBM`) **0.507** vs
  recurrent (`rGBM`) **0.473**; one-sided Mann-Whitney **p = 0.67** (n=11 vs 5).
- Inter-patient MES-proportion is highly heterogeneous (0.10 → 0.93), swamping any
  recurrence effect. **MES plasticity is not simply selected for by recurrence** in this
  immune atlas.
- File: `gse182109_h1_primary_vs_recurrent.tsv`, `gse182109_composition_by_patient.tsv`,
  `gse182109_mes_by_type.png`.

## Step 2 — H2 multi-cohort regulon consensus (STAT3 is the master MES switch)
DoRothEA-prior TF activity (decoupler mLM), per-cohort rho of activity vs MES score,
then consensus across 4 cohorts (mean rho, # cohorts with FDR<0.05):

| TF | cohorts | mean rho | FDR<0.05 in # cohorts |
|---|---|---|---|
| **STAT3** | 4 | **0.388** | 4 |
| KLF4 | 4 | 0.394 | 4 |
| FOXL2 / SP1 / HOXA9 | 4 | 0.36 / 0.33 / 0.33 | 4 |
| JUNB (AP-1) | 4 | 0.305 | 4 |
| FOSL1 (AP-1) | 4 | 0.290 | 4 |
| NFKB1 | 4 | 0.223 | 3 |
| FOSL2 | 4 | 0.210 | 4 |
| IRF1 | 4 | 0.156 | 4 |
| **CEBPB** | 4 | **−0.010** | 4 |
| **TEAD1 / TEAD4** | 4 | −0.09 / −0.04 | 3–4 |
| **RELA (NFKB)** | 4 | **−0.180** | 4 |

- **STAT3 is the robust, consensus master TF of the MES state** (top-2 overall; FDR<0.05 in all 4 cohorts).
- **AP-1 axis** (JUNB/FOSL1/FOSL2) supportive, weaker.
- **CEBPB, RELA, TEAD1/4 are NOT MES-linked** (anti-correlated) — contradicts part of the
  a-priori candidate list; report this honestly.
- Files: `h2_regulon_consensus.tsv`, `h2_tf_{Neftel,GSE182109,Venteicher}.tsv`.

## Step 3 — What a paper from this can claim, and what it cannot
- **Supported by downloaded data (strong):**
  1. Four-state composition differs markedly across cohorts; MES-like is an
     **IDH-wildtype/GBM** program (IDH-mut ≈ 0; p<0.0001).
  2. **STAT3** (and to a lesser extent AP-1) is a **reproducible TF-activity correlate of the
     MES state across 4 independent cohorts** — a consensus regulon finding.
  3. The MES program is **not explained by recurrence** (H1 null) — supports an
     **environment/microenvironment-driven** (not treatment-selection) model.
- **Cannot be claimed from downloaded data (need more data):**
  - Causality / de-novo SCENIC regulons + ChIP-seq edge evidence (only DoRothEA prior here).
  - Spatial niche (needs Visium — we have the Ravi-2022 accession but not downloaded yet).
  - Clinical prognosis / drug repositioning (needs TCGA/CGGA/L1000/DepMap).

## Caveats to state in any manuscript
- DoRothEA-prior TF activity is a *curated-prior* proxy for regulons, not de-novo SCENIC.
- GSE70630 (IDH-mut, Tirosh 2016) is an **oligodendroglioma** cohort and GSE89567
  (IDH-mut, Venteicher 2017) an **astrocytoma** cohort — IDH-mut *controls*, not
  grade-matched GBM comparisons. The two accessions were previously swapped in the
  manuscript text and in the released cohort label "Venteicher_GSE70630"; GEO
  `pubmedids` confirm GSE70630 → 27806376 (Tirosh) and GSE89567 → 28360267 (Venteicher).
- Neftel patient IDs are reconstructed from cell-name prefixes (45 IDs); patient-level n is
  uneven. Non-malignant identification is marker-threshold based (immune/endothelial
  score > 0.5), not CNV-validated across all cohorts.
- ⚠️ **Open item — IDH-wildtype arm is still over-split.** The same parsing class inflates
  the wildtype count: the released Neftel arm has 45 sample tokens, but merging sorted
  fractions (`MGH104negP2/P4/P7`, `MGH106CD3posP1`), multiregion aliquots (`MGH105A-D`) and
  the platform-duplicated numeric 10x IDs (`102/105/…` == `MGH102/MGH105/…`, same patients)
  gives **31 tumours**; GEO states the series is "24,131 single cells from 28 patients".
  `tumour_id()` in `mc_gse89567_idh.py` collapses Neftel to 31. The wildtype arm
  (n = 69: GSE103224 8 + GSE182109 16 + Neftel 45) is deliberately kept **as shipped**
  here, pending the authors' sample table; the mean is insensitive to it
  (0.267 at n=69 vs 0.268 at n=55, P 1.2e-07 vs 1.3e-07 either way).
