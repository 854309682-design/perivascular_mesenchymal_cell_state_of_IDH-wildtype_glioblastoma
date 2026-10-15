# Multi-cohort analysis — provenance of the per-patient tables

Repo: `/run/media/meng/Seqdata/glioma/Paper_2` · generated 2026-09-10 · IDH arm corrected 2026-10-10 · **Neftel platform correction 2026-10-14**
All analysis is Python (scanpy/decoupler/pySCENIC), seed = 42, single-cell + bulk-region only.

> **Read this first — which table is canonical.**
>
> | table | status |
> |---|---|
> | `idh_patients.tsv` | **canonical** patient-level table: five cohorts, corrected patient definition, `genotype` + `MES_prop`. Every IDH number and Figure 2 in the manuscript comes from here. |
> | `composition_all_cohorts.tsv` | **intermediate.** Its four-state columns (AC/MES/NPC/OPC proportions per patient) are the raw composition; its **IDH-mutant block is deliberately superseded** (it carries the pre-correction GSE70630 tokens and no GSE89567). Only its **wild-type rows** are consumed, by `scripts/mc_gse89567_idh.py`. |
> | `gse70630_composition.tsv`, `gse89567_composition.tsv` | per-cohort composition; the GSE89567 file is the corrected (10-patient) version, the GSE70630 file is the **pre-correction** one. |
> | `*_state_scores.tsv` | per-cell tables (patient, cohort, state, score_MES). Large; not committed to the public repo. |
> | `composition_all_cohorts.tsv`, `idh_wt_vs_mut.tsv`, `idh_patients.tsv` | **rebuilt 2026-10-14 by `scripts/rev23_neftel_ss2.py`** for the GSE131928 rows; the previous versions are in `superseded_2026-10-14/`. |

## Datasets used (all already on disk)

| Cohort | Platform | Role | Patients | Malignant cells | Mean MES-like |
|---|---|---|---|---|---|
| GSE103224 (Yuan 2018) | 10x, IDH-wt GBM | discovery | 8 | 20,643 | 0.265 |
| GSE131928 (Neftel 2019) | **Smart-seq2 only** (see §0) | discovery / gold standard | **35** | **6,701** | **0.216** |
| GSE182109 (Abdelfattah 2022) | 10x, primary + recurrent | recurrence arm (H1) | 16 | 142,128 | 0.496 |
| GSE70630 (Tirosh 2016) | Smart-seq2, IDH-mut oligodendroglioma | IDH-mut control | 6 | 4,082 | **0.002** |
| GSE89567 (Venteicher 2017) | Smart-seq2, IDH-mut astrocytoma | IDH-mut control | 10 | 5,187 | **0.011** |

## §0 — GSE131928 platform correction (2026-10-14)

The two public matrices for this cohort are delivered on **different scales**: Smart-seq2 is
log2(TPM+1) (nonzero values 0.02–10) while the 10x file is **raw TPM** (86–691). The original
analysis concatenated them unexamined — `mc_cohort_state.py` called both "log-scaled TPM
files" — which made this cohort's per-cell MES score and any pseudotime axis **predominantly a
platform axis**:

| variant | malignant | MES share | ρ(platform, MES) | ρ(platform, dpt) | ρ(dpt, MES) |
|---|---|---|---|---|---|
| A as shipped (mixed scale) | 15,809 | 0.243 | **+0.815** | **−0.898** | −0.818 |
| B both platforms log2(TPM+1) | 15,700 | 0.297 | +0.425 | −0.645 | −0.209 |
| **C Smart-seq2 only (used)** | **6,701** | **0.173** | n/a | n/a | **+0.039** |

Variant C is used: `rev23_neftel_ss2.py` re-scores the cohort on its Smart-seq2 block (the
platform that carries that study's per-patient design) and rewrites `neftel_composition.tsv`,
`Neftel_GSE131928_state_scores.tsv`, the GSE131928 rows of `composition_all_cohorts.tsv` and
`idh_patients.tsv`, and `idh_wt_vs_mut.tsv`. Consequences: **35 patients** (was 45), mean MES
0.216 (was 0.185), and the cohort **no longer supports the multi-cohort pseudotime claim**
(ρ = 0.04 on its own platform). `rev01_export_scenic.py` and `rev02_mlm_reimplement.py` default
to the same restriction; set `NEFTEL_SS2_ONLY=0` to reproduce the superseded behaviour. The
naive concatenated matrix cannot even be passed to scanpy's variance/expm1 paths (overflow →
inf → `ValueError`), which is how the defect surfaced.

## §1 — Four-state composition + IDH comparison

- **IDH-wt vs IDH-mut MES-like proportion:** IDH-wt **0.299** (n = 59) vs IDH-mut **0.008**
  (n = 16); Mann-Whitney **P = 6.6 × 10⁻⁸** — `idh_wt_vs_mut.tsv`. A **cohort-mean-centred**
  test (which removes cohort-level shifts) is not significant (**P = 0.25**), because genotype
  is fully confounded with cohort, histology, grade and platform; both are reported in the
  manuscript and in Supplementary Fig. S2e.
- Per-cohort means (IDH-wt): GSE103224 0.265, GSE182109 0.496, GSE131928 0.216;
  (IDH-mut): GSE70630 0.002, GSE89567 0.011 → `../../results/revision/rev15_idh_stratified.tsv`.
- ⚠️ **Superseded:** an earlier run pooled the GSE89567 arm with `barcode.split("_")[0]`.
  GSE89567 barcodes are dash-separated (`MGH64-P16-A01`), so 247 single cells each became
  their own "patient": n = 8 + 251 = **259**, mean 0.008, **P = 1.6e-57**. Two of the 240
  surviving single cells were MES-positive, so the P-value was inflated by n, not by biology.
  `mc_gse89567_idh.py` reproduces that artifact (`idh_wt_vs_mut_superseded.tsv`, function
  `_legacy_reference()`) and writes the corrected comparison. Parse patient IDs with
  `tumour_id()`.
- ⚠️ **Open item — the wild-type arm may still be over-split.** Neftel patient IDs are
  reconstructed from cell-name prefixes; after the platform correction the arm has **35
  Smart-seq2 sample tokens**, and merging sorted fractions (`MGH104negP2/P4/P7`,
  `MGH106CD3posP1`), multiregion aliquots (`MGH105A–D`) and numeric 10x aliases gives ~31
  tumours (GEO describes the series as 28 patients). `tumour_id()` in `mc_gse89567_idh.py`
  performs that collapse. The wild-type arm is deliberately kept **as shipped** pending the
  authors' sample table, because the mean is insensitive to it; if the collapse is applied,
  `idh_patients.tsv`, Figure 2a–c and Table 1 must be rebuilt together.
- **Interpretation:** the MES-like state is **essentially absent in IDH-mutant gliomas** and is
  an **IDH-wildtype (GBM) phenotype** — a clean, quantitative biological control whose
  direction survives every correction (superseded n = 259 → corrected n = 16 → Neftel
  platform-corrected n = 59/16).
- Files: `composition_all_cohorts.tsv`, `idh_wt_vs_mut.tsv`, `idh_patients.tsv`,
  `gse89567_composition.tsv`, `neftel_composition.tsv`, `gse70630_composition.tsv`,
  `state_composition_by_cohort.png`; figure: `../../figures/Figure2_main.pdf|jpeg`
  (built by `scripts/mc_fig2_idh_recurrence.py`, **not** by `rev08_figures.py`, whose Figure-2
  builder is a documented no-op after it clobbered the corrected figure once).

## §1b — H1 primary-vs-recurrence ("recurrence → MES")

- **Null result:** patient-level MES-like proportion primary (`ndGBM`) **0.507** vs recurrent
  (`rGBM`) **0.473**; one-sided Mann-Whitney **P = 0.67** (n = 11 vs 5); two-sided **P = 0.74**,
  mean difference −0.034 (bootstrap 95 % CI −0.354 … +0.285), Cliff's delta −0.127.
- **Power-limited:** with 5 recurrent patients a two-sided test at α = 0.05 has < 30 % power
  for the observed effect size — absence of evidence, not evidence of absence. The manuscript
  says so explicitly.
- Inter-patient MES proportion is highly heterogeneous (0.10 → 0.93), swamping any recurrence
  effect: **MES plasticity is not simply selected for by recurrence.**
- Files: `gse182109_h1_primary_vs_recurrent.tsv`, `gse182109_composition_by_patient.tsv`,
  `gse182109_mes_by_type.png`, `../../results/revision/rev02_recurrence_power.tsv`.

## §2 — H2 multi-cohort regulon consensus (prior route, re-implemented 2026-10-14)

DoRothEA-prior TF activity (multivariate linear model, re-implemented and validated against
the archived decoupler matrix: **median per-TF ρ = 0.998**), per-cohort ρ of activity vs the
MES score, consensus across the four cohorts in which that route was run. The values below are
**current** (`../../results/revision/rev02_tf_consensus_reimplemented.tsv`), i.e. after the
GSE131928 platform correction:

| TF | cohorts | mean ρ | FDR < 0.05 in # cohorts |
|---|---|---|---|
| **STAT3** | 4 | **+0.421** | 4 |
| SP1 | 4 | +0.395 | 4 |
| **FOSL1** (AP-1) | 4 | **+0.369** | 4 |
| HIF1A | 4 | +0.366 | 4 |
| KLF4 | 4 | +0.360 | 4 |
| PPARG / HOXA9 | 4 | +0.338 / +0.335 | 4 |
| ZEB2 | 4 | +0.303 | 4 |
| **FOSL2** (AP-1) | 4 | +0.235 | 4 |
| NFKB1 | 4 | +0.233 | 3 |
| **JUNB** (AP-1) | 4 | +0.210 | 4 |
| **CEBPD** | 4 | +0.124 | 3 |
| TEAD1 | 4 | +0.057 | 4 |
| **CEBPB** | 4 | **−0.122** | 4 |
| **RELA** (NF-κB) | 4 | **−0.167** | 4 |
| **TWIST1** | 4 | **−0.204** | 4 |

- **STAT3 is the robust consensus master correlate** of the MES state (top of the list, and
  FDR < 0.05 in 4/4 cohorts).
- The **AP-1 axis** (FOSL1 > FOSL2 > JUNB) is supportive.
- **CEBPB, RELA and TWIST1 are anti-correlated** in every cohort; **TEAD1 is positive** here
  (it was quoted as negative before the correction) and is therefore no longer placed in the
  switch. Patient-level inference attenuates everything: JUNB pooled +0.21 → per-patient
  median +0.05; STAT3 +0.42 → +0.31; FOSL1 +0.37 → +0.29
  (`../../results/revision/rev02_pooled_vs_patient.tsv`).
- The **prior-independent (pySCENIC) route** is separate: `../../results/revision/rev07_*.tsv`
  — three cohorts, 211/348/470 regulons, consensus headed by CEBPD, BATF3, EPAS1, FOSL2,
  CEBPB; top-30 overlap with the prior route **5/30 (FOSL1, FOSL2, JUNB, NFIC, STAT1), Fisher
  exact P = 0.0075**.
- Files: `h2_regulon_consensus.tsv`, `h2_tf_{Neftel,GSE182109,Venteicher}.tsv` (2026-09-10
  snapshots; the current values are the `rev02_*` tables).

## §3 — What this data can and cannot support

- **Supported (strong):** (1) MES-like is an IDH-wildtype/GBM program (IDH-mut ≈ 0;
  P = 6.6 × 10⁻⁸ pooled, not significant after cohort-mean centring); (2) STAT3 plus a smaller
  AP-1/C-EBP component is a reproducible TF-activity correlate of the MES state across four
  cohorts; (3) the MES program is not explained by recurrence (H1 null, power-limited).
- **Not supported by this data:** causality (regulon inference is correlative), the
  perivascular niche claim (the IvyGAP re-analysis places the signal in the peri-necrotic
  region instead — `../../results/revision/rev04_*`), and a multi-cohort *pseudotime* claim
  (only the discovery cohort supports it — `rev19_*`).
- **Needs more data:** spatial transcriptomics (the Ravi 2022 objects were not obtainable when
  the analysis ran), and ChIP-seq/perturbation for direct binding.

## Caveats to state in any manuscript

- DoRothEA-prior TF activity is a *curated-prior* proxy for regulons, not de-novo SCENIC.
- GSE70630 (Tirosh 2016) is an **oligodendroglioma** cohort and GSE89567 (Venteicher 2017) an
  **astrocytoma** cohort — IDH-mut *controls*, not grade-matched GBM comparisons. The two
  accessions were previously swapped in the manuscript text and in the released cohort label
  "Venteicher_GSE70630"; GEO `pubmedids` confirm GSE70630 → 27806376 (Tirosh) and
  GSE89567 → 28360267 (Venteicher).
- Non-malignant identification is marker-threshold based (immune/endothelial score > 0.5), not
  CNV-validated across all cohorts.
- Platform differs across cohorts (10x vs Smart-seq2), so genotype is confounded with platform;
  this is stated wherever the IDH comparison is made.
