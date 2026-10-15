# Release notes — v1.0.2 (ready to paste into GitHub)

**Who has to publish this:** the repository owner. As of 2026-10-15 the account
`LingzhangMeng` has only **read** access (`permissions: {"pull": true, "push": false}`) on
`854309682-design/perivascular_mesenchymal_cell_state_of_IDH-wildtype_glioblastoma`, so the
release API returns 404. The pushes that shipped v1.0.1/v1.0.2 went through the repository's
**SSH deploy key** (`.ssh/glioma_paper2`), which has write access but cannot call the API.
The owner (contact in the README: Chuanyu Li, 854309682@qq.com) — or anyone given write
access — should publish it.

**One-click link** (the tag `v1.0.2` already exists on `main`):

```
https://github.com/854309682-design/perivascular_mesenchymal_cell_state_of_IDH-wildtype_glioblastoma/releases/new?tag=v1.0.2
```

Set **Release title** and paste the body below. Publishing is what triggers Zenodo: Zenodo's
GitHub integration mints a **new version DOI** for a published release, and the concept DOI
`10.5281/zenodo.23262404` (already cited in the manuscript) then resolves to it.

---

## Release title

```
v1.0.2 — revision pipeline, sensitivity analysis and text-pass corrections
```

## Release body

Code and result tables behind the revised manuscript (review-response revision).

**What this release contains**

- the full revision pipeline: de-novo regulon inference (pySCENIC: GRNBoost2 → cisTarget →
  AUCell) and its reconciliation with the prior-based route; the re-implemented DoRothEA prior
  route; the IvyGAP region re-analysis (atlas samples only, composition-adjusted, NNLS
  deconvolution); TCGA/CGGA survival from raw data with a random-effects meta-analysis; base-R
  NicheNet ligand activity and receptor fractions; the per-cohort trajectory + PAGA code behind
  Figure 4; the Figure 2/S1/S2 builders; the submission-archive packager
- `rev27_consensus_sensitivity.py`: the prior-route consensus computed on the IDH-wildtype
  cohorts only and leave-one-cohort-out, with the rank of every candidate factor per variant and
  the per-cohort MES-score spread
- every result table quoted in the manuscript (`tables/revision/`, `tables/multicohort/`);
  per-cell intermediates are excluded because the scripts regenerate them

**Two methodological corrections** (stated in the manuscript)

- the de-novo arm is capped at 500 cells per state / 2,000 cells per cohort (GRNBoost2 cost is
  linear in cells × genes), and
- the GSE131928 cohort is scored on its **Smart-seq2 block only**, because its two public
  matrices are delivered on different scales (log2(TPM+1) versus raw TPM); the concatenated
  matrix made that cohort's MES score and pseudotime axis predominantly a *platform* axis.

**Supersedes:** v1.0.1 (adds `rev27`, the second-round text-pass corrections and the
concept-DOI citation).
