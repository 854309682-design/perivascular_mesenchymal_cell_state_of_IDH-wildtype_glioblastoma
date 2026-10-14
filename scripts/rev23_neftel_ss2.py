#!/usr/bin/env python3
"""rev23 — Neftel (GSE131928) re-scored on its Smart-seq2 block, and the Neftel-dependent
cross-cohort tables rebuilt.

Why
---
The two GEO matrices for this cohort are delivered on **different scales**: the Smart-seq2
file is log2(TPM+1) (nonzero values 0.02-10) while the 10x file is raw TPM (86-691). The
shipped analysis concatenated them unexamined — `release/scripts/mc_cohort_state.py` line 97
calls both "log-scaled TPM files" — which makes this cohort's per-cell MES score and any
pseudotime axis mostly a *platform* axis. Measured on the as-shipped matrices
(`/tmp/p2probe/neftel_scale2.log`, reproduced in logs/rev23.log):

    variant                     malignant   MES share   rho(dpt,MES)   rho(platform,MES)   rho(platform,dpt)
    A as shipped (mixed scale)     15,809       0.243         -0.818              +0.815              -0.898
    B 10x log2(TPM+1)              15,700       0.297         -0.209              +0.425              -0.645
    C Smart-seq2 only               6,701       0.173         +0.039                 n/a                 n/a

Option C is used here: the Smart-seq2 block is internally consistent, carries the
per-patient design of that study, and removes the platform confounding entirely. Under C the
MES-pseudotime correlation of this cohort is ~0, which is reported as such.

What it rewrites (results/multicohort/)
--------------------------------------
    neftel_composition.tsv               per-patient four-state composition
    Neftel_GSE131928_state_scores.tsv    per-cell patient/state/score_MES (malignant cells)
    composition_all_cohorts.tsv          Neftel rows replaced; other cohorts untouched
    idh_patients.tsv                     Neftel rows replaced
    idh_wt_vs_mut.tsv                    recomputed from the updated idh_patients.tsv

The previous versions are first copied to results/multicohort/superseded_2026-10-14/.
Cell scoring uses the same code path as every other cohort
(`rev01_export_scenic.score_states`, identical to the shipped
`mc_cohort_state.score_state` + `exclude_immune`: Neftel meta-modules with ctrl_size = 50,
marker-score exclusion at 0.5, state = argmax of the four pooled scores).

Downstream: re-run `scripts/rev15_idh_stratified_table1.py` (Table 1 + IDH arm),
`scripts/rev02_mlm_reimplement.py` (prior route; its Neftel loader is patched the same way),
`scripts/rev19_figure4_trajectory.py` (Figure 4) and `scripts/rev08_figures.py` (Figures 2/3/5).

Run:
  GLIOMA_PROJ=/run/media/meng/Seqdata/glioma/Paper_2 \
  GLIOMA_DATA=/run/media/meng/Seqdata/glioma/Datasets \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev23_neftel_ss2.py
"""
import os
import sys
import shutil

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES
import rev01_export_scenic as r1

MC = os.path.join(lib.PROJ, "results", "multicohort")
BAK = os.path.join(MC, "superseded_2026-10-14")
COHORT = "Neftel_GSE131928"
FILES = ["neftel_composition.tsv", f"{COHORT}_state_scores.tsv", "composition_all_cohorts.tsv",
         "idh_patients.tsv", "idh_wt_vs_mut.tsv"]


def backup():
    os.makedirs(BAK, exist_ok=True)
    for f in FILES:
        src = os.path.join(MC, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(BAK, f))
    print(f"backed up {len(FILES)} tables -> {BAK}")


def neftel_composition(sub):
    """Per-patient four-state composition of the Smart-seq2 malignant cells."""
    mal = r1.load_cohort(COHORT, sub)          # SS2-only by default (see rev01 docstring)
    comp = mal.obs.groupby("patient", observed=True)["state"].value_counts().unstack(fill_value=0)
    comp_n = comp.div(comp.sum(axis=1), axis=0)
    for s in STATE_NAMES:
        comp_n[f"{s}_prop"] = comp_n.get(s, 0)
    comp_n["cohort"] = COHORT
    comp_n.index.name = "patient"
    scores = mal.obs[["patient", "cohort", "state", "score_MES"]].copy()
    print(f"[rev23] Neftel Smart-seq2: {mal.n_obs} malignant cells, "
          f"{comp_n.shape[0]} patients; mean MES_prop {comp_n['MES_prop'].mean():.3f}")
    return comp_n, scores


def replace_rows(old_path, new_rows, cohort_col="cohort"):
    """Replace the cohort's rows in a table, keeping every other row byte-identical."""
    old = pd.read_csv(old_path, sep="\t")
    mask = (old[cohort_col] == COHORT).to_numpy()
    cols = [c for c in old.columns if c in new_rows.columns]
    if not mask.any():                       # cohort absent: append at the end
        return old, pd.concat([old, new_rows[cols]], ignore_index=True)
    pos = int(np.argmax(mask))
    keep = old[~mask].copy()
    new = pd.concat([keep.iloc[:pos], new_rows[cols], keep.iloc[pos:]], ignore_index=True)
    return old, new


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev23.log"))
    sub = four_state_tables(load_neftel_modules())[1]

    # before-numbers, for the record
    before = pd.read_csv(os.path.join(MC, "neftel_composition.tsv"), sep="\t")
    before_scores = pd.read_csv(os.path.join(MC, f"{COHORT}_state_scores.tsv"), sep="\t")
    print(f"[rev23] as-shipped Neftel: {before.shape[0]} patients, "
          f"{before_scores.shape[0]} malignant cells, mean MES_prop {before['MES_prop'].mean():.3f}")

    backup()
    comp_n, scores = neftel_composition(sub)

    # ---- 1. per-patient composition ----
    comp_n.to_csv(os.path.join(MC, "neftel_composition.tsv"), sep="\t")
    # ---- 2. per-cell scores ----
    scores.to_csv(os.path.join(MC, f"{COHORT}_state_scores.tsv"), sep="\t")

    # ---- 3. cohort-level composition table ----
    old, new = replace_rows(os.path.join(MC, "composition_all_cohorts.tsv"),
                            comp_n.reset_index()[["patient", "cohort", "MES_prop", "NPC_prop",
                                                  "AC_prop", "OPC_prop"]])
    new.to_csv(os.path.join(MC, "composition_all_cohorts.tsv"), sep="\t", index=False)
    print(f"[rev23] composition_all_cohorts.tsv: {old.shape[0]} -> {new.shape[0]} patients "
          f"({int((old['cohort'] == COHORT).sum())} Neftel rows -> "
          f"{int((new['cohort'] == COHORT).sum())})")

    # ---- 4. IDH patient table ----
    idh_rows = new[new["cohort"] == COHORT][["patient", "MES_prop"]].copy()
    idh_rows["cohort"] = COHORT
    idh_rows["genotype"] = "IDH-wt"
    idh_rows = idh_rows[["patient", "cohort", "genotype", "MES_prop"]]
    old_i, new_i = replace_rows(os.path.join(MC, "idh_patients.tsv"), idh_rows)
    new_i.to_csv(os.path.join(MC, "idh_patients.tsv"), sep="\t", index=False)
    print(f"[rev23] idh_patients.tsv: {old_i.shape[0]} -> {new_i.shape[0]} patients "
          f"({int((old_i['cohort'] == COHORT).sum())} Neftel rows -> "
          f"{int((new_i['cohort'] == COHORT).sum())})")

    # ---- 5. IDH-wt vs IDH-mut summary (same schema as the shipped table) ----
    from scipy.stats import mannwhitneyu
    wt = new_i.loc[new_i["genotype"] == "IDH-wt", "MES_prop"]
    mut = new_i.loc[new_i["genotype"] == "IDH-mut", "MES_prop"]
    stat, p = mannwhitneyu(wt, mut, alternative="two-sided")
    res = pd.DataFrame([{
        "comparison": "IDH_wt_vs_IDH_mut", "n_IDHwt": len(wt), "n_IDHmut": len(mut),
        "MES_IDHwt_mean": round(float(wt.mean()), 3), "MES_IDHmut_mean": round(float(mut.mean()), 3),
        "MES_IDHwt_median": round(float(wt.median()), 3),
        "MES_IDHmut_median": round(float(mut.median()), 3),
        "mannwhitney_U": float(stat), "p": float(p)}])
    res.to_csv(os.path.join(MC, "idh_wt_vs_mut.tsv"), sep="\t", index=False)
    print(f"[rev23] IDH-wt {wt.mean():.3f} (n={len(wt)}) vs IDH-mut {mut.mean():.3f} "
          f"(n={len(mut)}): MWU P = {p:.3g}")

    # ---- record ----
    cmp = pd.DataFrame([
        dict(variant="as-shipped (2 platforms, mixed scale)", patients=before.shape[0],
             malignant_cells=before_scores.shape[0], mean_MES_prop=round(float(before['MES_prop'].mean()), 3)),
        dict(variant="Smart-seq2 only (used)", patients=comp_n.shape[0],
             malignant_cells=scores.shape[0], mean_MES_prop=round(float(comp_n['MES_prop'].mean()), 3)),
    ])
    cmp.to_csv(os.path.join(lib.PROJ, "results", "revision", "rev23_neftel_before_after.tsv"),
               sep="\t", index=False)
    print(cmp.to_string(index=False))
    lib.log.info("DONE rev23:\n" + cmp.to_string(index=False))


if __name__ == "__main__":
    main()
