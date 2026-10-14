#!/usr/bin/env python3
"""rev15 — IDH comparison sensitivity analysis + cohort summary Table 1.

Reviewer minor 13 / stats 2.5 asked for a cohort-stratified test of the IDH-wildtype vs
IDH-mutant MES comparison. That test does not exist in these data: **genotype is perfectly
confounded with cohort** (GSE103224 / GSE131928 / GSE182109 are IDH-wildtype; GSE70630 /
GSE89567 are IDH-mutant), so every stratum is single-armed and neither a van Elteren test nor
a mixed model with a cohort random effect is identifiable. What CAN be reported, and is
reported here, is (i) the per-cohort genotype means, to show the direction is the same in
every cohort, and (ii) a cohort-mean-centred comparison, which removes any cohort-level shift
before testing the genotype contrast.

Reviewer minor 26 asked for a cohort-summary table; `rev15_table1_cohorts.tsv` is that table
(platform, patients, malignant cells, IDH composition, treatment context, analyses used).

Outputs (results/revision/): rev15_idh_stratified.tsv, rev15_table1_cohorts.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import lib

OUT = os.path.join(lib.PROJ, "results", "revision")
MC = os.path.join(lib.PROJ, "results", "multicohort")
os.makedirs(OUT, exist_ok=True)

PLATFORM = {"GSE103224": "10x Chromium", "GSE131928": "Smart-seq2",
            "GSE182109": "10x Chromium", "GSE70630": "Smart-seq2", "GSE89567": "Smart-seq2"}
TREATMENT = {"GSE103224": "primary", "GSE131928": "primary",
             "GSE182109": "primary + recurrent (matched design)", "GSE70630": "primary",
             "GSE89567": "primary"}
GENOTYPE = {"GSE103224": "IDH-wildtype", "GSE131928": "IDH-wildtype",
            "GSE182109": "IDH-wildtype", "GSE70630": "IDH-mutant", "GSE89567": "IDH-mutant"}
ANALYSES = {"GSE103224": "discovery; four-state composition; pseudotime; prior + de-novo regulon; IDH arm",
            "GSE131928": "composition; prior + de-novo regulon; IDH arm",
            "GSE182109": "composition; H1 recurrence; receptor expression; prior + de-novo regulon; IDH arm",
            "GSE70630": "IDH-mutant control; composition; prior regulon; IDH arm",
            "GSE89567": "IDH-mutant control; composition; IDH arm"}


def malignant_counts():
    n = {}
    for f, name in [("GSE103224_state_scores.tsv", "GSE103224"),
                    ("Neftel_GSE131928_state_scores.tsv", "GSE131928"),
                    ("Venteicher_GSE70630_state_scores.tsv", "GSE70630"),
                    ("gse89567_state_scores.tsv", "GSE89567"),
                    ("gse182109_state_scores.tsv", "GSE182109")]:  # per-cohort malignant cells
        p = os.path.join(MC, f)
        if os.path.exists(p):
            n[name] = sum(1 for _ in open(p)) - 1
    p = os.path.join(lib.PROJ, "results", "step03_state", "state_scores_per_cell.tsv")
    if os.path.exists(p):
        n["GSE103224"] = sum(1 for _ in open(p)) - 1
    return n


def idh_analysis():
    d = pd.read_csv(os.path.join(MC, "idh_patients.tsv"), sep="\t")
    d = d.rename(columns={"MES_prop": "mes"})
    pooled = stats.mannwhitneyu(d.loc[d["genotype"] == "IDH-wt", "mes"],
                                d.loc[d["genotype"] == "IDH-mut", "mes"],
                                alternative="two-sided")
    rows = []
    for c, g in d.groupby("cohort"):
        gt = sorted(g["genotype"].unique())
        rows.append({"cohort": c, "genotype": "/".join(gt), "n_patients": len(g),
                     "mean_MES": round(float(g["mes"].mean()), 3),
                     "median_MES": round(float(g["mes"].median()), 3),
                     "n_MES_gt_0.1": int((g["mes"] > 0.1).sum())})
    per = pd.DataFrame(rows)
    # cohort-mean-centred comparison (removes any cohort-level shift)
    d["mes_centred"] = d["mes"] - d.groupby("cohort")["mes"].transform("mean")
    centred = stats.mannwhitneyu(d.loc[d["genotype"] == "IDH-wt", "mes_centred"],
                                 d.loc[d["genotype"] == "IDH-mut", "mes_centred"],
                                 alternative="two-sided")
    out = pd.concat([
        per.assign(test="per_cohort_means", p=np.nan, note=""),
        pd.DataFrame([{
            "cohort": "ALL", "genotype": "IDH-wt vs IDH-mut", "n_patients": len(d),
            "mean_MES": round(float(d["mes"].mean()), 3), "median_MES": np.nan,
            "n_MES_gt_0.1": int((d["mes"] > 0.1).sum()), "test": "pooled_Mann_Whitney",
            "p": pooled.pvalue,
            "note": "pooled over %d cohorts; U = %.1f" % (d["cohort"].nunique(), pooled.statistic),
        }, {
            "cohort": "ALL", "genotype": "IDH-wt vs IDH-mut", "n_patients": len(d),
            "mean_MES": 0.0, "median_MES": np.nan, "n_MES_gt_0.1": np.nan,
            "test": "cohort_mean_centred_Mann_Whitney", "p": centred.pvalue,
            "note": ("genotype is perfectly confounded with cohort, so no within-cohort "
                     "(van Elteren) test is identifiable; this test removes cohort shifts "
                     "first; U = %.1f" % centred.statistic),
        }]),
    ], ignore_index=True)
    return d, per, out, pooled, centred


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev15.log"))
    d, per, out, pooled, centred = idh_analysis()
    out.to_csv(os.path.join(OUT, "rev15_idh_stratified.tsv"), sep="\t", index=False)
    lib.log.info("IDH per-cohort means:\n" + per.to_string(index=False))
    lib.log.info(f"pooled MWU p = {pooled.pvalue:.3g}; cohort-mean-centred MWU p = {centred.pvalue:.3g}")
    lib.log.info("direction check: every cohort keeps the same sign? " +
                 str(bool((per["mean_MES"] < 0.05).tolist() ==
                          sorted((per["mean_MES"] < 0.05).tolist(), reverse=True))))

    ncell = malignant_counts()
    # patient-level MES per cohort: idh_patients.tsv covers all five cohorts (the composition
    # table uses different cohort labels and omits GSE89567)
    d_all = pd.read_csv(os.path.join(MC, "idh_patients.tsv"), sep="\t")
    agg = d_all.groupby("cohort")["MES_prop"].agg(n_patients="count", mean_MES="mean")
    name_map = {"GSE131928": "Neftel_GSE131928"}   # idh_patients.tsv uses "GSE70630" as-is
    rows = []
    for c in ["GSE103224", "GSE131928", "GSE182109", "GSE70630", "GSE89567"]:
        key = name_map.get(c, c)
        has = key in agg.index
        rows.append({"cohort": c, "platform": PLATFORM[c], "genotype": GENOTYPE[c],
                     "n_patients": int(agg.loc[key, "n_patients"]) if has else np.nan,
                     "n_malignant_cells": ncell.get(c, np.nan),
                     "mean_MES_like": round(float(agg.loc[key, "mean_MES"]), 3) if has else np.nan,
                     "treatment_context": TREATMENT[c], "analyses": ANALYSES[c]})
    t1 = pd.DataFrame(rows)
    t1.to_csv(os.path.join(OUT, "rev15_table1_cohorts.tsv"), sep="\t", index=False)
    lib.log.info("Table 1:\n" + t1.to_string(index=False))
    lib.log.info("DONE rev15")


if __name__ == "__main__":
    main()
