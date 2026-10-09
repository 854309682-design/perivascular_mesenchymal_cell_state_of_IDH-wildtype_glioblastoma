#!/usr/bin/env python3
"""Multi-cohort four-state composition from already-downloaded data.

Processes GSE131928 (Neftel, log-scaled TPM) and GSE70630 (Venteicher, log-scaled
gene x cell) through the same meta-module state scoring used for GSE103224 /
GSE182109, excludes immune/endothelial (non-malignant) cells by marker score, and
assembles a 4-cohort per-patient state-composition table plus an IDH-wt vs IDH-mut
comparison.

Outputs:
  results/multicohort/composition_all_cohorts.tsv
  results/multicohort/neftel_composition.tsv  gse70630_composition.tsv
  results/multicohort/idh_wt_vs_mut.tsv
  results/multicohort/state_composition_by_cohort.png
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "multicohort")
os.makedirs(OUT, exist_ok=True)
NONMALIG_THR = 0.5
IMMUNE = ["PTPRC", "CD68", "AIF1", "C1QA", "C1QB", "C1QC", "LYZ", "CD3D", "CD3E", "NKG7", "IL7R"]
ENDO = ["PECAM1", "VWF", "CLDN5", "ENG", "FLT1"]


def score_state(ad, state_sub):
    for name in ["MES1", "MES2", "NPC1", "NPC2", "AC", "OPC"]:
        genes = state_sub.get(name, [])
        present = [g for g in genes if g in ad.var_names]
        if len(present) >= 5:
            sc.tl.score_genes(ad, present, score_name=f"state_{name}", ctrl_size=50)
    ad.obs["score_MES"] = np.maximum(ad.obs["state_MES1"], ad.obs["state_MES2"])
    ad.obs["score_NPC"] = np.maximum(ad.obs["state_NPC1"], ad.obs["state_NPC2"])
    ad.obs["score_AC"] = ad.obs["state_AC"]
    ad.obs["score_OPC"] = ad.obs["state_OPC"]
    scmat = ad.obs[["score_MES", "score_NPC", "score_AC", "score_OPC"]].values
    ad.obs["state"] = [STATE_NAMES[i] for i in scmat.argmax(axis=1)]
    return ad


def exclude_immune(ad):
    for cname, genes in [("immune", IMMUNE), ("endothelial", ENDO)]:
        present = [g for g in genes if g in ad.var_names]
        sc.tl.score_genes(ad, present, score_name=f"score_{cname}", ctrl_size=50)
    ad.obs["nonmalignant"] = ((ad.obs["score_immune"] > NONMALIG_THR) |
                              (ad.obs["score_endothelial"] > NONMALIG_THR)).astype(int)
    return ad


def load_log_matrix(path, cohort, patient_fn, source):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip("'\"").str.strip()
    ad = sc.AnnData(sp.csr_matrix(df.T.values.astype(np.float32)))
    ad.var_names = df.index.values
    ad.var_names_make_unique()
    ad.obs_names = df.columns.astype(str).values
    ad.obs["cohort"] = cohort
    ad.obs["source"] = source
    ad.obs["patient"] = [patient_fn(c) for c in ad.obs_names]
    ad.obs_names_make_unique()
    return ad


def compose(ad, cohort):
    ad = exclude_immune(ad)
    ad = score_state(ad, state_sub)
    mal = ad[ad.obs["nonmalignant"] == 0].copy()
    comp = mal.obs.groupby("patient")["state"].value_counts().unstack(fill_value=0)
    comp_n = comp.div(comp.sum(axis=1), axis=0)
    for s in STATE_NAMES:
        comp_n[f"{s}_prop"] = comp_n.get(s, 0)
    comp_n["cohort"] = cohort
    comp_n.index.name = "patient"
    mal.obs[["patient", "cohort", "state", "score_MES"]].to_csv(
        os.path.join(OUT, f"{cohort}_state_scores.tsv"), sep="\t")
    return comp_n


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "mc_cohort_state.log"))
    global state_sub
    state_sub = four_state_tables(load_neftel_modules())[1]

    # ---- GSE131928 (Neftel): two log-scaled TPM files (Smartseq2 + 10x) ----
    f1 = lib.NEFTEL_TPM_SS2
    f2 = lib.NEFTEL_TPM_10X
    a1 = load_log_matrix(f1, "Neftel_GSE131928", lambda c: c.split("-")[0], "Smartseq2")
    a2 = load_log_matrix(f2, "Neftel_GSE131928", lambda c: c.split("_")[0], "10x")
    neft = sc.concat([a1, a2], join="outer")
    neft.var_names_make_unique()
    lib.log.info(f"Neftel loaded: {neft.n_obs} cells x {neft.n_vars} genes")
    neft_comp = compose(neft, "Neftel_GSE131928")
    neft_comp.to_csv(os.path.join(OUT, "neftel_composition.tsv"), sep="\t")
    lib.log.info(f"Neftel patients: {neft_comp.shape[0]}; mean MES {neft_comp['MES_prop'].mean():.3f}")

    # ---- GSE70630 (Venteicher, IDH-mut control): log-scaled gene x cell ----
    vf = lib.GSE70630_TXT
    vent = load_log_matrix(vf, "Venteicher_GSE70630", lambda c: c.split("_")[0], "processed")
    lib.log.info(f"Venteicher loaded: {vent.n_obs} cells x {vent.n_vars} genes")
    vent_comp = compose(vent, "Venteicher_GSE70630")
    vent_comp.to_csv(os.path.join(OUT, "gse70630_composition.tsv"), sep="\t")
    lib.log.info(f"Venteicher patients: {vent_comp.shape[0]}; mean MES {vent_comp['MES_prop'].mean():.3f}")

    # ---- combine with GSE103224 (IDH-wt) + GSE182109 ----
    gse103 = pd.read_csv(os.path.join(OUT, "..", "step03_state", "state_composition_by_sample.tsv"),
                         sep="\t", index_col=0)
    gse103["patient"] = gse103.index
    gse103["cohort"] = "GSE103224"
    gse103["MES_prop"] = gse103.get("MES", 0)
    gse103["NPC_prop"] = gse103.get("NPC", 0)
    gse103["AC_prop"] = gse103.get("AC", 0)
    gse103["OPC_prop"] = gse103.get("OPC", 0)
    gse103 = gse103[["patient", "cohort", "MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]]

    g182 = pd.read_csv(os.path.join(OUT, "gse182109_composition_by_patient.tsv"),
                       sep="\t", index_col=[0, 1])
    g182 = g182.reset_index()
    g182 = g182[["patient", "type", "MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]]
    g182["cohort"] = "GSE182109"
    g182 = g182[["patient", "cohort", "MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]].drop_duplicates("patient")

    allc = pd.concat([gse103, g182, neft_comp.reset_index()[["patient", "cohort", "MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]],
                      vent_comp.reset_index()[["patient", "cohort", "MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]]],
                     ignore_index=True)
    allc.to_csv(os.path.join(OUT, "composition_all_cohorts.tsv"), sep="\t", index=False)
    lib.log.info(f"Combined {allc.shape[0]} patients across {allc['cohort'].nunique()} cohorts")

    # ---- IDH-wt vs IDH-mut (MES-like proportion, patient level) ----
    idh_wt = allc[allc["cohort"].isin(["GSE103224", "GSE182109", "Neftel_GSE131928"])]["MES_prop"]
    idh_mut = allc[allc["cohort"] == "Venteicher_GSE70630"]["MES_prop"]
    from scipy.stats import mannwhitneyu
    stat, p = mannwhitneyu(idh_wt, idh_mut, alternative="two-sided")
    res = {"comparison": "IDH_wt_vs_IDH_mut", "n_IDHwt": len(idh_wt), "n_IDHmut": len(idh_mut),
           "MES_IDHwt_mean": round(float(idh_wt.mean()), 3), "MES_IDHmut_mean": round(float(idh_mut.mean()), 3),
           "MES_IDHwt_median": round(float(idh_wt.median()), 3), "MES_IDHmut_median": round(float(idh_mut.median()), 3),
           "mannwhitney_U": float(stat), "p": float(p)}
    pd.DataFrame([res]).to_csv(os.path.join(OUT, "idh_wt_vs_mut.tsv"), sep="\t", index=False)
    lib.log.info(f"IDH-wt {idh_wt.mean():.3f} vs IDH-mut {idh_mut.mean():.3f} (MWU p={p:.4f})")

    # ---- per-cohort mean MES figure ----
    try:
        means = allc.groupby("cohort")["MES_prop"].mean().sort_values()
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(means.index, means.values, color="#C0392B")
        ax.set_ylabel("mean MES-like proportion (patient)"); ax.set_title("Four-state MES-like by cohort")
        plt.xticks(rotation=30, ha="right"); plt.tight_layout()
        plt.savefig(os.path.join(OUT, "state_composition_by_cohort.png"), dpi=150)
        plt.close("all")
    except Exception as e:
        lib.log.warning(f"figure error: {e}")

    lib.log.info("\n" + allc.groupby("cohort")["MES_prop"].agg(["count", "mean", "median"]).round(3).to_string())


if __name__ == "__main__":
    main()
