#!/usr/bin/env python3
"""Multi-cohort H1: GSE182109 (Abdelfattah 2022) primary-vs-recurrent MES-like shift.

Lean pipeline (no clustering/PCA/UMAP): load 44 Cell Ranger samples, QC, normalize
(log CPM/10k), score immune/endothelial markers to exclude non-malignant cells,
score the Neftel meta-modules, assign 4-state labels, and quantify the MES-like
proportion per patient (GBM only) aggregated by type (ndGBM=primary, rGBM=recurrent).
Tests H1 (recurrent > primary) at the patient level.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, fontsize=8, frameon=False, facecolor="white")

BASE = os.path.join(lib.COHORTS, "gse182109")
OUT = os.path.join(lib.PROJ, "results", "multicohort")
os.makedirs(OUT, exist_ok=True)
SEED = 42
MIN_GENES = 200
NONMALIG_THR = 0.5
IMMUNE = ["PTPRC", "CD68", "AIF1", "C1QA", "C1QB", "C1QC", "LYZ", "CD3D", "CD3E", "NKG7", "IL7R"]
ENDO = ["PECAM1", "VWF", "CLDN5", "ENG", "FLT1"]


def parse_name(fname):
    base = fname.split("_", 1)[1]
    name = base.replace("_matrix.mtx.gz", "")
    tok = name.split("-")
    ctype = tok[0]
    patient = "-".join(tok[:2])
    region = "-".join(tok[2:]) if len(tok) > 2 else "NA"
    return name, ctype, patient, region


def load_one(mtxf):
    sfx = mtxf[:-len("_matrix.mtx.gz")]
    feat = pd.read_csv(sfx + "_features.tsv.gz", sep="\t", header=None)
    genes = feat[1].astype(str).values
    m = scipy.io.mmread(mtxf)
    m = sp.csr_matrix(m.T).astype(np.float32)
    bcs = (pd.read_csv(sfx + "_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).values)
    ad = sc.AnnData(m)
    ad.obs_names = [f"{os.path.basename(sfx).split('_',1)[1]}_{b}" for b in bcs]
    ad.var_names = genes
    ad.var_names_make_unique()
    return ad


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "mc_gse182109.log"))
    state_sub = four_state_tables(load_neftel_modules())[1]

    samples = sorted(glob.glob(os.path.join(BASE, "*_matrix.mtx.gz")))
    lib.log.info(f"Found {len(samples)} samples")

    adatas, meta_rows = [], []
    for f in samples:
        name, ctype, patient, region = parse_name(os.path.basename(f))
        ad = load_one(f)
        ad.obs["sample"] = name; ad.obs["type"] = ctype
        ad.obs["patient"] = patient; ad.obs["region"] = region
        adatas.append(ad)
        meta_rows.append({"sample": name, "type": ctype, "patient": patient,
                          "region": region, "n_cells": ad.n_obs})
    pd.DataFrame(meta_rows).to_csv(os.path.join(OUT, "gse182109_sample_meta.tsv"), sep="\t", index=False)

    ad = sc.concat(adatas, join="outer")
    ad.var_names_make_unique()
    lib.log.info(f"Concatenated: {ad.n_obs} cells x {ad.n_vars} genes")

    sc.pp.filter_cells(ad, min_genes=MIN_GENES)
    sc.pp.filter_genes(ad, min_cells=10)
    lib.log.info(f"After QC: {ad.n_obs} x {ad.n_vars}")

    # normalize
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # immune / endothelial marker scores (identify non-malignant)
    for cname, genes in [("immune", IMMUNE), ("endothelial", ENDO)]:
        present = [g for g in genes if g in ad.var_names]
        sc.tl.score_genes(ad, present, score_name=f"score_{cname}", ctrl_size=50)
    ad.obs["nonmalignant"] = ((ad.obs["score_immune"] > NONMALIG_THR) |
                              (ad.obs["score_endothelial"] > NONMALIG_THR)).astype(int)
    lib.log.info(f"Non-malignant (immune/endo) fraction: {ad.obs['nonmalignant'].mean():.3f}")

    # score six meta-modules
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

    # malignant cells for composition
    mal = ad[ad.obs["nonmalignant"] == 0].copy()
    # GBM only (ndGBM + rGBM; LGG excluded from the H1 test but reported)
    gbm = mal[mal.obs["type"].isin(["ndGBM", "rGBM"])]
    comp = gbm.obs.groupby(["patient", "type"])["state"].value_counts().unstack(fill_value=0)
    comp_n = comp.div(comp.sum(axis=1), axis=0)
    comp_n["MES_prop"] = comp_n.get("MES", 0)
    comp_n["NPC_prop"] = comp_n.get("NPC", 0)
    comp_n["AC_prop"] = comp_n.get("AC", 0)
    comp_n["OPC_prop"] = comp_n.get("OPC", 0)
    comp_n["type"] = comp_n.index.map(lambda p: p[1])
    comp_n.to_csv(os.path.join(OUT, "gse182109_composition_by_patient.tsv"), sep="\t")

    # H1: recurrent (rGBM) vs primary (ndGBM) MES-like proportion, patient level
    prim = comp_n[comp_n["type"] == "ndGBM"]["MES_prop"]
    rec = comp_n[comp_n["type"] == "rGBM"]["MES_prop"]
    stat, p = mannwhitneyu(rec, prim, alternative="greater")
    res = {"comparison": "recurrent_rGBM_vs_primary_ndGBM",
           "n_primary": len(prim), "n_recurrent": len(rec),
           "MES_prop_primary_mean": round(float(prim.mean()), 3),
           "MES_prop_recurrent_mean": round(float(rec.mean()), 3),
           "MES_prop_primary_median": round(float(prim.median()), 3),
           "MES_prop_recurrent_median": round(float(rec.median()), 3),
           "mannwhitney_U": float(stat), "p_alt_greater": float(p)}
    pd.DataFrame([res]).to_csv(os.path.join(OUT, "gse182109_h1_primary_vs_recurrent.tsv"),
                               sep="\t", index=False)
    lib.log.info(f"H1: recurrent {rec.mean():.3f} vs primary {prim.mean():.3f} "
                 f"(one-sided MWU p={p:.4f}); n_primary={len(prim)}, n_recurrent={len(rec)}")
    lib.log.info("\n" + comp_n[["MES_prop", "NPC_prop", "AC_prop", "OPC_prop"]].round(3).to_string())

    # figure
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        order = ["ndGBM", "rGBM"]
        xs = np.arange(len(order))
        means = [comp_n[comp_n["type"] == t]["MES_prop"].mean() for t in order]
        sems = [comp_n[comp_n["type"] == t]["MES_prop"].std() /
                np.sqrt(max(1, len(comp_n[comp_n["type"] == t]))) for t in order]
        ax.bar(xs, means, yerr=np.array(sems) * 1.96, capsize=5, color=["#2980B9", "#C0392B"])
        ax.set_xticks(xs); ax.set_xticklabels(["primary\n(ndGBM)", "recurrent\n(rGBM)"])
        ax.set_ylabel("MES-like proportion (patient mean)"); ax.set_title("H1: recurrent vs primary (GSE182109)")
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "gse182109_mes_by_type.png"), dpi=150)
        plt.close("all")
    except Exception as e:
        lib.log.warning(f"figure error: {e}")

    mal.obs[["patient", "type", "state", "score_MES"]].to_csv(
        os.path.join(OUT, "gse182109_state_scores.tsv"), sep="\t")
    lib.log.info("Done.")


if __name__ == "__main__":
    main()
