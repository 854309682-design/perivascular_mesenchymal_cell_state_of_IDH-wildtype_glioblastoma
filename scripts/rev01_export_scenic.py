#!/usr/bin/env python3
"""rev01 — Export state-stratified malignant-cell looms for the de-novo pySCENIC arm.

Reviewer item M2/DA2/DA3: the prior-independent arm behind Fig. 3b and Results §3.3
could not be reproduced (the Sep-11..15 workspace is missing). This script rebuilds
its *inputs* in a documented, reproducible way, for the same four cohorts in which
the manuscript ran both inference routes:

    GSE103224       (Yuan 2018,   IDH-wt GBM, 10x)
    Neftel_GSE131928(Neftel 2019, IDH-wt GBM, **Smart-seq2 block only** — see below)
    GSE182109       (Abdelfattah 2022, IDH-wt GBM, 10x, primary + recurrent)
    GSE70630        (Tirosh 2016, IDH-mut oligodendroglioma)

Cell selection mirrors the original Methods text: malignant cells only (immune /
endothelial marker-score exclusion), stratified subsample across the four Neftel
states so that rare MES cells are not lost, capped per state and in total.

**Neftel is restricted to its Smart-seq2 block (2026-10-14).** The two GEO matrices for
that cohort are on different scales — the Smart-seq2 file is log2(TPM+1) (nonzero values
0.02-10) while the 10x file is raw TPM (86-691) — so concatenating them (which the
submitted analysis did; see release/scripts/mc_cohort_state.py) makes the per-cell MES
score and any pseudotime axis mostly a *platform* axis. Measured on the as-shipped
matrices: rho(platform, MES score) = +0.82 and rho(platform, pseudotime) = -0.90; with the
Smart-seq2 block alone the MES-pseudotime correlation is ~0 (rho = +0.04). Set
`NEFTEL_SS2_ONLY=0` to reproduce the superseded two-platform behaviour.

**Documented deviations from the submitted Methods text** (the submitted run's code no
longer exists, so this re-run's configuration is stated explicitly and is reported in
the revised manuscript): the subsample is capped at 500 cells per state / 2,000 cells
in total rather than 10,000 — GRNBoost2 fits one model per *target gene*, so its
wall-clock cost is linear in cells x genes — genes detected in fewer than 5% of the
subsampled cells are dropped before export, and Neftel is single-platform. These changes
affect the tractability of the network search and the platform confounding of one cohort.

Outputs (results/scenic/):
    rev01_<cohort>.loom        genes x cells, log-normalized (float32)
    rev01_<cohort>_cellmeta.tsv CellID, state, patient, score_MES, score_NPC/AC/OPC
    rev01_subsample.tsv        cells per cohort x state

Run (Python 3.11 env with scanpy + loompy + pySCENIC):
  GLIOMA_PROJ=/run/media/meng/Seqdata/glioma/Paper_2 \
  GLIOMA_DATA=/run/media/meng/Seqdata/glioma/Datasets \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev01_export_scenic.py
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
import anndata as ad_mod
import loompy as lp

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "scenic")
os.makedirs(OUT, exist_ok=True)

SEED = 42
PER_STATE_CAP = 500
TOTAL_CAP = 2000
MIN_CELL_FRAC = 0.05   # keep genes detected in >=5% of the subsampled cells (GRN tractability;
                       # documented deviation - see the module docstring and the Methods rewrite)
NONMALIG_THR = 0.5
MIN_GENES = 200
IMMUNE = ["PTPRC", "CD68", "AIF1", "C1QA", "C1QB", "C1QC", "LYZ", "CD3D", "CD3E", "NKG7", "IL7R"]
ENDO = ["PECAM1", "VWF", "CLDN5", "ENG", "FLT1"]


def score_states(ad, state_sub):
    """Marker-score exclusion + Neftel six-module state scoring (same as mc_cohort_state)."""
    for cname, genes in [("immune", IMMUNE), ("endothelial", ENDO)]:
        present = [g for g in genes if g in ad.var_names]
        sc.tl.score_genes(ad, present, score_name=f"score_{cname}", ctrl_size=50)
    ad.obs["nonmalignant"] = ((ad.obs["score_immune"] > NONMALIG_THR) |
                              (ad.obs["score_endothelial"] > NONMALIG_THR)).astype(int)
    for name in ["MES1", "MES2", "NPC1", "NPC2", "AC", "OPC"]:
        present = [g for g in state_sub.get(name, []) if g in ad.var_names]
        if len(present) >= 5:
            sc.tl.score_genes(ad, present, score_name=f"state_{name}", ctrl_size=50)
    ad.obs["score_MES"] = np.maximum(ad.obs["state_MES1"], ad.obs["state_MES2"])
    ad.obs["score_NPC"] = np.maximum(ad.obs["state_NPC1"], ad.obs["state_NPC2"])
    ad.obs["score_AC"] = ad.obs["state_AC"]
    ad.obs["score_OPC"] = ad.obs["state_OPC"]
    m = ad.obs[["score_MES", "score_NPC", "score_AC", "score_OPC"]].values
    ad.obs["state"] = [STATE_NAMES[i] for i in m.argmax(axis=1)]
    return ad


def load_log_matrix(path, cohort, patient_fn, source):
    """gene x cell log-scale matrix (same reader as mc_cohort_state.load_log_matrix)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip("'\"").str.strip()
    a = sc.AnnData(sp.csr_matrix(df.T.values.astype(np.float32)))
    a.var_names = df.index.values
    a.var_names_make_unique()
    a.obs_names = df.columns.astype(str).values
    a.obs["cohort"] = cohort
    a.obs["source"] = source
    a.obs["patient"] = [patient_fn(c) for c in a.obs_names]
    a.obs_names_make_unique()
    return a


def load_cohort(name, state_sub):
    """Return the malignant-cell AnnData (log-normalized) for one cohort."""
    if name == "GSE103224":
        p = os.path.join(lib.PROC_DIR, "gse103224_regulon.h5ad")
        a = sc.read_h5ad(p)
        a.obs["cohort"] = "GSE103224"
        a.obs["patient"] = a.obs["sample"].astype(str)
        # 05_regulon.py saved the four state scores but not the argmax label: re-derive it
        # with the project's documented rule (max of the four pooled meta-module scores).
        m = a.obs[["score_MES", "score_NPC", "score_AC", "score_OPC"]].values
        a.obs["state"] = pd.Series([STATE_NAMES[i] for i in m.argmax(axis=1)],
                                   index=a.obs_names, dtype=object)
        return a                                   # already malignant + log-normalized

    if name == "Neftel_GSE131928":
        # 2026-10-14: the two GEO matrices for this cohort are on DIFFERENT scales — the
        # Smart-seq2 file is log2(TPM+1) (nonzero 0.02-10) while the 10x file is raw TPM
        # (86-691). Concatenating them makes the per-cell MES score and any pseudotime axis
        # mostly a platform axis (measured: rho(platform, MES) = +0.82, rho(platform, dpt) =
        # -0.90; see handoff.md and the probe in logs/rev01_neftel_scale.log). The cohort is
        # therefore restricted to the Smart-seq2 block, which is internally consistent and is
        # also the platform that carries the per-patient design. Set NEFTEL_SS2_ONLY=0 to
        # restore the old two-platform concatenation (not used in the revised manuscript).
        a1 = load_log_matrix(lib.NEFTEL_TPM_SS2, name, lambda c: c.split("-")[0], "Smartseq2")
        if os.environ.get("NEFTEL_SS2_ONLY", "1") == "0":
            a2 = load_log_matrix(lib.NEFTEL_TPM_10X, name, lambda c: c.split("_")[0], "10x")
            a = sc.concat([a1, a2], join="outer")
            a.var_names_make_unique()
            lib.log.warning("[Neftel_GSE131928] two-platform concatenation re-enabled "
                            "(mixed scales) — not used in the revised manuscript")
        else:
            a = a1
            lib.log.info("[Neftel_GSE131928] restricted to the Smart-seq2 block "
                         "(the 10x file is on a different scale; see the module docstring)")

    elif name == "GSE182109":
        base = os.path.join(lib.COHORTS, "gse182109")
        files = sorted(glob.glob(os.path.join(base, "*_matrix.mtx.gz")))
        ads = []
        for f in files:
            sfx = f[:-len("_matrix.mtx.gz")]
            genes = pd.read_csv(sfx + "_features.tsv.gz", sep="\t", header=None)[1].astype(str).values
            m = sp.csr_matrix(scipy.io.mmread(f).T).astype(np.float32)
            bcs = pd.read_csv(sfx + "_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).values
            smp = os.path.basename(sfx).split("_", 1)[1]
            x = sc.AnnData(m)
            x.obs_names = [f"{smp}_{b}" for b in bcs]
            x.var_names = genes
            x.var_names_make_unique()
            tok = smp.split("-")
            x.obs["cohort"] = name
            x.obs["sample"] = smp
            x.obs["type"] = tok[0]
            x.obs["patient"] = "-".join(tok[:2])
            ads.append(x)
        a = sc.concat(ads, join="outer")
        a.var_names_make_unique()
        sc.pp.filter_cells(a, min_genes=MIN_GENES)
        sc.pp.filter_genes(a, min_cells=10)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

    elif name == "GSE70630":
        a = load_log_matrix(lib.GSE70630_TXT, name, lambda c: c.split("_")[0], "processed")

    else:
        raise ValueError(name)

    a = score_states(a, state_sub)
    mal = a[a.obs["nonmalignant"] == 0].copy()
    lib.log.info(f"[{name}] {a.n_obs} cells -> {mal.n_obs} malignant")
    return mal


def stratified_subsample(mal, rng):
    picked = []
    for st, grp in mal.obs.groupby("state", observed=True):
        idx = grp.index.to_numpy()
        if len(idx) > PER_STATE_CAP:
            idx = rng.choice(idx, PER_STATE_CAP, replace=False)
        picked.extend(idx.tolist())
    if len(picked) > TOTAL_CAP:
        picked = rng.choice(np.array(picked), TOTAL_CAP, replace=False).tolist()
    return list(dict.fromkeys(picked))


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev01_export_scenic.log"))
    state_sub = four_state_tables(load_neftel_modules())[1]
    rng = np.random.default_rng(SEED)
    rows = []

    cohorts = sys.argv[1:] or ["GSE103224", "GSE182109", "Neftel_GSE131928", "GSE70630"]
    for cohort in cohorts:
        mal = load_cohort(cohort, state_sub)
        picked = stratified_subsample(mal, rng)
        sub = mal[picked].copy()
        x = sub.X
        x = np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x)
        x = np.nan_to_num(x.astype(np.float32))
        # gene filter: GRNBoost2 fits one model per target gene, so the wall-clock cost is
        # linear in the gene count. Genes detected in fewer than MIN_CELL_FRAC of the
        # subsampled cells cannot contribute a stable co-expression module.
        keep = (x > 0).sum(axis=0) >= MIN_CELL_FRAC * x.shape[0]
        x = x[:, keep]
        sub = sub[:, keep].copy()
        lo = os.path.join(OUT, f"rev01_{cohort}.loom")
        lp.create(lo, x.T,
                  row_attrs={"Gene": np.array(sub.var_names, dtype=str)},
                  col_attrs={"CellID": np.array(sub.obs_names, dtype=str)})
        meta = sub.obs[["state", "patient", "score_MES", "score_NPC", "score_AC", "score_OPC"]].copy()
        meta.index.name = "CellID"
        meta.to_csv(os.path.join(OUT, f"rev01_{cohort}_cellmeta.tsv"), sep="\t")
        vc = sub.obs["state"].value_counts().to_dict()
        lib.log.info(f"[{cohort}] loom {x.shape[1]} genes x {x.shape[0]} cells -> {lo} | states {vc}")
        for st, n in vc.items():
            rows.append({"cohort": cohort, "state": st, "n_cells": n,
                         "n_malignant_total": int(mal.n_obs)})
    # merge into the audit table instead of overwriting it, so a single-cohort re-export
    # (e.g. Neftel after the platform correction) does not erase the other cohorts' rows
    out = os.path.join(OUT, "rev01_subsample.tsv")
    new = pd.DataFrame(rows)
    if os.path.exists(out):
        old = pd.read_csv(out, sep="\t")
        old = old[~old["cohort"].isin(new["cohort"].unique())]
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(out, sep="\t", index=False)
    lib.log.info("DONE rev01")


if __name__ == "__main__":
    main()
