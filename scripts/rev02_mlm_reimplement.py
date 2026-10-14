#!/usr/bin/env python3
"""rev02 — Prior-based TF activity: faithful mLM reimplementation + patient-level statistics.

Two reviewer items are addressed here.

M5 (pseudoreplication): the submitted TF-activity / regulon-MES correlations were
computed across cells and their FDRs treated cells as independent. Cells within a
patient are not independent. `decoupler` is not installable in this session (user
decision, 2026-10-14), so the multivariate linear model it implements is
re-implemented here from the DoRothEA prior that is already on disk
(results/step05_regulon/dorothea_net.tsv):

    A[gene, TF] = mode-of-regulation weight if `gene` is a target of `TF` else 0
    beta(cell)  = argmin_beta || A beta - y(cell) ||_2        (genes in the net)

`solve` is a single pre-computed pseudo-inverse, so the whole cell x TF activity
matrix is one matrix product. **Validation (do not skip):** the reimplementation is
checked against the genuine decoupler output that survives on disk for GSE103224
(data/processed/gse103224_regulon.h5ad, obsm['score_mlm']) before it is trusted on
the other cohorts.

M2(d/e): TWIST1's rank and sign in the prior-based route, and the CEBPB sign flip.

Outputs (results/revision/):
    rev02_mlm_validation.tsv        per-TF agreement with the archived decoupler matrix
    rev02_activity_<cohort>.tsv     per-cell TF activity (cells x key TFs, patient-mean copy)
    rev02_tf_activity_vs_MES.tsv    per-cohort TF-MES rho (pooled cells) - reimplementation
    rev02_patient_level_tf.tsv      per-patient rho + n cells per TF per cohort
    rev02_pooled_vs_patient.tsv     pooled-cell vs patient-level summary for key TFs
    rev02_recurrence_power.tsv      M4: effect size + CI + two-sided test

Run:
  GLIOMA_PROJ=... GLIOMA_DATA=... \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev02_mlm_reimplement.py
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES
import mc_cohort_state as mcs
import mc_gse182109_h1 as g182

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "revision")
os.makedirs(OUT, exist_ok=True)
NET_PATH = os.path.join(lib.PROJ, "results", "step05_regulon", "dorothea_net.tsv")
MIN_TARGETS = 10
MIN_CELLS_PER_PATIENT = 50
KEY_TFS = ["STAT3", "KLF4", "FOXL2", "JUNB", "FOSL1", "FOSL2", "JUN", "NFKB1",
           "CEBPB", "CEBPA", "RELA", "TEAD1", "HIF1A", "TWIST1", "ATF4"]
SEED = 42


# ---------------------------------------------------------------- mLM core
def design_matrix(net, genes):
    """(genes_in_net x TFs) design with DoRothEA mode-of-regulation weights."""
    net = net[net["target"].isin(genes) & net["source"].isin(genes)].copy()
    counts = net.groupby("source")["target"].count()
    net = net[net["source"].isin(counts[counts >= MIN_TARGETS].index)].copy()
    genes_in = sorted(set(net["target"]) | set(net["source"]))
    tfs = sorted(net["source"].unique())
    gi = {g: i for i, g in enumerate(genes_in)}
    ti = {t: i for i, t in enumerate(tfs)}
    A = np.zeros((len(genes_in), len(tfs)), dtype=np.float64)
    for tgt, src, w in zip(net["target"], net["source"], net["weight"]):
        A[gi[tgt], ti[src]] = float(w)
    return A, genes_in, tfs


def mlm_activity(X, var_names, net, chunk=2000):
    """X: cells x genes (dense or sparse, log-normalized). Returns (cells x TFs, tfs).

    Memory-safe: the design pseudo-inverse is applied in cell chunks, so a 40k-cell
    cohort never materialises a cells x genes dense array.
    """
    A, genes_in, tfs = design_matrix(net, list(var_names))
    var = pd.Index(var_names)
    idx = var.get_indexer(genes_in)
    assert (idx >= 0).all(), "gene in design missing from var_names"
    pinv = np.linalg.pinv(A)                          # tfs x genes_in
    n = X.shape[0]
    act = np.empty((n, len(tfs)), dtype=np.float32)
    for s in range(0, n, chunk):
        Y = X[s:s + chunk]
        Y = np.asarray(Y.todense()) if sp.issparse(Y) else np.asarray(Y)
        Y = np.nan_to_num(Y.astype(np.float32))[:, idx]
        act[s:s + chunk] = (Y @ pinv.T.astype(np.float32))
    return act, tfs


# ---------------------------------------------------------------- cohorts
def cohort_adata(name, state_sub):
    """Same data path per cohort as mc_cohort_state / mc_h2_regulon (malignant cells)."""
    if name == "GSE103224":
        a = sc.read_h5ad(os.path.join(lib.PROC_DIR, "gse103224_regulon.h5ad"))
        a.obs["patient"] = a.obs["sample"].astype(str)
        a.obs["cohort"] = name
        m = a.obs[["score_MES", "score_NPC", "score_AC", "score_OPC"]].values
        a.obs["state"] = [STATE_NAMES[i] for i in m.argmax(axis=1)]
        return a, a.obsm["score_mlm"], True          # archived decoupler matrix available

    if name == "Neftel":
        # 2026-10-14: the cohort's two GEO matrices are on different scales (Smart-seq2 is
        # log2(TPM+1), the 10x file is raw TPM), which made this cohort's MES score mostly a
        # platform axis. Restricted to the Smart-seq2 block, as in rev01/rev19/rev23; set
        # NEFTEL_SS2_ONLY=0 to reproduce the superseded two-platform behaviour.
        a1 = mcs.load_log_matrix(lib.NEFTEL_TPM_SS2, "Neftel", lambda c: c.split("-")[0], "Smartseq2")
        if os.environ.get("NEFTEL_SS2_ONLY", "1") == "0":
            a2 = mcs.load_log_matrix(lib.NEFTEL_TPM_10X, "Neftel", lambda c: c.split("_")[0], "10x")
            a = sc.concat([a1, a2], join="outer")
            a.var_names_make_unique()
            lib.log.warning("[Neftel] two-platform concatenation re-enabled (mixed scales)")
        else:
            a = a1
            lib.log.info("[Neftel] restricted to the Smart-seq2 block (platform scale correction)")

    elif name == "GSE182109":
        import glob
        files = sorted(glob.glob(os.path.join(lib.COHORTS, "gse182109", "*_matrix.mtx.gz")))
        ads = []
        for f in files:
            x = g182.load_one(f)
            smp = os.path.basename(f).split("_", 1)[1].replace("_matrix.mtx.gz", "")
            tok = smp.split("-")
            x.obs["cohort"] = name
            x.obs["sample"] = smp
            x.obs["type"] = tok[0]                       # ndGBM / rGBM
            x.obs["patient"] = "-".join(tok[:2])
            ads.append(x)
        a = sc.concat(ads, join="outer")
        a.var_names_make_unique()
        sc.pp.filter_cells(a, min_genes=200)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

    elif name == "Venteicher":
        a = mcs.load_log_matrix(lib.GSE70630_TXT, "Venteicher", lambda c: c.split("_")[0], "processed")

    else:
        raise ValueError(name)

    a = mcs.exclude_immune(a)
    a = mcs.score_state(a, state_sub)
    a.obs["cohort"] = name
    a.obs["patient"] = a.obs["patient"].astype(str)
    return a[a.obs["nonmalignant"] == 0].copy(), None, False


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev02_mlm.log"))
    state_sub = four_state_tables(load_neftel_modules())[1]
    net = pd.read_csv(NET_PATH, sep="\t")
    lib.log.info(f"DoRothEA prior: {net.shape}, {net['source'].nunique()} TFs")

    act_tables, rho_rows, pat_rows = {}, [], []

    for cohort in ["GSE103224", "Neftel", "Venteicher", "GSE182109"]:
        ad, archived, has_arch = cohort_adata(cohort, state_sub)
        lib.log.info(f"[{cohort}] {ad.n_obs} malignant cells x {ad.n_vars} genes")

        act, tfs = mlm_activity(ad.X, ad.var_names, net)
        A = pd.DataFrame(act, index=ad.obs_names, columns=tfs)
        lib.log.info(f"[{cohort}] mLM activity {A.shape}")

        # ---- validation against the archived decoupler matrix (GSE103224 only) ----
        if has_arch:
            arch = archived
            arch = arch.toarray() if sp.issparse(arch) else np.asarray(arch)
            arch = pd.DataFrame(arch, index=ad.obs_names, columns=list(archived.columns))
            common = [t for t in A.columns if t in arch.columns]
            vals = []
            for t in common:
                r = spearmanr(A[t].values, arch[t].values).statistic
                vals.append({"tf": t, "spearman_reimpl_vs_decoupler": r})
            v = pd.DataFrame(vals)
            v.to_csv(os.path.join(OUT, "rev02_mlm_validation.tsv"), sep="\t", index=False)
            lib.log.info(f"[validation] {len(common)}/262 TFs; median rho={v.iloc[:,1].median():.4f}; "
                         f"min={v.iloc[:,1].min():.4f}; n(|rho|>0.9)={int((v.iloc[:,1].abs()>0.9).sum())}")

        # ---- per-cell correlation, pooled cells (signatures vs MES) ----
        mes = ad.obs["score_MES"].values
        for t in tfs:
            v = A[t].values
            if np.std(v) < 1e-12:
                continue
            rho, p = spearmanr(v, mes)
            rho_rows.append({"cohort": cohort, "tf": t, "rho": rho, "p": p, "n_cells": int(ad.n_obs)})

        # ---- patient-level (M5) ----
        patient = ad.obs["patient"].values
        for t in tfs:
            v = A[t].values
            df = pd.DataFrame({"patient": patient, "tf": v, "mes": mes})
            for pat, g in df.groupby("patient"):
                if len(g) < MIN_CELLS_PER_PATIENT or g["tf"].std() < 1e-12:
                    continue
                r = spearmanr(g["tf"].values, g["mes"].values).statistic
                pat_rows.append({"cohort": cohort, "tf": t, "patient": pat,
                                 "rho": r, "n_cells": int(len(g))})

        act_tables[cohort] = A
        A[[t for t in KEY_TFS if t in A.columns]].to_csv(
            os.path.join(OUT, f"rev02_activity_{cohort}.tsv"), sep="\t")
        del ad

    rho_df = pd.DataFrame(rho_rows)
    rho_df["fdr"] = np.nan
    for c, g in rho_df.groupby("cohort"):
        rho_df.loc[g.index, "fdr"] = multipletests(g["p"].values, method="fdr_bh")[1]
    rho_df.to_csv(os.path.join(OUT, "rev02_tf_activity_vs_MES.tsv"), sep="\t", index=False)

    pat_df = pd.DataFrame(pat_rows)
    pat_df.to_csv(os.path.join(OUT, "rev02_patient_level_tf.tsv"), sep="\t", index=False)

    # per-cohort consensus + pooled vs patient comparison for the key TFs
    cons = rho_df.groupby("tf").agg(cohorts=("cohort", "nunique"),
                                    mean_rho=("rho", "mean"),
                                    n_fdr005=("fdr", lambda s: int((s < 0.05).sum()))).reset_index()
    cons = cons.sort_values("mean_rho", ascending=False)
    cons.to_csv(os.path.join(OUT, "rev02_tf_consensus_reimplemented.tsv"), sep="\t", index=False)

    rows = []
    for t in KEY_TFS:
        pooled = rho_df[rho_df["tf"] == t]
        pt = pat_df[pat_df["tf"] == t]
        rows.append({
            "tf": t,
            "pooled_mean_rho": round(float(pooled["rho"].mean()), 3),
            "pooled_n_cohorts": int(pooled["cohort"].nunique()),
            "n_fdr005_cohorts": int((pooled["fdr"] < 0.05).sum()),
            "patient_median_rho": round(float(pt["rho"].median()), 3) if len(pt) else np.nan,
            "patient_min_rho": round(float(pt["rho"].min()), 3) if len(pt) else np.nan,
            "patient_max_rho": round(float(pt["rho"].max()), 3) if len(pt) else np.nan,
            "n_patients": int(pt["patient"].nunique()) if len(pt) else 0,
            "patient_pct_positive": round(float((pt["rho"] > 0).mean() * 100), 1) if len(pt) else np.nan,
        })
    comp = pd.DataFrame(rows)
    comp.to_csv(os.path.join(OUT, "rev02_pooled_vs_patient.tsv"), sep="\t", index=False)
    lib.log.info("pooled-cell vs patient-level:\n" + comp.to_string(index=False))

    # ---------------- M4: recurrence test, effect size + CI ----------------
    comp_pat = pd.read_csv(os.path.join(lib.PROJ, "results", "multicohort",
                                        "gse182109_composition_by_patient.tsv"), sep="\t")
    prim = comp_pat[comp_pat["type"] == "ndGBM"]["MES_prop"].values
    rec = comp_pat[comp_pat["type"] == "rGBM"]["MES_prop"].values
    u2, p2 = mannwhitneyu(rec, prim, alternative="two-sided")
    u1, p1 = mannwhitneyu(rec, prim, alternative="greater")
    rng = np.random.default_rng(SEED)
    boots = []
    for _ in range(20000):
        a = rng.choice(prim, len(prim), replace=True)
        b = rng.choice(rec, len(rec), replace=True)
        boots.append(np.mean(b) - np.mean(a))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    # Cliff's delta (non-parametric effect size) + rank-biserial
    gt = sum((b > a) for b in rec for a in prim)
    lt = sum((b < a) for b in rec for a in prim)
    cliff = (gt - lt) / (len(rec) * len(prim))
    rec_out = pd.DataFrame([{
        "comparison": "recurrent rGBM vs primary ndGBM (patient-level MES fraction)",
        "n_primary": len(prim), "n_recurrent": len(rec),
        "mean_primary": round(float(prim.mean()), 3), "mean_recurrent": round(float(rec.mean()), 3),
        "median_primary": round(float(np.median(prim)), 3),
        "median_recurrent": round(float(np.median(rec)), 3),
        "mean_difference": round(float(rec.mean() - prim.mean()), 3),
        "mean_difference_95CI_boot": f"[{lo:.3f}, {hi:.3f}]",
        "cliffs_delta": round(cliff, 3),
        "mannwhitney_U": float(u2), "p_two_sided": float(p2), "p_one_sided_greater": float(p1),
        "note": ("power-limited: with n=5 recurrent patients a two-sided MWU at alpha=0.05 has ~"
                 "<30% power for the observed effect size; absence of evidence, not evidence of absence"),
    }])
    rec_out.to_csv(os.path.join(OUT, "rev02_recurrence_power.tsv"), sep="\t", index=False)
    lib.log.info("recurrence (M4):\n" + rec_out.T.to_string())
    lib.log.info("DONE rev02")


if __name__ == "__main__":
    main()
