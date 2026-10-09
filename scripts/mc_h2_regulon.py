#!/usr/bin/env python3
"""H2: multi-cohort regulon consensus — which TFs' activity tracks the MES state.

Runs DoRothEA-prior TF activity (decoupler mLM) on each cohort's malignant cells
(subsampled for tractability), correlates each TF's activity with the MES score,
and reports the per-cohort rho + a consensus (mean rho across cohorts, # cohorts
with FDR<0.05). Merges with the GSE103224 result.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables
import mc_cohort_state as mcs
import mc_gse182109_h1 as g182

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "multicohort")
NET = os.path.join(lib.PROJ, "results", "step05_regulon", "dorothea_net.tsv")
MAX_CELLS = 15000
SEED = 42
LIT = ["STAT3", "NFKB1", "FOSL2", "FOSL1", "JUN", "JUNB", "CEBPB", "TEAD1", "TEAD4", "RELA", "IRF1"]


def tf_activity(ad, cohort):
    """ad must be log-normalized, already contain state_* + score_MES."""
    import decoupler as dc
    mal = ad[ad.obs["nonmalignant"] == 0].copy()
    if mal.n_obs > MAX_CELLS:
        mal = mal[np.random.RandomState(SEED).choice(mal.n_obs, MAX_CELLS, replace=False)].copy()
    # sanitize non-finite expression values
    if sp.issparse(mal.X):
        mal.X.data = np.nan_to_num(mal.X.data)
    else:
        mal.X = np.nan_to_num(mal.X)
    net = pd.read_csv(NET, sep="\t")
    net = net[net["target"].isin(mal.var_names) & net["source"].isin(mal.var_names)].copy()
    tgt = net.groupby("source")["target"].count()
    net = net[net["source"].isin(tgt[tgt >= 10].index)].copy()
    dc.mt.decouple(mal, net, methods=["mlm"], verbose=False)
    est = mal.obsm["score_mlm"]
    est_np = est.values if hasattr(est, "values") else est
    tfs = list(est.columns) if hasattr(est, "columns") else list(range(est.shape[1]))
    mes = mal.obs["score_MES"].values
    rows = []
    for j, tf in enumerate(tfs):
        act = est_np[:, j]
        if np.std(act) < 1e-9:
            continue
        rho, p = spearmanr(act, mes)
        rows.append({"tf": tf, "cohort": cohort, "rho": rho, "p": p})
    df = pd.DataFrame(rows)
    df["fdr"] = multipletests(df["p"].values, method="fdr_bh")[1]
    df = df.sort_values("rho", ascending=False)
    df.to_csv(os.path.join(OUT, f"h2_tf_{cohort}.tsv"), sep="\t", index=False)
    lib.log.info(f"[{cohort}] {len(df)} TFs; top: {df.head(3)[['tf','rho']].to_dict('records')}")
    return df


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "mc_h2.log"))
    state_sub = four_state_tables(load_neftel_modules())[1]

    dfs = {}
    # Neftel (IDH-wt GBM)
    a1 = mcs.load_log_matrix(lib.NEFTEL_TPM_SS2,
                             "Neftel", lambda c: c.split("-")[0], "Smartseq2")
    a2 = mcs.load_log_matrix(lib.NEFTEL_TPM_10X,
                             "Neftel", lambda c: c.split("_")[0], "10x")
    neft = sc.concat([a1, a2], join="outer"); neft.var_names_make_unique()
    neft = mcs.exclude_immune(neft); neft = mcs.score_state(neft, state_sub)
    dfs["Neftel"] = tf_activity(neft, "Neftel")

    # Venteicher (IDH-mut control)
    vent = mcs.load_log_matrix(lib.GSE70630_TXT,
                               "Venteicher", lambda c: c.split("_")[0], "processed")
    vent = mcs.exclude_immune(vent); vent = mcs.score_state(vent, state_sub)
    dfs["Venteicher"] = tf_activity(vent, "Venteicher")

    # GSE182109 (sample cells, load mtx samples)
    import glob
    samples = sorted(glob.glob(os.path.join(lib.COHORTS, "gse182109", "*_matrix.mtx.gz")))
    adatas = [g182.load_one(f) for f in samples]
    ad = sc.concat(adatas, join="outer"); ad.var_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    ad = mcs.exclude_immune(ad); ad = mcs.score_state(ad, state_sub)
    dfs["GSE182109"] = tf_activity(ad, "GSE182109")

    # merge with GSE103224 (from step05 tf_activity_vs_MES.tsv)
    g103 = pd.read_csv(os.path.join(lib.PROJ, "results", "step05_regulon", "tf_activity_vs_MES.tsv"), sep="\t")
    g103 = g103[["tf", "rho", "fdr"]]; g103["cohort"] = "GSE103224"
    dfs["GSE103224"] = g103

    allc = pd.concat(dfs.values(), ignore_index=True)
    # consensus: mean rho across cohorts that scored the TF, + # cohorts with FDR<0.05
    cons = allc.groupby("tf").agg(cohorts=("cohort", "nunique"),
                                  mean_rho=("rho", "mean"),
                                  n_fdr005=("fdr", lambda s: int((s < 0.05).sum()))).reset_index()
    cons = cons.sort_values("mean_rho", ascending=False)
    cons.to_csv(os.path.join(OUT, "h2_regulon_consensus.tsv"), sep="\t", index=False)
    lib.log.info("Consensus top master-TF candidates:\n" + cons.head(15)[["tf", "cohorts", "mean_rho", "n_fdr005"]].to_string(index=False))
    lit = cons[cons["tf"].isin(LIT)].sort_values("mean_rho", ascending=False)
    lib.log.info("Literature-prior TFs:\n" + lit[["tf", "cohorts", "mean_rho", "n_fdr005"]].to_string(index=False))


if __name__ == "__main__":
    main()
