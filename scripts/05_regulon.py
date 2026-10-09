#!/usr/bin/env python3
"""Step 5 — Regulon / TF-activity inference for the MES transition (GSE103224).

pySCENIC (GRNBoost2 + cisTarget) is not runnable in this session: `pyscenic` is
not installed, the hg38 cisTarget motif-rankings database is a multi-GB download
(blocked by the network throttle), and GRNBoost2 on ~20k cells is heavy. This
step uses the strategy's alternative regulon evidence channel — a curated
TF->target prior (DoRothEA, confidence A/B/C) with a multivariate linear model
(mLM) to derive per-cell TF activity (regulon-activity analogue). It then ranks
TFs by the correlation of their activity with the MES module score across the
malignant cells, and cross-checks the top hits against literature-prior master
TFs (FOSL2/JUN, STAT3, CEBPB, TEAD1/3, NFKB1).

Outputs:
  results/step05_regulon/dorothea_net.tsv          (cached prior)
  results/step05_regulon/tf_activity_vs_MES.tsv    (per-TF rho + rank)
  results/step05_regulon/top_master_tf.tsv
  data/processed/gse103224_regulon.h5ad
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib

sc.settings.verbosity = 0
STATES = os.path.join(lib.PROC_DIR, "gse103224_states.h5ad")
QC = os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad")
OUT = os.path.join(lib.PROJ, "results", "step05_regulon")
os.makedirs(OUT, exist_ok=True)
SEED = 42

LITERATURE_TFS = ["FOSL2", "FOS", "JUN", "JUNB", "JUND", "STAT3", "CEBPB", "CEBPA",
                  "TEAD1", "TEAD3", "TEAD4", "NFKB1", "RELA", "BATF", "IRF1", "SPI1"]
DOROTHEA_LEVELS = ["A", "B", "C"]


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step05.log"))
    # ---- fetch / cache DoRothEA net ----
    net_path = os.path.join(OUT, "dorothea_net.tsv")
    if os.path.exists(net_path):
        net = pd.read_csv(net_path, sep="\t")
        lib.log.info(f"Loaded cached DoRothEA net: {net.shape}")
    else:
        import decoupler as dc
        net = dc.op.dorothea(organism="human", levels=DOROTHEA_LEVELS)
        net.to_csv(net_path, sep="\t", index=False)
        lib.log.info(f"Fetched DoRothEA net: {net.shape}")

    # ---- full-gene normalized expression of malignant cells ----
    qc = sc.read_h5ad(QC)
    if not sp.issparse(qc.X):
        qc.X = sp.csr_matrix(qc.X.astype(np.float32))
    sc.pp.filter_genes(qc, min_cells=3)
    sc.pp.normalize_total(qc, target_sum=1e4)
    sc.pp.log1p(qc)

    states = sc.read_h5ad(STATES)
    # map MES score etc. from step03 per-cell scores (indexed by barcode)
    scores = pd.read_csv(os.path.join(lib.PROJ, "results", "step03_state",
                                      "state_scores_per_cell.tsv"), sep="\t", index_col=0)
    mal = qc[qc.obs_names.isin(states.obs_names[states.obs["malignant"] == 1])].copy()
    mal.obs = mal.obs.join(scores[["score_MES", "score_NPC", "score_AC", "score_OPC"]], how="left")
    lib.log.info(f"Malignant cells for regulon: {mal.n_obs} x {mal.n_vars}")

    # ---- filter net to expressed genes + limit to TFs with enough targets ----
    net = net[net["target"].isin(mal.var_names) & net["source"].isin(mal.var_names)].copy()
    lib.log.info(f"DoRothEA net (present genes): {net.shape}")
    tgt_counts = net.groupby("source")["target"].count()
    net = net[net["source"].isin(tgt_counts[tgt_counts >= 10].index)].copy()
    lib.log.info(f"After min 10 targets per TF: {net.shape}; "
                 f"{net['source'].nunique()} TFs")

    # ---- TF activity via mLM (decoupler) ----
    import decoupler as dc
    dc.mt.decouple(mal, net, methods=["mlm"], verbose=False)
    est = mal.obsm["score_mlm"]  # cells x TFs
    if hasattr(est, "values"):
        tfs = list(est.columns)
        est_np = est.values
    else:
        est_np = est
        tfs = list(est_np.columns) if hasattr(est_np, "columns") else list(range(est_np.shape[1]))
    lib.log.info(f"mLM TF activity matrix: {est_np.shape}")

    # ---- rank TFs by correlation of activity with MES score ----
    mes = mal.obs["score_MES"].values
    rows = []
    for j, tf in enumerate(tfs):
        act = est_np[:, j]
        if np.std(act) < 1e-9:
            continue
        rho, pval = spearmanr(act, mes)
        rows.append({"tf": tf, "rho": rho, "p": pval})
    tf_df = pd.DataFrame(rows).sort_values("rho", ascending=False).reset_index(drop=True)
    from statsmodels.stats.multitest import multipletests
    tf_df["fdr"] = multipletests(tf_df["p"].values, method="fdr_bh")[1]
    tf_df.to_csv(os.path.join(OUT, "tf_activity_vs_MES.tsv"), sep="\t", index=False)
    lib.log.info(f"Ranked {len(tf_df)} TFs; top 10:\n"
                 + tf_df.head(10)[["tf", "rho", "fdr"]].to_string(index=False))

    # ---- cross-check literature master TFs ----
    lit = tf_df[tf_df["tf"].isin(LITERATURE_TFS)].sort_values("rho", ascending=False)
    lit.to_csv(os.path.join(OUT, "top_master_tf.tsv"), sep="\t", index=False)
    lib.log.info("Literature-prior TF ranks:\n" + lit[["tf", "rho", "fdr"]].to_string(index=False))

    mal.write_h5ad(os.path.join(lib.PROC_DIR, "gse103224_regulon.h5ad"))
    lib.log.info("Saved gse103224_regulon.h5ad")


if __name__ == "__main__":
    main()
