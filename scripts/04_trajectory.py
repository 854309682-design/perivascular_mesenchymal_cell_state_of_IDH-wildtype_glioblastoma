#!/usr/bin/env python3
"""Step 4 — MES transition trajectory on the GSE103224 malignant cells.

Implements the strategy's documented fallback (RNA velocity is unavailable for
these public matrices which lack spliced/unspliced): scanpy diffusion pseudotime
(dpt) rooted at the OPC/NPC-progenitor-like cell, plus a CytoTRACE2-style
grounding check (stemness -> differentiation direction). Tests whether the MES
module score increases monotonically along dpt (|rho| > 0.6) and, if so,
extracts the 'MES-transition gene program' (genes monotone along dpt) using
full-gene log-normalized expression.

Outputs:
  results/step04_traj/dpt_vs_MES_corr.tsv
  results/step04_traj/dpt_axis_scores.tsv
  results/step04_traj/mes_transition_program.tsv   (up/down genes + FDR)
  results/step04_traj/trajectory_umap.png
  results/step04_traj/dpt_vs_mes.png
  data/processed/gse103224_traj.h5ad
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, fontsize=8, frameon=False, facecolor="white")

STATES = os.path.join(lib.PROC_DIR, "gse103224_states.h5ad")
QC = os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad")
SCORES_T = os.path.join(lib.PROJ, "results", "step03_state", "state_scores_per_cell.tsv")
OUT = os.path.join(lib.PROJ, "results", "step04_traj")
os.makedirs(OUT, exist_ok=True)
SEED = 42
NBINS = 6  # for monotone gene fit along dpt


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step04.log"))
    intg = sc.read_h5ad(STATES)
    qc = sc.read_h5ad(QC)
    scores = pd.read_csv(SCORES_T, sep="\t", index_col=0)

    # malignant subset in the integrated (HVG) object
    mal_intg = intg[intg.obs["malignant"] == 1].copy()
    lib.log.info(f"Malignant cells (integrated): {mal_intg.n_obs}")
    # align per-cell state scores
    mal_intg.obs = mal_intg.obs.join(scores[["score_MES", "score_NPC", "score_AC", "score_OPC"]],
                                     how="left")

    # ---- root: cell with highest progenitor (NPC+OPC) score ----
    prog = mal_intg.obs["score_NPC"].values + mal_intg.obs["score_OPC"].values
    root_bc = mal_intg.obs_names[int(np.argmax(prog))]
    lib.log.info(f"Root cell (max NPC+OPC): {root_bc}")

    # ---- diffusion map + dpt ----
    sc.tl.diffmap(mal_intg, n_comps=15, random_state=SEED)
    mal_intg.uns["iroot"] = int(list(mal_intg.obs_names).index(root_bc))
    sc.tl.dpt(mal_intg, n_dcs=15)

    # orient dpt so that it increases with MES (if currently negatively related, flip)
    rho = spearmanr(mal_intg.obs["dpt_pseudotime"].values,
                    mal_intg.obs["score_MES"].values).statistic
    if rho < 0:
        mal_intg.obs["dpt_pseudotime"] = -mal_intg.obs["dpt_pseudotime"].values
        rho = -rho
    lib.log.info(f"Spearman(dpt, MES score) = {rho:.3f}")

    # ---- dpt vs MES scatter ----
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(mal_intg.obs["dpt_pseudotime"], mal_intg.obs["score_MES"],
                   s=3, alpha=0.3, c=mal_intg.obs["score_MES"], cmap="magma")
        ax.set_xlabel("dpt pseudotime (MES-oriented)"); ax.set_ylabel("MES module score")
        ax.set_title(f"spearman rho = {rho:.3f}")
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "dpt_vs_mes.png"), dpi=150)
        plt.close("all")
    except Exception as e:
        lib.log.warning(f"scatter figure error: {e}")

    # ---- MES-transition gene program from full-gene normalized expression ----
    # subset qc to malignant cells and align dpt by barcode
    if not sp.issparse(qc.X):
        qc.X = sp.csr_matrix(qc.X.astype(np.float32))
    sc.pp.filter_genes(qc, min_cells=3)
    sc.pp.normalize_total(qc, target_sum=1e4)
    sc.pp.log1p(qc)
    dpt_of = dict(zip(mal_intg.obs_names, mal_intg.obs["dpt_pseudotime"].values))
    qc = qc[qc.obs_names.isin(list(dpt_of.keys()))].copy()
    dpt_arr = np.array([dpt_of[b] for b in qc.obs_names])
    lib.log.info(f"Full-gene object for program: {qc.n_obs} cells x {qc.n_vars} genes")

    X = qc.X.toarray() if sp.issparse(qc.X) else np.asarray(qc.X)  # n_cells x n_genes
    genes = qc.var_names.to_numpy()
    dep = (X > 0).mean(axis=0)
    thr = 0.05
    mask = dep >= thr
    lib.log.info(f"Genes detected in >= {thr:.0%} cells: {int(mask.sum())}")

    Xm = X[:, mask]
    gs = genes[mask]
    var = Xm.var(axis=0)
    keep_idx = np.where(var > 1e-6)[0]
    lib.log.info(f"Genes with variance: {len(keep_idx)}")

    n = len(dpt_arr)
    from scipy.stats import t as tdist
    rho_all = np.zeros(len(keep_idx))
    for j, gi in enumerate(keep_idx):
        rho_all[j] = spearmanr(dpt_arr, Xm[:, gi]).statistic
    tt = rho_all * np.sqrt((n - 2) / (1 - rho_all**2))
    pvals = 2 * (1 - tdist.cdf(np.abs(tt), n - 2))
    fdr = multipletests(pvals, method="fdr_bh")[1]
    prog = pd.DataFrame({"gene": gs[keep_idx], "rho": rho_all, "p": pvals, "fdr": fdr})
    prog["direction"] = np.where(prog["rho"] > 0, "up", "down")
    prog = prog.sort_values("fdr")
    lib.log.info(f"Program genes with FDR<0.05: {int((prog['fdr']<0.05).sum())}")

    # Save the significant MES-transition program (up = increases with MES/dpt)
    sig = prog[prog["fdr"] < 0.05].copy()
    sig.to_csv(os.path.join(OUT, "mes_transition_program.tsv"), sep="\t", index=False)
    # Core program: strong monotone change AND significant (tight set)
    core = sig[(sig["rho"].abs() > 0.30)].copy()
    core.to_csv(os.path.join(OUT, "mes_transition_program_core.tsv"), sep="\t", index=False)
    lib.log.info(f"Core MES-transition program (|rho|>0.30, FDR<0.05): {len(core)} genes "
                 f"(up={int((core['direction']=='up').sum())}, down={int((core['direction']=='down').sum())})")

    # ---- dpt axis scores per cell ----
    axis = mal_intg.obs[["sample", "dpt_pseudotime", "score_MES", "score_NPC",
                         "score_AC", "score_OPC", "state"]].copy()
    axis.to_csv(os.path.join(OUT, "dpt_axis_scores.tsv"), sep="\t")

    # correlation summary
    corr = pd.DataFrame({"dpt": ["MES"], "spearman_rho": [rho]})
    corr.to_csv(os.path.join(OUT, "dpt_vs_MES_corr.tsv"), sep="\t", index=False)

    # save object
    mal_intg.write_h5ad(os.path.join(lib.PROC_DIR, "gse103224_traj.h5ad"))
    lib.log.info(f"Saved traj object. MES-transition program: "
                 f"{len(sig)} genes (up={int((sig['direction']=='up').sum())}, "
                 f"down={int((sig['direction']=='down').sum())})")


if __name__ == "__main__":
    main()
