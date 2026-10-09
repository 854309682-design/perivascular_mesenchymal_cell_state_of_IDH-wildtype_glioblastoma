#!/usr/bin/env python3
"""Step 3 — Four-state meta-module scoring + state composition (single cohort).

Scores the Neftel 2019 meta-modules (MES1/MES2/AC/OPC/NPC1/NPC2) on the
GSE103224 malignant cells (full-gene normalized expression), assigns each cell to
its highest-scoring state (MES/NPC/AC/OPC), and quantifies per-sample state
composition, MES-like proportion and state entropy (= an index of plasticity).

NOTE (single-cohort scope): the primary-vs-recurrence and IDH-wt/mut *shifts*
(H1) require the other cohorts (Neftel GSE131928, Abdelfattah GSE182109) which
cannot be downloaded in this session (network throttle). This step therefore
produces the state-composition baseline + MES-enriched cell subset used by
Steps 4-5, and the shift test is deferred to the multi-cohort scaffold.

Outputs:
  results/step03_state/state_scores_per_cell.tsv
  results/step03_state/state_composition_by_sample.tsv
  results/step03_state/state_composition.png
  results/step03_state/state_umap.png
  data/processed/gse103224_states.h5ad
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, fontsize=8, frameon=False, facecolor="white")

QC = os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad")
INTG = os.path.join(lib.PROC_DIR, "gse103224_integrated.h5ad")
OUT = os.path.join(lib.PROJ, "results", "step03_state")
os.makedirs(OUT, exist_ok=True)
SEED = 42


def shannon_entropy(probs):
    p = np.asarray(probs, dtype=float)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if p.size else 0.0


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step03.log"))
    modules = load_neftel_modules()
    state_genes, sub_genes = four_state_tables(modules)

    qc = sc.read_h5ad(QC)
    intg = sc.read_h5ad(INTG)
    lib.log.info(f"QC: {qc.n_obs} x {qc.n_vars}; integrated: {intg.n_obs} x {intg.n_vars}")

    # ---- full-gene normalization of QC object ----
    if not sp.issparse(qc.X):
        qc.X = sp.csr_matrix(qc.X.astype(np.float32))
    sc.pp.filter_genes(qc, min_cells=3)
    sc.pp.normalize_total(qc, target_sum=1e4)
    sc.pp.log1p(qc)
    lib.log.info(f"Normalized full-gene object: {qc.n_obs} x {qc.n_vars}")

    # ---- align malignant flag + sample from integrated object by barcode ----
    meta = pd.DataFrame({"sample": intg.obs["sample"].values,
                         "malignant": intg.obs["malignant"].values.astype(int)},
                        index=intg.obs_names)
    qc.obs["malignant"] = meta.loc[qc.obs_names, "malignant"].values
    qc.obs["sample"] = meta.loc[qc.obs_names, "sample"].values

    # ---- score modules (per submodule) ----
    for name, genes in sub_genes.items():
        present = [g for g in genes if g in qc.var_names]
        if len(present) < 5:
            lib.log.warning(f"module {name}: only {len(present)}/ genes present; skipping")
        sc.tl.score_genes(qc, present, score_name=f"state_{name}", ctrl_size=50)
        lib.log.info(f"scored state_{name} with {len(present)} genes")

    # four-state scores (MES = max of MES1/MES2, NPC = max of NPC1/NPC2)
    qc.obs["score_MES"] = np.maximum(qc.obs["state_MES1"].values, qc.obs["state_MES2"].values)
    qc.obs["score_NPC"] = np.maximum(qc.obs["state_NPC1"].values, qc.obs["state_NPC2"].values)
    qc.obs["score_AC"] = qc.obs["state_AC"].values
    qc.obs["score_OPC"] = qc.obs["state_OPC"].values

    # state label = argmax over the four states
    scores = qc.obs[["score_MES", "score_NPC", "score_AC", "score_OPC"]].values
    idx = np.argmax(scores, axis=1)
    qc.obs["state"] = [STATE_NAMES[i] for i in idx]
    qc.obs["state_maxscore"] = scores[np.arange(len(idx)), idx]

    # ---- analyse on malignant cells only ----
    mal = qc[qc.obs["malignant"] == 1].copy()
    lib.log.info(f"Malignant cells: {mal.n_obs}")

    comp = mal.obs.groupby(["sample", "state"]).size().unstack(fill_value=0)
    comp_norm = comp.div(comp.sum(axis=1), axis=0)
    # MES proportion + MES1/MES2 split + entropy
    comp_norm["MES_prop"] = comp_norm.get("MES", 0)
    comp_norm["NPC_prop"] = comp_norm.get("NPC", 0)
    # MES1/MES2 split among MES cells
    mes1 = (mal.obs["state_MES1"] >= mal.obs["state_MES2"]).groupby(mal.obs["sample"]).mean()
    mes1 = mes1.reindex(comp_norm.index).fillna(0)
    comp_norm["MES1_fraction_of_MES"] = mes1
    comp_norm["entropy"] = comp_norm[STATE_NAMES].apply(shannon_entropy, axis=1)
    comp_norm.index.name = "sample"
    comp_norm.to_csv(os.path.join(OUT, "state_composition_by_sample.tsv"), sep="\t")

    # per-cell state scores
    cell = mal.obs[["sample", "state",
                    "score_MES", "score_NPC", "score_AC", "score_OPC",
                    "state_MES1", "state_MES2", "state_NPC1", "state_NPC2",
                    "state_maxscore"]].copy()
    cell.to_csv(os.path.join(OUT, "state_scores_per_cell.tsv"), sep="\t")
    lib.log.info(f"\n{comp_norm[['MES_prop','NPC_prop','entropy']].round(3).to_string()}")

    # overall state composition
    overall = mal.obs["state"].value_counts(normalize=True).round(3)
    lib.log.info("Overall state composition:\n" + overall.to_string())

    # ---- figures ----
    # map state labels onto the integrated object (has the UMAP embedding)
    state_by_bc = dict(zip(mal.obs_names, mal.obs["state"].values))
    intg.obs["state"] = [state_by_bc.get(b, "not_malignant") for b in intg.obs_names]
    intg.obs["state"] = intg.obs["state"].astype(str)
    lib.log.info("state on integrated:", intg.obs["state"].value_counts().to_dict())
    try:
        fig, ax = plt.subplots(figsize=(8, 4))
        comp_norm[STATE_NAMES].plot(kind="bar", stacked=True, ax=ax,
                                    color=["#C0392B", "#2980B9", "#27AE60", "#8E44AD"])
        ax.set_ylabel("proportion of malignant cells"); ax.set_xlabel("sample")
        ax.set_title("Cell-state composition (Neftel meta-modules)")
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "state_composition.png"), dpi=150)
        plt.close("all")

        with plt.rc_context({"figure.figsize": (8, 6)}):
            sc.pl.umap(intg, color="state", show=False, palette="tab10", frameon=False)
            plt.savefig(os.path.join(OUT, "state_umap.png"), dpi=150, bbox_inches="tight")
            plt.close("all")
    except Exception as e:
        lib.log.warning(f"figure error: {e}")

    intg.write_h5ad(os.path.join(lib.PROC_DIR, "gse103224_states.h5ad"))
    lib.log.info("Saved gse103224_states.h5ad")


if __name__ == "__main__":
    main()
