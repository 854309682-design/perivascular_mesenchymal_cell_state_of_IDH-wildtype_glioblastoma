#!/usr/bin/env python3
"""rev18 — Supplementary Figure S1 rebuilt with a legible, full-width UMAP (review item 4.1).

The submitted S1c packed the integration UMAP into a small panel whose labels were illegible
at any print size. This rebuild keeps the same content — QC distributions, the integration
UMAP, the malignant/non-malignant split, per-sample mixing and the inferCNV support — but
gives the UMAP a full-width row with readable points and an outside legend, and states the
cell funnel (23,793 raw -> 23,785 after thresholds -> 42 scrublet doublets removed ->
23,743 retained; results/step01_qc/qc_summary.tsv, logs/step01.log).

Outputs: figures/FigureS1_supp.{pdf,jpeg} (JDEG + vector PDF, no grid lines).
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__))
import lib

FIG = os.path.join(lib.PROJ, "figures")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": False, "font.family": "DejaVu Sans"})
PAL = plt.get_cmap("tab10")
SAMPLES = ["PJ016", "PJ017", "PJ018", "PJ025", "PJ030", "PJ032", "PJ035", "PJ048"]


def main():
    lib.setup_logging()
    a = sc.read_h5ad(os.path.join(lib.PROC_DIR, "gse103224_integrated.h5ad"))
    um = np.asarray(a.obsm["X_umap"])
    # clip axes to a robust range so outlier cells cannot squash the map
    lo, hi = np.percentile(um, [0.3, 99.7], axis=0)
    o = a.obs

    fig = plt.figure(figsize=(13, 11))
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1.15, 1.15, 1.15], hspace=0.42, wspace=0.28)

    # ---- QC distributions (post-filter), thresholds marked ----
    for j, (col, lab, thr) in enumerate([("nCount", "UMI counts per cell", 500),
                                         ("nGene", "genes per cell", 300),
                                         ("pctMT", "mitochondrial %", 15)]):
        ax = fig.add_subplot(gs[0, j])
        data = [o.loc[o["sample"] == s, col].values for s in SAMPLES]
        ax.violinplot(data, showextrema=False, widths=0.85)
        ax.set_xticks(range(1, len(SAMPLES) + 1))
        ax.set_xticklabels([s[-3:] for s in SAMPLES], rotation=0, fontsize=7)
        if col != "pctMT":
            ax.set_yscale("log")
        ax.axhline(thr, color="#C0392B", ls="--", lw=0.9)
        ax.text(0.02, thr, f" threshold {thr}", color="#C0392B", fontsize=6.5, va="bottom",
                transform=ax.get_yaxis_transform())
        ax.set_ylabel(lab)
        ax.set_title(f"({chr(97+j)}) {lab} after filtering", fontsize=9, loc="left")
    fig.text(0.01, 0.985, "23,793 cells raw  →  23,785 kept by thresholds  →  42 scrublet "
             "doublets removed  →  23,743 cells retained", fontsize=10)

    # ---- full-width UMAP by cell type ----
    ax = fig.add_subplot(gs[1, :])
    mal = (o["broad_celltype"] == "malignant_candidate").values
    ax.scatter(um[mal, 0], um[mal, 1], s=2.0, c="#D5DBDB", linewidths=0, label="malignant candidate")
    nonmal = ~mal
    ax.scatter(um[nonmal, 0], um[nonmal, 1], s=3.5, c="#C0392B", linewidths=0, label="non-malignant")
    for i, cl in enumerate(["11", "16", "2"]):
        sel = (o["leiden"].astype(str) == cl).values
        if sel.sum():
            ax.annotate(f"cluster {cl} (excluded, n={int(sel.sum())})",
                        (um[sel, 0].mean(), um[sel, 1].mean()), fontsize=8, color="#7B241C",
                        ha="center")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, fontsize=9, loc="upper right", markerscale=4)
    ax.set_title("(d) Harmony-integrated UMAP of 23,743 cells: 20,643 malignant candidates "
                 "(grey) and 3,100 non-malignant cells in clusters 11/16/2 (red)", fontsize=9, loc="left")

    # ---- UMAP by sample ----
    ax = fig.add_subplot(gs[2, 0])
    for i, s in enumerate(SAMPLES):
        sel = (o["sample"] == s).values
        ax.scatter(um[sel, 0], um[sel, 1], s=1.4, color=PAL(i), linewidths=0, label=s)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, fontsize=7, markerscale=4, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.22))
    ax.set_title("(e) Batch mixing: cells by patient", fontsize=9, loc="left")

    # ---- inferCNV support ----
    ax = fig.add_subplot(gs[2, 1])
    order = o.groupby("leiden")["cnv_score"].mean().sort_values().index.astype(str)
    data = [o.loc[o["leiden"].astype(str) == cl, "cnv_score"].values for cl in order]
    ax.boxplot(data, tick_labels=order, showfliers=False, widths=0.7)
    ax.axhline(0, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("Leiden cluster"); ax.set_ylabel("inferCNV score")
    plt.setp(ax.get_xticklabels(), fontsize=6, rotation=90)
    ax.set_title("(f) Genome-wide CNV score by cluster", fontsize=9, loc="left")

    # ---- cell composition ----
    ax = fig.add_subplot(gs[2, 2])
    comp = pd.crosstab(o["sample"], o["broad_celltype"], normalize="index")
    bottom = np.zeros(len(comp))
    for i, c in enumerate(comp.columns):
        ax.bar(range(len(comp)), comp[c].values, bottom=bottom, color=["#D5DBDB", "#C0392B"][i],
               label=c)
        bottom += comp[c].values
    ax.set_xticks(range(len(comp)))
    ax.set_xticklabels([s[-3:] for s in comp.index], fontsize=7)
    ax.set_ylabel("fraction of cells")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    ax.set_title("(g) Malignant/non-malignant split per tumour", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[3, :]); ax.axis("off")
    ax.text(0.0, 0.9, "(h) Malignant-cell identification. Cells were filtered at nUMI ≥ 500, "
            "nGene ≥ 300 and mitochondrial fraction < 15%, and 42 scrublet doublets were removed, "
            "leaving 23,743 cells. Harmony integration (batch = sample) with Leiden clustering "
            "gave 18 clusters; canonical marker scores and genome-wide inferCNV profiles "
            "identified clusters 11, 16 and 2 as non-malignant (immune, endothelial, "
            "oligodendrocyte), leaving 20,643 malignant cells for state analysis. Panel (d) is "
            "the full-width integration UMAP: the malignant compartment spans the map and the "
            "excluded non-malignant clusters form discrete islands.", fontsize=9, va="top")

    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"FigureS1_supp.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote FigureS1_supp.pdf/.jpeg")
    print("cells:", a.n_obs, "| malignant candidates:", int(mal.sum()),
          "| non-malignant:", int(nonmal.sum()))


if __name__ == "__main__":
    main()
