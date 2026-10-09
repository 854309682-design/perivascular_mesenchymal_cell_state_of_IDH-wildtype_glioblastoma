#!/usr/bin/env python3
"""Step 1 — QC & standardization for GSE103224 (Paper 2).

Reads data/processed/gse103224_raw.h5ad, filters to protein-coding + lncRNA
genes (GENCODE v44 basic annotation), computes QC metrics, applies standard 10x
thresholds (nCount>=500, nGene>=300, MT%%<15), removes doublets with scrublet
(per sample), records per-sample QC summary, and saves the QC'd object.

Outputs:
  results/step01_qc/qc_summary.tsv
  results/step01_qc/qc_violin_<sample>.png
  results/step01_qc/qc_violin_all.png
  data/processed/gse103224_qc.h5ad
"""
import os
import sys
import re
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, fontsize=8, frameon=False, facecolor="white")

RAW = os.path.join(lib.PROC_DIR, "gse103224_raw.h5ad")
OUT = os.path.join(lib.PROJ, "results", "step01_qc")
os.makedirs(OUT, exist_ok=True)

GENCODE_DIR = os.path.join(lib.RAW_DIR, "gencode")
GENCODE_GTF = os.path.join(GENCODE_DIR, "gencode.v44.basic.annotation.gtf.gz")

# QC thresholds
MIN_COUNTS = 500
MIN_GENES = 300
MAX_PCT_MT = 15.0
RANDOM_SEED = 42


def parse_gencode_gtf(path, gene_types=None):
    """Return {gene_id: gene_type} for 'gene' features. Some genes have multiple
    gene_type values (rare); take the first non-'gene' type."""
    keep = {"protein_coding", "lncRNA"} if gene_types is None else set(gene_types)
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                f = line.split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attr = f[8]
                m = re.search(r'gene_id "([^"]+)"', attr)
                mt = re.search(r'gene_type "([^"]+)"', attr)
                gb = re.search(r'gene_biotype "([^"]+)"', attr)
                if not m:
                    continue
                gtype = (mt.group(1) if mt else (gb.group(1) if gb else None))
                gid = m.group(1).split(".")[0]
                if gtype is None:
                    continue
                # store, prefer protein_coding/lncRNA if duplicated
                if gid not in out or (gtype in keep and out[gid] not in keep):
                    out[gid] = gtype
    return out


def main():
    lib.setup_logging(path=os.path.join(OUT, "..", "step01.log"))
    adata = sc.read_h5ad(RAW)
    lib.log.info(f"Loaded raw: {adata.n_obs} cells x {adata.n_vars} genes")

    # ---- gene-type filter (protein-coding + lncRNA) via GENCODE v44 ----
    if os.path.exists(GENCODE_GTF):
        try:
            lib.log.info("Parsing GENCODE v44 basic annotation for gene biotype...")
            gtype = parse_gencode_gtf(GENCODE_GTF, {"protein_coding", "lncRNA"})
            # map by gene_id (var.gene_id); fall back to symbol for missing
            gid = adata.var["gene_id"].values
            mapped = np.array([gtype.get(g, "") for g in gid])
            adata.var["gene_type"] = mapped
            keep_genes = (adata.var["gene_type"].isin(["protein_coding", "lncRNA"])).values
            adata = adata[:, keep_genes].copy()
            lib.log.info(
                f"Kept {adata.n_vars} protein-coding+lncRNA genes "
                f"({keep_genes.sum()}); {np.sum(~keep_genes)} removed")
        except (EOFError, OSError, ValueError) as e:
            lib.log.warning(f"GENCODE gtf unusable/partial ({e}); skipping gene-type filter")
    else:
        lib.log.warning("GENCODE gtf not found; skipping gene-type filter")

    # ---- QC metrics ----
    adata.var["mt"] = adata.var_names.str.startswith("MT-").astype(int)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None,
                               log1p=False)
    adata.obs["nCount"] = adata.obs["total_counts"]
    adata.obs["nGene"] = adata.obs["n_genes_by_counts"]
    adata.obs["pctMT"] = adata.obs["pct_counts_mt"]

    # record pre-filter counts
    pre = adata.obs.groupby("sample").agg(n_cells=("sample", "size"),
                                          median_counts=("nCount", "median"),
                                          median_genes=("nGene", "median"),
                                          median_pctMT=("pctMT", "median"))

    # ---- QC violin (all samples, pre-filter) ----
    sc.pl.violin(adata, ["nCount", "nGene", "pctMT"], groupby="sample", rotation=90,
                 multi_panel=True, ax=None, show=False)
    plt.savefig(os.path.join(OUT, "qc_violin_all.png"), bbox_inches="tight")
    plt.close("all")

    # ---- apply thresholds ----
    lib.log.info("Applying QC thresholds: nCount>=%d, nGene>=%d, pctMT<%.1f%%",
                 MIN_COUNTS, MIN_GENES, MAX_PCT_MT)
    keep = (adata.obs["nCount"] >= MIN_COUNTS) & (adata.obs["nGene"] >= MIN_GENES) & \
           (adata.obs["pctMT"] < MAX_PCT_MT)
    n_before = adata.n_obs
    adata = adata[keep].copy()
    lib.log.info(f"Threshold filter kept {adata.n_obs}/{n_before} cells "
                 f"({adata.n_obs / n_before * 100:.1f}%)")

    # ---- scrublet doublet removal (per sample) ----
    doublet_counts = {}
    for sample in adata.obs["sample"].unique():
        sub = adata[adata.obs["sample"] == sample]
        if sub.n_obs < 200:
            lib.log.info(f"{sample}: only {sub.n_obs} cells, skipping scrublet")
            doublet_counts[sample] = np.zeros(sub.n_obs, dtype=bool)
            continue
        lib.log.info(f"Scrublet on {sample} ({sub.n_obs} cells)")
        sc.pp.scrublet(sub, random_state=RANDOM_SEED)
        doublet_counts[sample] = sub.obs["predicted_doublet"].values
    doublet = np.concatenate([doublet_counts[s] for s in adata.obs["sample"].unique()]) \
        if adata.obs["sample"].unique().size else np.array([], dtype=bool)
    adata.obs["predicted_doublet"] = doublet
    n_doublet = int(doublet.sum())
    adata = adata[~doublet].copy()
    lib.log.info(f"Scrublet removed {n_doublet} doublets; remaining {adata.n_obs}")

    # ---- QC summary table ----
    post = adata.obs.groupby("sample").agg(n_cells=("sample", "size"),
                                           median_counts=("nCount", "median"),
                                           median_genes=("nGene", "median"),
                                           median_pctMT=("pctMT", "median"))
    summary = pre.join(post, rsuffix="_post", how="left")
    summary["filter_rate_pct"] = (1 - summary["n_cells_post"] / summary["n_cells"]) * 100
    summary.to_csv(os.path.join(OUT, "qc_summary.tsv"), sep="\t")
    lib.log.info("\n" + summary.to_string())

    # per-sample violin of surviving cells
    sc.pl.violin(adata, ["nCount", "nGene", "pctMT"], groupby="sample", rotation=90,
                 multi_panel=True, show=False)
    plt.savefig(os.path.join(OUT, "qc_violin_post.png"), bbox_inches="tight")
    plt.close("all")

    adata.write_h5ad(os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad"))
    lib.log.info(f"Saved QC'd object: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
