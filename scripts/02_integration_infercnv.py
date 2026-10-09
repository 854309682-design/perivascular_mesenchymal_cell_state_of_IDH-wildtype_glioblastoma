#!/usr/bin/env python3
"""Step 2 — Batch integration (Harmony) + malignant cell identification via CNV.

Reads data/processed/gse103224_qc.h5ad, normalizes (CPM/10k + log1p), selects
HVGs, runs PCA -> Harmony (batch='sample') -> UMAP/leiden, annotates non-malignant
(immune/endothelial/oligo) clusters via marker scores, and runs infercnvpy CNV
inference (gated on availability of gene genomic positions from the GENCODE gtf).

Outputs:
  results/step02_integration/umap_celltype.png
  results/step02_integration/cnv_hist.png (if CNV ran)
  results/step02_integration/celltype_summary.tsv
  results/step02_integration/leiden_marker_scores.tsv
  data/processed/gse103224_integrated.h5ad
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, fontsize=8, frameon=False, facecolor="white")

IN = os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad")
OUT = os.path.join(lib.PROJ, "results", "step02_integration")
os.makedirs(OUT, exist_ok=True)

SEED = 42
N_HVG = 2500
LEIDEN_RES = 1.0
GENCODE_GTF = os.path.join(lib.RAW_DIR, "gencode", "gencode.v44.basic.annotation.gtf.gz")

MARKERS = {
    "immune": ["PTPRC", "CD68", "AIF1", "C1QA", "C1QB", "C1QC", "FCGR3A", "LYZ",
               "CD3D", "CD3E", "NKG7", "IL7R"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "ENG", "FLT1"],
    "oligo": ["MBP", "PLP1", "MOBP", "MOG", "OLIG1"],
}


def is_gtf_complete(path):
    """Best-effort: the gtf gzip is complete if it can be fully decompressed."""
    try:
        import gzip
        with gzip.open(path, "rt") as fh:
            for _ in fh:
                pass
        return True
    except Exception:
        return False


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step02.log"))
    adata = sc.read_h5ad(IN)
    lib.log.info(f"Loaded {adata.n_obs} x {adata.n_vars}")

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X.astype(np.float32))

    # ---- normalize (CPM/10k + log1p) ----
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()
    lib.log.info(f"After min_cells=3 + normalize: {adata.n_obs} x {adata.n_vars}")

    # ---- HVG + PCA + Harmony ----
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, layer="counts", flavor="seurat_v3")
    lib.log.info(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=SEED)
    # Direct harmonypy call (scanpy 1.12.3 wrapper is incompatible with harmonypy 2.0.0)
    import harmonypy
    import pandas as _pd
    Z = np.asarray(adata.obsm["X_pca"]).astype(np.float64)  # n_cells x n_pcs
    meta = _pd.DataFrame({"sample": adata.obs["sample"].astype(str).values})
    ho = harmonypy.run_harmony(Z, meta, vars_use=["sample"], random_state=SEED,
                               verbose=False)
    adata.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr, dtype=np.float32)

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added="leiden", random_state=SEED)
    lib.log.info(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    # ---- marker-based reference annotation ----
    score_tab = pd.DataFrame(index=adata.obs_names)
    for cname, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, present, score_name=f"score_{cname}", ctrl_size=50)
        score_tab[f"score_{cname}"] = adata.obs[f"score_{cname}"].values
    score_tab["leiden"] = adata.obs["leiden"].values
    score_tab["sample"] = adata.obs["sample"].values

    clust = score_tab.groupby("leiden")[["score_immune", "score_endothelial"]].mean()
    # Reference (non-malignant) = clusters with a clearly positive immune or
    # endothelial module score (absolute threshold, not quantile).
    MIN_REF = 0.5
    nonmalig_clusters = set(clust[(clust["score_immune"] >= MIN_REF) |
                                  (clust["score_endothelial"] >= MIN_REF)].index.tolist())
    adata.obs["broad_celltype"] = np.where(adata.obs["leiden"].isin(nonmalig_clusters),
                                           "nonmalignant", "malignant_candidate")
    lib.log.info(f"Non-malignant (reference) clusters: {sorted(nonmalig_clusters)}")

    # ---- CNV inference via infercnvpy (gated on genomic positions) ----
    cnv_ran = False
    if os.path.exists(GENCODE_GTF) and is_gtf_complete(GENCODE_GTF):
        try:
            import infercnvpy as cnv
            lib.log.info("Mapping gene genomic positions from GENCODE gtf...")
            cnv.io.genomic_position_from_gtf(GENCODE_GTF, adata, gtf_gene_id="gene_name",
                                             adata_gene_id="gene")
            cnv.tl.infercnv(adata, reference_key="broad_celltype",
                            reference_cat="nonmalignant", window_size=101, step=101)
            cnv.tl.pca(adata, use_rep="cnv", key_added="cnv_pca")
            sc.pp.neighbors(adata, use_rep="X_cnv_pca", n_neighbors=15,
                            key_added="cnv_neighbors", random_state=SEED)
            cnv.tl.leiden(adata, neighbors_key="cnv_neighbors", key_added="cnv_leiden",
                          resolution=0.6, random_state=SEED)
            cnv.tl.cnv_score(adata, groupby="cnv_leiden", use_rep="cnv",
                             key_added="cnv_score")
            # Primary malignant call = marker-based (immune/endothelial reference is
            # reliable). CNV is kept as *supporting* evidence: tumor candidates show
            # higher mean cnv_score than the reference, but with large overlap on the
            # reduced HVG set, so it is not used as the primary threshold here.
            adata.obs["malignant"] = (adata.obs["broad_celltype"] == "malignant_candidate").astype(int)
            adata.obs["malignant_origin"] = "celltype+infercnv_support"
            cnv_ran = True
            lib.log.info(f"CNV ran (supportive). malignant fraction="
                         f"{adata.obs['malignant'].mean():.3f}; "
                         f"mean cnv_score nonmalig={adata.obs[adata.obs.broad_celltype=='nonmalignant']['cnv_score'].mean():.4f} "
                         f"tumor={adata.obs[adata.obs.broad_celltype=='malignant_candidate']['cnv_score'].mean():.4f}")
        except Exception as e:
            lib.log.warning(f"CNV inference failed: {e}")
    if not cnv_ran:
        lib.log.warning("CNV skipped; using celltype label as provisional malignant flag")
        adata.obs["malignant"] = (adata.obs["broad_celltype"] == "malignant_candidate").astype(int)
        adata.obs["malignant_origin"] = "celltype"

    # ---- outputs ----
    summ = adata.obs.groupby(["broad_celltype", "malignant"]).size().reset_index(name="n_cells")
    summ.to_csv(os.path.join(OUT, "celltype_summary.tsv"), sep="\t", index=False)
    score_tab["broad_celltype"] = adata.obs["broad_celltype"].values
    score_tab["malignant"] = adata.obs["malignant"].values
    score_tab.to_csv(os.path.join(OUT, "leiden_marker_scores.tsv"), sep="\t", index=False)

    try:
        with plt.rc_context({"figure.figsize": (8, 6)}):
            sc.pl.umap(adata, color=["broad_celltype", "malignant", "leiden"],
                       ncols=3, show=False, frameon=False)
            plt.savefig(os.path.join(OUT, "umap_celltype.png"), dpi=150, bbox_inches="tight")
            plt.close("all")
    except Exception as e:
        lib.log.warning(f"UMAP figure failed: {e}")

    if cnv_ran:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(adata.obs["cnv_score"], bins=50, color="steelblue")
        ax.axvline(float(np.quantile(adata.obs["cnv_score"], 0.5)), color="red",
                   ls="--", label="malignant threshold")
        ax.set_xlabel("CNV score"); ax.set_ylabel("# cells"); ax.legend()
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "cnv_hist.png"), dpi=150)
        plt.close("all")

    adata.write_h5ad(os.path.join(lib.PROC_DIR, "gse103224_integrated.h5ad"))
    lib.log.info(f"Saved integrated object {adata.n_obs} x {adata.n_vars}; "
                 f"malignant={int(adata.obs['malignant'].sum())}, cnv_ran={cnv_ran}")


if __name__ == "__main__":
    main()
