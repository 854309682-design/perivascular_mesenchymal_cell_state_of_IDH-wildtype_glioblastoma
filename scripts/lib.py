#!/usr/bin/env python3
"""Shared data-loading utilities for Paper 2.

Reads the GSE103224 (Yuan et al. 2018, Genome Medicine) per-sample UMI matrices
delivered as 'genes x cells' txt.gz files. Layout per row: ENSG-ID(with version)
in column 0, HGNC symbol in column 1, then the cell UMI counts. No header row.
Barcodes are synthesized as '<sample>_col<N>'.
"""
import os
import gzip
import logging
import sys
import numpy as np
import pandas as pd
import scanpy as sc

# ----------------------------------------------------------------------------
# Paths
# Portable: override with the environment variables GLIOMA_PROJ (repo root) and
# GLIOMA_DATA (directory holding the downloaded public data). The defaults derive
# from this file's own location, so the scripts run from any clone that sits next
# to a `Datasets/` directory. See DATA_MANIFEST.md for the expected layout.
# ----------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))                  # <repo>/scripts
PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(_HERE)      # <repo>


def _pick_dir(*candidates):
    """First existing directory, else the last candidate (keeps messages readable)."""
    for cand in candidates:
        if cand and os.path.isdir(cand):
            return cand
    return candidates[-1]


DATASETS = os.environ.get("GLIOMA_DATA") or _pick_dir(
    os.path.join(os.path.dirname(PROJ), "Datasets"),                # authors' layout
    os.path.join(os.path.dirname(os.path.dirname(PROJ)), "Datasets"),  # <project>/<paper>/<release>
    os.path.join(PROJ, "Datasets"),
    os.path.join(PROJ, "data"),
)
DATA = DATASETS                                        # alias used by the multi-cohort scripts
RAW_DIR = os.path.join(DATASETS, "raw")                # TCGA STAR counts, clinical XML, GENCODE
PROC_DIR = os.path.join(DATASETS, "processed")         # GSE103224 h5ad objects
COHORTS = os.path.join(DATASETS, "cohorts")            # per-cohort extracted matrices
GSE103224_DIR = os.path.join(DATASETS, "GSE103224")
GSE131928_DIR = os.path.join(DATASETS, "GSE131928")
NEFTEL_TPM_SS2 = os.path.join(GSE131928_DIR, "GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
NEFTEL_TPM_10X = os.path.join(GSE131928_DIR, "GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz")
GSE70630_TXT = os.path.join(DATASETS, "GSE70630_OG_processed_data_v2.txt.gz")
GSE89567_TXT = os.path.join(DATASETS, "GSE89567_IDH_A_processed_data.txt.gz")

GSE103224_SAMPLES = ["PJ016", "PJ017", "PJ018", "PJ025", "PJ030",
                     "PJ032", "PJ035", "PJ048"]


def setup_logging(name="paper2", path=None):
    handlers = [logging.StreamHandler(sys.stdout)]
    if path:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        handlers.append(logging.FileHandler(path))
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=handlers)


log = logging.getLogger("paper2")


def gse103224_matrix_path(sample):
    d = GSE103224_DIR
    for f in os.listdir(d):
        if f"{sample}.filtered.matrix.txt.gz" in f:
            return os.path.join(d, f)
    raise FileNotFoundError(f"No matrix for {sample} in {d}")


def read_sample_matrix(sample):
    """Return (counts_df[gene_id x cells], symbol_Series[gene_id])."""
    path = gse103224_matrix_path(sample)
    log.info(f"Reading {sample} from {os.path.basename(path)}")
    df = pd.read_csv(path, sep="\t", header=None, index_col=0)
    gene_id = df.index.map(lambda s: s.split(".")[0])
    symbol = df.iloc[:, 0].astype(str)
    counts = df.iloc[:, 1:].copy()
    counts.index = pd.Index(gene_id, name="gene_id")
    symbol.index = counts.index
    counts.columns = pd.Index([f"{sample}_col{i}" for i in range(counts.shape[1])],
                              name="barcode")
    return counts, symbol


def dedup_by_symbol(counts, symbol):
    """Collapse gene_id rows to one row per gene symbol, keeping the row with the
    highest mean UMI. Returns counts indexed by symbol + a gene_id map Series."""
    gene_id = counts.index.to_numpy()
    sym = symbol.reindex(counts.index).to_numpy()
    # mean per gene row
    m = counts.mean(axis=1).to_numpy()
    df_inner = pd.DataFrame({"gene_id": gene_id, "symbol": sym, "mean": m}, index=counts.index)
    # keep one row per symbol (highest mean)
    keep = df_inner.groupby("symbol")["mean"].idxmax()
    c = counts.loc[keep].copy()
    c.index = pd.Index(df_inner.loc[keep, "symbol"].values, name="gene")
    idmap = pd.Series(df_inner.loc[keep, "gene_id"].values, index=c.index, name="gene_id")
    return c, idmap


def load_gse103224(samples=None):
    """Concatenate all GSE103224 samples into one df (genes x cells) indexed by symbol."""
    samples = samples or GSE103224_SAMPLES
    # Use symbol map from the first sample as canonical (genes rows are identical across samples)
    counts_all, idmap = None, None
    symbol_series = None
    for i, s in enumerate(samples):
        counts, symbol = read_sample_matrix(s)
        if i == 0:
            counts_all = counts
            symbol_series = symbol
        else:
            counts_all = pd.concat([counts_all, counts], axis=1)
    c, idmap = dedup_by_symbol(counts_all, symbol_series)
    log.info(f"Concatenated {len(samples)} samples: {c.shape[0]} genes x {c.shape[1]} cells")
    return c, idmap


def to_anndata(df, idmap=None):
    """Convert a genes x cells DataFrame (index=gene) to AnnData with obs['sample']."""
    adata = sc.AnnData(df.T.astype(np.float32))
    adata.obs["barcode"] = df.columns.to_numpy()
    adata.obs["sample"] = adata.obs["barcode"].map(lambda b: b.rsplit("_", 1)[0]).to_numpy()
    adata.var["gene"] = df.index.to_numpy()
    if idmap is not None:
        adata.var["gene_id"] = idmap.reindex(df.index).to_numpy()
    adata.var_names = adata.var["gene"].values
    adata.var_names_make_unique()
    return adata
