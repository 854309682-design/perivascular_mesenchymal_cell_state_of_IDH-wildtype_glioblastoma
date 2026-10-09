#!/usr/bin/env python3
"""Step 0 — build the concatenated GSE103224 raw AnnData object.

The numbered pipeline starts from `<Datasets>/processed/gse103224_raw.h5ad`
(<GLIOMA_DATA>/processed/gse103224_raw.h5ad). This bootstrap creates that file
from the per-sample UMI matrices (genes x cells, 'txt.gz') downloaded from
GSE103224 — see DATA_MANIFEST.md. It composes the existing helpers
`lib.load_gse103224()` + `lib.to_anndata()` and is idempotent: it does nothing if
the output already exists unless `--force` is passed.

Usage:
    python3 scripts/00_build_gse103224_raw.py [--force]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib


def main():
    out = os.path.join(lib.PROC_DIR, "gse103224_raw.h5ad")
    if os.path.exists(out) and "--force" not in sys.argv:
        print(f"[step00] {out} already exists — skipping (pass --force to rebuild)")
        return

    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step00_build_raw.log"))
    os.makedirs(lib.PROC_DIR, exist_ok=True)
    lib.log.info(f"GSE103224 matrices: {lib.GSE103224_DIR}")
    counts, idmap = lib.load_gse103224()
    adata = lib.to_anndata(counts, idmap)
    adata.write_h5ad(out)
    lib.log.info(f"wrote {out}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
