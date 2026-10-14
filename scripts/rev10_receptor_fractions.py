#!/usr/bin/env python3
"""rev10 — Receptor expression on malignant cells (review M3 minor 6 / item 2.8).

The submitted legend stated "CSF1R 7%" but CSF1R was not plotted, and the text called
CSF1 the top ligand while placing its "downstream program in tumour cells" — incoherent
unless a myeloid intermediary is posited. This script reports the fraction of malignant
GSE182109 cells with non-zero counts for every receptor named in the text, with n.

Output: results/revision/rev10_receptor_fractions.tsv
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
import lib
import mc_cohort_state as mcs
import mc_gse182109_h1 as g182

sc.settings.verbosity = 0
RECEPTORS = ["IL6ST", "IL6R", "LRP1", "ITGB1", "EGFR", "ERBB2", "BMPR2", "BMPR1A",
             "CD44", "FGFR1", "TGFBR1", "TGFBR2", "CSF1R", "CSF2RA", "ACVR2A", "ACVR2B",
             "NOTCH1", "PTPRC", "SDC1", "NRP1", "KDR", "PECAM1"]


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev10_receptor.log"))
    files = sorted(glob.glob(os.path.join(lib.COHORTS, "gse182109", "*_matrix.mtx.gz")))
    ads = []
    for f in files:
        x = g182.load_one(f)
        smp = os.path.basename(f).split("_", 1)[1].replace("_matrix.mtx.gz", "")
        x.obs["type"] = smp.split("-")[0]
        x.obs["patient"] = "-".join(smp.split("-")[:2])
        ads.append(x)
    a = sc.concat(ads, join="outer")
    a.var_names_make_unique()
    sc.pp.filter_cells(a, min_genes=200)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)              # marker scoring must run on normalised data (Mann-Whitney
                                # style raw-count scoring mis-ranks the compartments)
    # malignant identification exactly as in the pipeline (marker-score exclusion)
    for cname, genes in [("immune", mcs.IMMUNE), ("endothelial", mcs.ENDO)]:
        present = [g for g in genes if g in a.var_names]
        sc.tl.score_genes(a, present, score_name=f"score_{cname}", ctrl_size=50)
    a.obs["nonmalignant"] = ((a.obs["score_immune"] > 0.5) | (a.obs["score_endothelial"] > 0.5)).astype(int)
    mal = a[a.obs["nonmalignant"] == 0]
    n = mal.n_obs
    rows = []
    for r in RECEPTORS:
        if r not in mal.var_names:
            rows.append({"receptor": r, "pct_malignant_expressing": np.nan, "n_malignant": n,
                         "note": "not in matrix"})
            continue
        v = mal[:, r].layers["counts"].toarray().ravel() if "counts" in mal.layers else \
            np.asarray(mal[:, r].X.todense()).ravel()
        rows.append({"receptor": r, "pct_malignant_expressing": round(float((v > 0).mean() * 100), 1),
                     "n_malignant": n, "note": ""})
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(lib.PROJ, "results", "revision", "rev10_receptor_fractions.tsv"),
             sep="\t", index=False)
    lib.log.info(f"malignant GSE182109 cells: {n}\n" + d.to_string(index=False))
    lib.log.info("DONE rev10")


if __name__ == "__main__":
    main()
