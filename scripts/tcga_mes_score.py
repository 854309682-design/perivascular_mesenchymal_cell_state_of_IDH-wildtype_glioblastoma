#!/usr/bin/env python3
"""Build a TCGA MES-like signature score (Neftel MES1+MES2 meta-module, z-scored
mean of log1p-TPM) and merge with clinical (survival/age/grade) + IDH status.
Output: results/step07_clinical/tcga_mes_clinical.tsv (one row per case).
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules

sc.settings.verbosity = 0
TCGA = os.path.join(lib.RAW_DIR, "tcga")
OUT = os.path.join(lib.PROJ, "results", "step07_clinical")
os.makedirs(OUT, exist_ok=True)


def main():
    modules = load_neftel_modules()
    mes_genes = sorted(set(modules.get("MES1", []) + modules.get("MES2", [])))
    print("MES meta-module genes:", len(mes_genes), flush=True)

    man = pd.read_csv(os.path.join(TCGA, "manifest.tsv"), sep="\t")
    # map file_name -> sample submitter_id + project
    fn2sample = dict(zip(man["file_name"], man["sample"]))
    fn2proj = dict(zip(man["file_name"], man["project"]))

    expr_files = sorted(glob.glob(os.path.join(TCGA, "*rna_seq.augmented_star_gene_counts.tsv")))
    print("expression files:", len(expr_files), flush=True)

    # gather log1p TPM for MES genes per sample
    rows = []
    for f in expr_files:
        fname = os.path.basename(f)
        sample = fn2sample.get(fname)
        proj = fn2proj.get(fname)
        if sample is None:
            continue
        df = pd.read_csv(f, sep="\t", comment="#", usecols=["gene_name", "tpm_unstranded"])
        df = df[df["gene_name"].isin(mes_genes)]
        vals = {g: v for g, v in zip(df["gene_name"], df["tpm_unstranded"])}
        row = {"sample": sample, "project": proj}
        # log1p TPM; missing -> NaN
        for g in mes_genes:
            row[f"gene_{g}"] = np.log1p(vals.get(g, 0.0))
        rows.append(row)
    xs = pd.DataFrame(rows)
    print("cases with MES expression:", xs.shape, flush=True)

    # z-score each gene across samples, then mean -> MES score
    gene_cols = [c for c in xs.columns if c.startswith("gene_")]
    z = (xs[gene_cols] - xs[gene_cols].mean()) / xs[gene_cols].std()
    xs["MES_score"] = z.mean(axis=1)
    xs["n_genes_present"] = xs[gene_cols].notna().sum(axis=1)
    xs = xs[["sample", "project", "MES_score", "n_genes_present"]]
    print("MES score range:", xs["MES_score"].min().round(3), "-", xs["MES_score"].max().round(3),
          "| median", xs["MES_score"].median().round(3), flush=True)

    # merge clinical
    clin = pd.read_csv(os.path.join(TCGA, "clinical", "tcga_clinical.tsv"), sep="\t")
    clin["days_to_death"] = pd.to_numeric(clin["days_to_death"], errors="coerce")
    clin["days_to_last_follow_up"] = pd.to_numeric(clin["days_to_last_follow_up"], errors="coerce")
    clin["age_at_index"] = pd.to_numeric(clin["age_at_index"], errors="coerce")
    # overall survival time + event (use vital_status)
    def os_time(r):
        if str(r["vital_status"]).strip().lower() == "dead" and pd.notna(r["days_to_death"]):
            return r["days_to_death"], 1
        return r["days_to_last_follow_up"], 0
    clin[["os_time", "os_event"]] = clin.apply(lambda r: pd.Series(os_time(r)), axis=1)
    clin = clin[["submitter_id", "project", "age_at_index", "gender", "tumor_grade",
                 "histological_type", "os_time", "os_event"]]
    clin = clin.rename(columns={"submitter_id": "sample"})

    # IDH status (full aliquot barcode -> short case id)
    idh = pd.read_csv(os.path.join(TCGA, "mutation", "idh_status.tsv"), sep="\t")
    idh["gene"] = idh["gene"].astype(str)
    idh["case"] = idh["sample"].str.split("-").str[:3].str.join("-")
    idh_map = dict(zip(idh["case"], idh["gene"]))
    clin["IDH_mutated"] = clin["sample"].map(lambda s: 1 if idh_map.get(s) else 0)

    merged = xs.merge(clin.drop(columns=["project"]), on="sample", how="left").merge(
        idh[["case", "gene"]].rename(columns={"case": "sample", "gene": "IDH_genes"}),
        on="sample", how="left")
    merged.to_csv(os.path.join(OUT, "tcga_mes_clinical.tsv"), sep="\t", index=False)
    print("merged:", merged.shape, "| cases with os_time:", merged["os_time"].notna().sum(),
          "| cases IDH_mutated:", merged["IDH_mutated"].sum(), flush=True)
    print(merged.groupby("project")["MES_score"].agg(["count", "mean", "median"]).round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
