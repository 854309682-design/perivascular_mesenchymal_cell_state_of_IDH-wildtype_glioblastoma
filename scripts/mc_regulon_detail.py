#!/usr/bin/env python3
"""Regulon detail for GSE70630 (IDH-mut) and Neftel (IDH-wt): per-cohort master-TF
ranking + candidate-TF comparison, from the saved per-cohort DoRothEA TF-activity
results (h2_tf_*.tsv).
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib

OUT = os.path.join(lib.PROJ, "results", "multicohort")
CAND = ["STAT3", "FOSL1", "FOSL2", "JUNB", "JUN", "NFKB1", "CEBPB", "TEAD1", "TEAD4",
        "RELA", "IRF1", "HIF1A", "ATF4", "SP1", "FOXL2"]
COHORTS = [("Neftel", "h2_tf_Neftel.tsv"), ("Venteicher", "h2_tf_Venteicher.tsv")]


def load(path):
    df = pd.read_csv(path, sep="\t")
    return df


def main():
    frames = {}
    for name, f in COHORTS:
        p = os.path.join(OUT, f)
        if not os.path.exists(p):
            print(f"missing {p}", flush=True)
            continue
        df = load(p)
        frames[name] = df
    if not frames:
        return

    # per-cohort top master-TFs (by |rho|)
    detail_lines = ["# Regulon detail — GSE70630 (IDH-mut) vs Neftel (IDH-wt)"]
    for name, df in frames.items():
        top = df.reindex(df["rho"].abs().sort_values(ascending=False).index).head(10)
        detail_lines.append(f"\n## {name} — top 10 master-TF candidates (|rho| desc)")
        detail_lines.append(top[["tf", "rho", "fdr"]].to_string(index=False))
        top.to_csv(os.path.join(OUT, f"regulon_top_{name}.tsv"), sep="\t", index=False)

    # candidate comparison (mean rho per cohort)
    cand_rows = []
    for name, df in frames.items():
        sub = df[df["tf"].isin(CAND)].set_index("tf")["rho"]
        for tf in CAND:
            cand_rows.append({"tf": tf, "cohort": name, "rho": sub.get(tf, np.nan)})
    cand = pd.DataFrame(cand_rows).pivot(index="tf", columns="cohort", values="rho").round(3)
    cand.to_csv(os.path.join(OUT, "regulon_candidate_TFs.tsv"), sep="\t")
    detail_lines.append("\n## Candidate master-TF activity (rho vs MES)")
    detail_lines.append(cand.to_string())
    open(os.path.join(OUT, "regulon_detail.md"), "w").write("\n".join(detail_lines))
    print("\n".join(detail_lines), flush=True)

    # figure: grouped bar of candidate TFs rho per cohort
    try:
        cols = list(cand.columns)
        x = np.arange(len(cand.index))
        w = 0.35
        fig, ax = plt.subplots(figsize=(9, 5))
        for i, c in enumerate(cols):
            ax.bar(x + (i - 0.5) * w, cand[c].values, width=w, label=c,
                   color=["#C0392B", "#2980B9"][i % 2])
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(cand.index, rotation=45, ha="right")
        ax.set_ylabel("rho (TF activity vs MES score)"); ax.legend()
        ax.set_title("Candidate master-TF activity: IDH-wt (Neftel) vs IDH-mut (GSE70630)")
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "regulon_detail.png"), dpi=150)
        plt.close("all")
        print("figure saved", flush=True)
    except Exception as e:
        print("figure error", e, flush=True)


if __name__ == "__main__":
    main()
