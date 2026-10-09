#!/usr/bin/env python3
"""H3 clinical test: does the TCGA MES score predict survival, and is it
independent of IDH? IDH-stratified Kaplan-Meier (MES high/low) + multivariable Cox.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

sys.path.insert(0, os.path.dirname(__file__))
import lib

IN = os.path.join(lib.PROJ, "results", "step07_clinical", "tcga_mes_clinical.tsv")
OUT = os.path.join(lib.PROJ, "results", "step07_clinical")
os.makedirs(OUT, exist_ok=True)


def main():
    df = pd.read_csv(IN, sep="\t")
    df["os_time"] = pd.to_numeric(df["os_time"], errors="coerce")
    df["os_event"] = pd.to_numeric(df["os_event"], errors="coerce")
    df["age_at_index"] = pd.to_numeric(df["age_at_index"], errors="coerce")
    d = df[(df["os_time"] > 0) & df["os_time"].notna()].copy()
    print("cases with survival:", len(d), "| events:", int(d["os_event"].sum()), flush=True)

    # MES high/low by median (within whole cohort)
    d["MES_group"] = np.where(d["MES_score"] >= d["MES_score"].median(), "high", "low")

    # ---- log-rank (all, and IDH-stratified) ----
    def lr(sub, name):
        if sub["MES_group"].nunique() < 2:
            return
        g = sub.groupby("MES_group")
        res = logrank_test(g.get_group("high")["os_time"], g.get_group("low")["os_time"],
                           g.get_group("high")["os_event"], g.get_group("low")["os_event"])
        print(f"[{name}] logrank p={res.p_value:.4f} (n_high={len(g.get_group('high'))}, "
              f"n_low={len(g.get_group('low'))})", flush=True)
        return res

    lr(d, "all")
    for idhv, label in [(1, "IDH-mut"), (0, "IDH-wt")]:
        sub = d[d["IDH_mutated"] == idhv]
        lr(sub, label)

    # ---- KM plot, IDH-stratified ----
    try:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for col, idhv, label in [(0, 1, "IDH-mut"), (1, 0, "IDH-wt")]:
            sub = d[d["IDH_mutated"] == idhv]
            ax[col].set_title(label)
            for grp, c in [("high", "#C0392B"), ("low", "#2980B9")]:
                g = sub[sub["MES_group"] == grp]
                if len(g) < 3:
                    continue
                kmf = KaplanMeierFitter()
                kmf.fit(g["os_time"], g["os_event"], label=f"MES {grp} (n={len(g)})")
                kmf.plot_survival_function(ax=ax[col], c=c, ci_show=True)
            ax[col].set_xlabel("days"); ax[col].set_ylabel("survival probability")
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "tcga_km.png"), dpi=150)
        plt.close("all")
        print("KM figure saved", flush=True)
    except Exception as e:
        print("KM fig error:", e, flush=True)

    # ---- multivariable Cox ----
    cdf = d[(d["age_at_index"].notna())].copy()
    cdf["IDH"] = cdf["IDH_mutated"].astype(float)
    cols = [c for c in ["MES_score", "age_at_index", "IDH", "tumor_grade"] if c in cdf.columns]
    print("Cox n:", len(cdf), "| covariates:", cols, flush=True)
    cph = CoxPHFitter()
    cdf2 = cdf[["os_time", "os_event"] + cols].dropna()
    cph.fit(cdf2, "os_time", "os_event")
    summ = cph.summary
    summ.to_csv(os.path.join(OUT, "tcga_cox.tsv"), sep="\t")
    print("=== Cox PH (multivariable) ===", flush=True)
    print(summ[["coef", "exp(coef)", "p", "exp(coef) lower 95%", "exp(coef) upper 95%"]].round(3).to_string(), flush=True)
    # also cox with MES only
    cph2 = CoxPHFitter(); cph2.fit(d[["os_time", "os_event", "MES_score"]].dropna(), "os_time", "os_event")
    print("MES-only Cox HR:", round(cph2.summary.loc["MES_score", "exp(coef)"], 3),
          "p:", round(cph2.summary.loc["MES_score", "p"], 4), flush=True)

    # per-cohort (GBM vs LGG) log-rank
    for proj in ["TCGA-GBM", "TCGA-LGG"]:
        sub = d[d["project"] == proj]
        lr(sub, proj)


if __name__ == "__main__":
    main()
