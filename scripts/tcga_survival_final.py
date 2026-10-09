#!/usr/bin/env python3
"""Valid censored TCGA survival: MES score vs overall survival, IDH-stratified.
Merges the MES score table with per-case XML survival (proper censoring), then
Kaplan-Meier (IDH-wt / IDH-mut), log-rank, and multivariable Cox (MES + age + IDH + grade).
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

CLIN = os.path.join(lib.PROJ, "results", "step07_clinical")
OUT = CLIN


def load():
    mes = pd.read_csv(os.path.join(CLIN, "tcga_mes_clinical.tsv"), sep="\t")
    xml = pd.read_csv(os.path.join(lib.RAW_DIR, "tcga", "clinical_xml", "clinical_xml_survival.tsv"), sep="\t")
    xml = xml.rename(columns={"submitter_id": "sample"})
    d = mes.merge(xml, on="sample", how="inner", suffixes=("", "_xml"))
    d["days_to_death"] = pd.to_numeric(d["days_to_death"], errors="coerce")
    d["days_to_last_followup"] = pd.to_numeric(d["days_to_last_followup"], errors="coerce")
    d["vital_status"] = d["vital_status"].astype(str)
    # OS: event = death, time = death or last follow-up
    d["os_event"] = (d["vital_status"].str.lower() == "dead").astype(int)
    d["os_time"] = np.where(d["os_event"] == 1, d["days_to_death"], d["days_to_last_followup"])
    d = d[(d["os_time"] > 0) & d["os_time"].notna()].copy()
    return d


def main():
    d = load()
    print(f"survival cases: {len(d)} | events(death): {int(d['os_event'].sum())} | "
          f"censored: {int((d['os_event']==0).sum())}", flush=True)
    d["MES_group"] = np.where(d["MES_score"] >= d["MES_score"].median(), "high", "low")

    def lr(sub, name):
        if sub["MES_group"].nunique() < 2:
            return
        g = sub.groupby("MES_group")
        r = logrank_test(g.get_group("high")["os_time"], g.get_group("low")["os_time"],
                         g.get_group("high")["os_event"], g.get_group("low")["os_event"])
        print(f"[{name}] logrank p={r.p_value:.4f} "
              f"(high n={len(g.get_group('high'))}, low n={len(g.get_group('low'))})", flush=True)

    lr(d, "all")
    lr(d[d["IDH_mutated"] == 0], "IDH-wt")
    lr(d[d["IDH_mutated"] == 1], "IDH-mut")
    lr(d[d["project"] == "TCGA-GBM"], "TCGA-GBM")
    lr(d[d["project"] == "TCGA-LGG"], "TCGA-LGG")

    # multivariable Cox
    cdf = d.dropna(subset=["os_time", "os_event", "MES_score", "age_at_index"]).copy()
    cdf["IDH"] = cdf["IDH_mutated"].astype(float)
    cols = ["MES_score", "age_at_index", "IDH"]
    cph = CoxPHFitter()
    cph.fit(cdf[["os_time", "os_event"] + cols], "os_time", "os_event")
    cph.summary.to_csv(os.path.join(OUT, "tcga_cox.tsv"), sep="\t")
    print("=== Cox (MES + age + IDH) ===\n" +
          cph.summary[["coef", "exp(coef)", "p", "exp(coef) lower 95%", "exp(coef) upper 95%"]].round(3).to_string(), flush=True)

    # save merged survival table (before plotting)
    d[["sample", "project", "MES_score", "MES_group", "IDH_mutated", "age_at_index",
       "os_time", "os_event", "vital_status"]].to_csv(os.path.join(OUT, "tcga_survival.tsv"),
                                                      sep="\t", index=False)

    # KM figure (manual step curves to avoid lifelines/matplotlib hang)
    try:
        def km_curve(times, events):
            order = np.argsort(times)
            t = np.asarray(times)[order]; e = np.asarray(events)[order]
            n = len(t); surv = 1.0; xs = [0]; ys = [1.0]
            for i in range(n):
                if e[i] == 1:
                    surv *= (1 - 1.0 / (n - i))
                    xs.append(t[i]); ys.append(surv)
            return xs, ys
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        for col, idhv, label in [(0, 1, "IDH-mutant"), (1, 0, "IDH-wildtype")]:
            sub = d[d["IDH_mutated"] == idhv]
            ax[col].set_title(f"{label} (n={len(sub)})")
            for grp, c in [("high", "#C0392B"), ("low", "#2980B9")]:
                g = sub[sub["MES_group"] == grp]
                if len(g) < 3:
                    continue
                xs, ys = km_curve(g["os_time"].values, g["os_event"].values)
                ax[col].step(xs, ys, where="post", c=c, label=f"MES {grp} (n={len(g)})")
            ax[col].set_xlabel("days"); ax[col].set_ylabel("overall survival"); ax[col].legend()
        plt.tight_layout(); plt.savefig(os.path.join(OUT, "tcga_km.png"), dpi=150)
        plt.close("all")
        print("KM figure saved", flush=True)
    except Exception as e:
        print("KM fig error:", e, flush=True)

    print("DONE", flush=True)


if __name__ == "__main__":
    main()
