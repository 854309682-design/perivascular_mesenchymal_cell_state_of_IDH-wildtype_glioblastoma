#!/usr/bin/env python3
"""rev16 — Supplementary Figure S2 rebuilt from the revised survival analysis (review M6).

Panels
  (a) cohort hazard ratios per SD of the MES score, unadjusted and grade-adjusted, plus the
      pooled random-effects (primary) and fixed-effects estimates
  (b) Kaplan-Meier, TCGA (median cut-off), with log-rank P and numbers at risk
  (c) Kaplan-Meier, CGGA-325 + CGGA-693 pooled
  (d) MES score by WHO grade (TCGA) with the Kruskal-Wallis statistic
  (e) IDH-wildtype vs IDH-mutant MES-like proportion per cohort, with the pooled and the
      cohort-mean-centred tests
  (f) overlap of the two inference routes' top-30 MES-linked lists with its Fisher exact test
      (review M2/DA3: the agreement of the two routes is quantified, not asserted). Added
      2026-10-14 together with the de-novo arm; reads results/revision/rev07_route_overlap.tsv.

Kaplan-Meier, the log-rank test and the at-risk table are computed here directly (lifelines is
not installed in this environment), so the panel matches the numbers reported in the text.
Outputs: figures/FigureS2_supp.{pdf,jpeg} (JPEG + vector PDF, no grid lines).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import lib

RES = os.path.join(lib.PROJ, "results", "revision")
FIG = os.path.join(lib.PROJ, "figures")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": False, "font.family": "DejaVu Sans"})
C1, C2 = "#C0392B", "#2C3E50"


def km(d, tcol="time", ecol="event"):
    """Kaplan-Meier estimate. Returns (times, survival, at-risk dict)."""
    d = d.dropna(subset=[tcol, ecol]).sort_values(tcol)
    t = d[tcol].to_numpy(float)
    e = d[ecol].to_numpy(int)
    times, surv = [0.0], [1.0]
    s, n = 1.0, len(t)
    for i, ti in enumerate(t):
        if e[i] == 1:
            s *= (1 - 1 / (n - i))
            times.append(ti); surv.append(s)
    return np.array(times), np.array(surv)


def logrank(d, gcol, tcol="time", ecol="event"):
    """Two-group log-rank test."""
    g = d[gcol].values
    groups = np.unique(g)
    times = np.unique(d.loc[d[ecol] == 1, tcol].dropna())
    O = {k: 0.0 for k in groups}
    E = {k: 0.0 for k in groups}
    V = 0.0
    for ti in times:
        at = (d[tcol] >= ti).values
        n = at.sum()
        d_i = ((d[tcol] == ti) & (d[ecol] == 1)).values.sum()
        n1 = (at & (g == groups[0])).sum()
        d1 = ((d[tcol] == ti) & (d[ecol] == 1) & (g == groups[0])).values.sum()
        O[groups[0]] += d1
        e1 = d_i * n1 / n
        E[groups[0]] += e1
        if n > 1:
            V += (n1 / n) * (1 - n1 / n) * d_i * (n - d_i) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi2 = (O[groups[0]] - E[groups[0]]) ** 2 / V
    return chi2, 1 - stats.chi2.cdf(chi2, 1)


def atrisk(d, tcol, ecol):
    d = d.dropna(subset=[tcol])
    tmax = np.nanmax(d[tcol])
    pts = np.linspace(0, tmax, 5)
    return [(int(np.sum(d[tcol] >= p))) for p in pts], [f"{int(p)}" for p in pts]


def main():
    lib.setup_logging()
    t = pd.read_csv(os.path.join(RES, "rev05_tcga_cases.tsv"), sep="\t")
    c3 = pd.read_csv(os.path.join(RES, "rev05_cgga325_cases.tsv"), sep="\t")
    c6 = pd.read_csv(os.path.join(RES, "rev05_cgga693_cases.tsv"), sep="\t")
    cox = pd.read_csv(os.path.join(RES, "rev06_cox_models.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(RES, "rev06_meta_analysis.tsv"), sep="\t")

    t = t.rename(columns={"os_time_days": "time"})
    t = t[np.isfinite(t["time"]) & (t["time"] > 0)]
    c3 = c3.rename(columns={"os_days": "time"})
    c6 = c6.rename(columns={"os_days": "time"})
    for c in (c3, c6):
        c["time"] = pd.to_numeric(c["time"], errors="coerce")
        c["event"] = pd.to_numeric(c["event"], errors="coerce")
    cgga = pd.concat([c3.assign(set="CGGA-325"), c6.assign(set="CGGA-693")])
    cgga = cgga[np.isfinite(cgga["time"]) & (cgga["time"] > 0)]

    fig = plt.figure(figsize=(12.5, 10.2))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.15, 1.0, 1.0, 0.85],
                          hspace=0.62, wspace=0.28)

    # (a) forest
    ax = fig.add_subplot(gs[0, :])
    rows = cox[cox["model"].str.contains("age \\+ IDH|WHO grade", regex=True)].copy()
    pooled = [("pooled, random effects (primary)", meta.loc[0, "HR"], meta.loc[0, "lo"],
               meta.loc[0, "hi"], "RE"),
              ("pooled, fixed effects", meta.loc[1, "HR"], meta.loc[1, "lo"], meta.loc[1, "hi"], "FE")]
    labels = list(rows["model"]) + [p[0] for p in pooled]
    vals = list(zip(rows["HR_per_SD"], rows["lo"], rows["hi"], ["cohort"] * len(rows))) + \
           [(p[1], p[2], p[3], p[4]) for p in pooled]
    y = np.arange(len(vals))[::-1]
    for yi, (hr, lo, hi, kind) in zip(y, vals):
        c = {"cohort": C1, "RE": "#8E44AD", "FE": "#16A085"}[kind]
        ax.errorbar(hr, yi, xerr=[[hr - lo], [hi - hr]], fmt="o", ms=4, lw=1.2, capsize=2.5, color=c)
        ax.text(hi * 1.03, yi, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", va="center", fontsize=7.5)
    ax.axvline(1, color="grey", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("hazard ratio per SD of the MES score (95% CI)")
    ax.set_title("(a) Prognostic value across cohorts and models (grade-adjusted shown)", fontsize=9, loc="left")

    # (b) KM TCGA
    for j, (dat, name, panel) in enumerate([(t, "TCGA", "(b) TCGA"),
                                            (cgga, "CGGA-325 + CGGA-693", "(c) pooled CGGA")]):
        ax = fig.add_subplot(gs[1, j])
        dat = dat.copy()
        med = float(np.nanmedian(dat["MES_Neftel"]))
        dat["grp"] = np.where(dat["MES_Neftel"] > med, "MES high", "MES low")
        for grp, col in [("MES high", C1), ("MES low", C2)]:
            sub = dat[dat["grp"] == grp]
            xs, ys = km(sub)
            ax.step(xs, ys, where="post", color=col, lw=1.4,
                    label=f"{grp} (n = {len(sub)})")
        chi2, p = logrank(dat, "grp")
        ar, labs = atrisk(dat, "time", "event")
        ax.set_xlabel("time (days)")
        ax.set_ylabel("overall survival")
        ax.set_ylim(0, 1.02)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax.set_title(f"{panel}: {name}\nmedian split; log-rank chi2 = {chi2:.1f}, "
                     f"P = {p:.2g} (partly grade-confounded; see (a))", fontsize=8.5, loc="left")
        ax.text(0.02, -0.30, "n at risk: " + "  ".join(f"{l}:{a}" for l, a in zip(labs, ar)),
                transform=ax.transAxes, fontsize=7, color="grey")

    # (d) grade trend
    ax = fig.add_subplot(gs[2, 0])
    g = t[np.isfinite(t["grade_num"]) & np.isfinite(t["MES_Neftel"])]
    groups = [gg["MES_Neftel"].values for _, gg in g.groupby("grade_num")]
    kt = stats.kruskal(*groups)
    ax.boxplot(groups, tick_labels=[f"G{int(x)}" for x in sorted(g['grade_num'].unique())],
               showfliers=False, widths=0.5)
    ax.scatter(np.concatenate([np.full(len(v), i + 1) + np.random.default_rng(i).normal(0, 0.06, len(v))
                               for i, v in enumerate(groups)]),
               np.concatenate(groups), s=5, color=C1, alpha=0.5, linewidths=0)
    ax.set_xlabel("WHO grade"); ax.set_ylabel("MES score")
    ax.set_title(f"(d) MES score rises with grade (TCGA)\nKruskal-Wallis P = {kt.pvalue:.2g}",
                 fontsize=9, loc="left")

    # (e) IDH per-cohort means
    ax = fig.add_subplot(gs[2, 1])
    idh = pd.read_csv(os.path.join(RES, "rev15_idh_stratified.tsv"), sep="\t")
    per = idh[idh["test"] == "per_cohort_means"]
    cols = [C2 if "wt" in gg else C1 for gg in per["genotype"]]
    ax.bar(range(len(per)), per["mean_MES"], color=cols, alpha=0.9)
    for i, (_, r) in enumerate(per.iterrows()):
        ax.text(i, r["mean_MES"] + 0.012, f"n={int(r['n_patients'])}", ha="center", fontsize=7)
    ax.set_xticks(range(len(per)))
    ax.set_xticklabels([c.replace("_", "\n") for c in per["cohort"]], fontsize=7)
    ax.set_ylabel("mean MES-like proportion")
    pooled_p = float(idh.loc[idh["test"] == "pooled_Mann_Whitney", "p"].iloc[0])
    cent_p = float(idh.loc[idh["test"] == "cohort_mean_centred_Mann_Whitney", "p"].iloc[0])
    ax.set_title("(e) IDH-wildtype (dark) vs IDH-mutant (red) per cohort\n"
                 f"pooled P = {pooled_p:.1g}; cohort-mean-centred P = {cent_p:.2g} "
                 "(genotype is confounded with cohort)", fontsize=9, loc="left")

    # (f) overlap of the two routes' top-30 lists (review M2/DA3)
    ax = fig.add_subplot(gs[3, :])
    ov = pd.read_csv(os.path.join(RES, "rev07_route_overlap.tsv"), sep="\t").iloc[0]
    cats = ["prior route only", "shared", "de-novo route only"]
    vals = [int(ov["prior_only"]), int(ov["shared"]), int(ov["denovo_only"])]
    cols = [C2, "#8E44AD", C1]
    bars = ax.barh(range(3), vals, color=cols, alpha=0.9, height=0.62)
    for i, (b, v) in enumerate(zip(bars, vals)):
        ax.text(v + 0.4, i, str(v), va="center", fontsize=9)
    ax.set_yticks(range(3)); ax.set_yticklabels(cats, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("number of top-30 MES-linked transcription factors")
    ax.set_xlim(0, max(vals) * 1.15)
    ax.set_title(f"(f) Overlap of the two routes' top-30 MES-linked lists: "
                 f"{int(ov['shared'])} shared of {int(ov['universe'])} scored TFs "
                 f"(Fisher exact P = {float(ov['fisher_p_greater']):.2f}, odds ratio "
                 f"{float(ov['odds_ratio']):.2f})\n"
                 f"shared factors: {ov['shared_factors']}", fontsize=9, loc="left")

    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"FigureS2_supp.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote FigureS2_supp.pdf/.jpeg")
    print(f"TCGA log-rank: {logrank(t, 'grp')[1]:.3g}; pooled CGGA: {logrank(cgga, 'grp')[1]:.3g}"
          if False else "")
    print("DONE rev16")


if __name__ == "__main__":
    main()
