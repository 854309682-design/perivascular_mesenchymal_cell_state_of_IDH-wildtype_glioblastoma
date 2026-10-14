#!/usr/bin/env python3
"""Figure 2 — the MES program is IDH-wildtype-specific and recurrence-independent.

Rebuilds the six-panel figure from the per-patient tables (JPEG + vector PDF, no
grid lines), replacing the panel (c) of the earlier, unpublished figure script
that plotted the IDH-mutant arm as 240 one-cell "patients" and printed the
inverted P (1.6e-57).  Panels:

  (a) mean MES-like proportion per cohort
  (b) per-patient MES-like proportion by cohort (boxplots, patients overlaid)
  (c) IDH-wildtype vs IDH-mutant mean MES-like proportion (corrected, n = 69 vs 16)
  (d) H1: recurrent vs primary GSE182109 (null)
  (e) GSE182109 per-patient MES-like proportion by tumour type
  (f) cohort x IDH-genotype heatmap of mean MES-like proportion

Inputs
  results/multicohort/idh_patients.tsv                    (written by mc_gse89567_idh.py)
  results/multicohort/gse182109_composition_by_patient.tsv
  results/multicohort/gse182109_h1_primary_vs_recurrent.tsv

Outputs
  figures/Figure2_main.pdf + figures/Figure2_main.jpeg
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

MC = os.path.join(lib.PROJ, "results", "multicohort")
FIGDIR = os.path.join(lib.PROJ, "figures")
RED, BLUE = "#C0392B", "#2E86C1"
BOX, MEDIAN = "#A9C7E8", "#E67E22"

# display names, keyed by the cohort column of idh_patients.tsv
DISPLAY = {
    "GSE103224": "GSE103224",
    "GSE182109": "GSE182109",
    "Neftel_GSE131928": "GSE131928",
    "GSE70630": "GSE70630",
    "GSE89567": "GSE89567",
}
# cohort -> IDH genotype, for panel (f)
GENOTYPE = {
    "GSE103224": "IDH-wt", "GSE182109": "IDH-wt", "Neftel_GSE131928": "IDH-wt",
    "GSE70630": "IDH-mut", "GSE89567": "IDH-mut",
}


def load():
    pat = pd.read_csv(os.path.join(MC, "idh_patients.tsv"), sep="\t")
    pat["cohort"] = pat["cohort"].map(DISPLAY)
    h1 = pd.read_csv(os.path.join(MC, "gse182109_h1_primary_vs_recurrent.tsv"),
                     sep="\t").iloc[0]
    g182 = pd.read_csv(os.path.join(MC, "gse182109_composition_by_patient.tsv"),
                       sep="\t")
    return pat, h1, g182


def panel_a(ax, pat, order):
    means = [pat.loc[pat["cohort"] == c, "MES_prop"].mean() for c in order]
    ax.bar(order, means, color=RED)
    ax.set_ylabel("mean MES prop")
    ax.set_title("(a) mean MES-like per cohort", fontsize=12)
    ax.tick_params(axis="x", rotation=30)
    for lab in ax.get_xticklabels():
        lab.set_ha("right")
    ax.set_ylim(0, max(means) * 1.08)


def panel_b(ax, pat, order):
    data = [pat.loc[pat["cohort"] == c, "MES_prop"].values for c in order]
    bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.6,
                    boxprops=dict(facecolor=BOX, edgecolor="black"),
                    medianprops=dict(color=MEDIAN, linewidth=2),
                    whiskerprops=dict(color="black"), capprops=dict(color="black"),
                    flierprops=dict(marker="none"))
    del bp
    for i, vals in enumerate(data, start=1):
        x = np.random.default_rng(42).normal(i, 0.06, len(vals))
        ax.plot(x, vals, "o", ms=3, color="#7F8C8D", alpha=0.6, zorder=3)
    ax.axhline(0, ls=":", lw=1, color="black")
    ax.set_ylabel("MES prop")
    ax.set_title("(b) per-patient MES-like by cohort", fontsize=12)
    ax.tick_params(axis="x", rotation=30)
    for lab in ax.get_xticklabels():
        lab.set_ha("right")


def panel_c(ax, pat, pval):
    """Group means WITH the per-patient points behind them (review minor 20)."""
    wt = pat.loc[pat["genotype"] == "IDH-wt", "MES_prop"].values
    mut = pat.loc[pat["genotype"] == "IDH-mut", "MES_prop"].values
    ax.bar([0, 1], [wt.mean(), mut.mean()], color=[RED, BLUE], width=0.55)
    for i, (vals, col) in enumerate([(wt, RED), (mut, BLUE)]):
        x = np.random.default_rng(7).normal(i, 0.055, len(vals))
        ax.plot(x, vals, "o", ms=3, color="#2C3E50", alpha=0.55, zorder=3)
        ax.text(i, vals.mean() + 0.01, f"n = {len(vals)}", ha="center", fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["IDH-wt", "IDH-mut"])
    ax.set_ylabel("MES prop")
    ax.set_title(f"(c) IDH-wt vs mut MES (p={pval:.1e})", fontsize=12)


def panel_d(ax, h1, g182):
    groups = [g182.loc[g182["type"] == t, "MES_prop"].values
              for t in ("ndGBM", "rGBM")]
    means = [g.mean() for g in groups]
    sds = [g.std(ddof=1) for g in groups]
    ax.bar(["primary", "recurrent"], means, yerr=sds, color=[BLUE, RED],
           capsize=0, error_kw=dict(ecolor="black", lw=2))
    ax.set_ylabel("mean MES prop")
    ax.set_title(f"(d) H1: recurrent vs primary (p={h1['p_alt_greater']:.2f})",
                 fontsize=12)


def panel_e(ax, g182):
    for t, colour in (("ndGBM", BLUE), ("rGBM", RED)):
        vals = np.sort(g182.loc[g182["type"] == t, "MES_prop"].values)
        x = np.random.default_rng(0).normal(
            0 if t == "ndGBM" else 1, 0.05, len(vals))
        ax.scatter(x, vals, color=colour, s=42, zorder=3)
    ax.set_xticks([0, 1], ["ndGBM", "rGBM"])
    ax.set_ylabel("MES prop")
    ax.set_title("(e) GSE182109 per-patient MES by type", fontsize=12)


def panel_f(ax, pat):
    cohorts = list(GENOTYPE)
    cols = ["IDH-mut", "IDH-wt"]
    mat = np.full((len(cohorts), 2), np.nan)
    for i, c in enumerate(cohorts):
        for j, g in enumerate(cols):
            if GENOTYPE[c] == g:
                mat[i, j] = pat.loc[pat["cohort"] == DISPLAY[c], "MES_prop"].mean()
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.1, vmax=0.6, aspect="auto")
    # mark the "no such genotype in this cohort" cells so they cannot be read as near-zero
    # values (review minor 21)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isnan(mat[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="0.85",
                                           edgecolor="0.4", hatch="///", lw=0.6))
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9, color="#444444")
            else:
                ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", fontsize=9,
                        color="black")
    ax.set_xticks(range(2), cols)
    ax.set_yticks(range(len(cohorts)), [DISPLAY[c] for c in cohorts])
    ax.set_title("(f) cohort x IDH mean MES", fontsize=12)
    ax.figure.colorbar(im, ax=ax)


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "fig2.log"))
    os.makedirs(FIGDIR, exist_ok=True)
    pat, h1, g182 = load()

    means = pat.groupby("cohort")["MES_prop"].mean().sort_values()
    order = list(means.index)
    lib.log.info(f"cohort order (ascending mean MES): "
                 f"{[(c, round(float(m), 4)) for c, m in means.items()]}")
    wt = pat.loc[pat["genotype"] == "IDH-wt", "MES_prop"]
    mut = pat.loc[pat["genotype"] == "IDH-mut", "MES_prop"]
    from scipy.stats import mannwhitneyu
    pval = mannwhitneyu(wt, mut, alternative="two-sided")[1]

    fig, axes = plt.subplots(2, 3, figsize=(13, 9.87))
    panel_a(axes[0, 0], pat, order)
    panel_b(axes[0, 1], pat, order)
    panel_c(axes[0, 2], pat, pval)
    panel_d(axes[1, 0], h1, g182)
    panel_e(axes[1, 1], g182)
    panel_f(axes[1, 2], pat)

    fig.suptitle("Figure 2 | The MES program is IDH-wildtype-specific and "
                 "recurrence-independent", fontsize=18, fontweight="bold")
    for ax in axes.ravel():
        ax.grid(False)                      # house style: no grid lines
    fig.tight_layout(rect=(0.01, 0.01, 0.99, 0.94))

    pdf = os.path.join(FIGDIR, "Figure2_main.pdf")
    jpg = os.path.join(FIGDIR, "Figure2_main.jpeg")
    fig.savefig(pdf)                        # vector
    fig.savefig(jpg, dpi=300)               # raster for on-screen review
    plt.close(fig)
    lib.log.info(f"wrote {pdf} and {jpg} "
                 f"(IDH-wt n={len(wt)} mean={wt.mean():.3f} vs "
                 f"IDH-mut n={len(mut)} mean={mut.mean():.3f}, P={pval:.2e})")


if __name__ == "__main__":
    main()
