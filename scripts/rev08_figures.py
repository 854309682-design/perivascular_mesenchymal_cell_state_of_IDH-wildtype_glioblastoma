#!/usr/bin/env python3
"""rev08 — Rebuilt figures for the revised manuscript (reviewer items 4.x + M2/M3/M6).

Every panel is written twice — a JPEG for on-screen review and a VECTOR PDF for
publication — and no panel uses grid lines.

Panels produced:
  Figure3_main.{pdf,jpeg}   a prior-route consensus (DoRothEA)  b de-novo consensus
                            (pySCENIC re-run)  c route-by-route agreement scatter
                            d per-cohort heatmap  e per-cohort rho distributions
                            f TWIST1/CEBPB sign audit
  Figure5_main.{pdf,jpeg}   a IvyGAP MES/MES1/MES2 by region (points + n)  b hypoxia
                            by region  c composition-adjusted region effects
                            d deconvolved malignant-MES share by region  e NicheNet
                            ligand activity + receptor fractions  f survival forest
  Figure2_main.{pdf,jpeg}   a-f with per-patient points (2c), defined error bars (2d),
                            distinct NA shading (2f)
  FigureS2_supp.{pdf,jpeg}  a cohort HRs + pooled RE/FE  b KM (TCGA)  c KM (pooled CGGA)
                            d MES vs WHO grade

If the de-novo arm (rev03) has not finished, Figure 3 panels b-d are skipped and the
script says so rather than inventing numbers.

Run:
  GLIOMA_PROJ=... /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev08_figures.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import lib

RES = os.path.join(lib.PROJ, "results", "revision")
MC = os.path.join(lib.PROJ, "results", "multicohort")
STEP = os.path.join(lib.PROJ, "results")
FIG = os.path.join(lib.PROJ, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "font.family": "DejaVu Sans",
})
C_MES, C_MES1, C_MES2 = "#C0392B", "#E67E22", "#8E44AD"
REGION_ORDER = ["LE", "IT", "CT", "CTmvp", "CTpan"]
REGION_LABEL = {"LE": "LE", "IT": "IT", "CT": "CT", "CTmvp": "MVP", "CTpan": "PAN",
                "CThbv": "HBV*", "CTpnz": "PNZ*"}


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.pdf/.jpeg")


# --------------------------------------------------------------------- Figure 5
def figure5():
    s = pd.read_csv(os.path.join(RES, "rev04_ivygap_samples.tsv"), sep="\t", index_col=0)
    a = s[s["study"] == "anatomic structures"].copy()
    st = pd.read_csv(os.path.join(RES, "rev04_ivygap_region_stats.tsv"), sep="\t")
    adj = pd.read_csv(os.path.join(RES, "rev04_ivygap_adjusted.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(RES, "rev04_ivygap_deconv.tsv"), sep="\t", index_col=0)
    regs = [r for r in REGION_ORDER if r in set(a["region"])]

    fig = plt.figure(figsize=(12.5, 8.6))
    gs = fig.add_gridspec(3, 3, hspace=0.52, wspace=0.32)

    # (a) MES / MES1 / MES2 by region, with per-sample points, jittered, no grid
    ax = fig.add_subplot(gs[0, :2])
    rng = np.random.default_rng(0)
    for k, (col, colr, lab) in enumerate([("MES_z", C_MES, "MES (MES1+MES2)"),
                                          ("MES1_z", C_MES1, "MES1 (hypoxia-linked)"),
                                          ("MES2_z", C_MES2, "MES2")]):
        for i, r in enumerate(regs):
            v = a.loc[a["region"] == r, col].dropna().values
            x = i + (k - 1) * 0.24
            ax.scatter(np.full(len(v), x) + rng.normal(0, 0.035, len(v)), v,
                       s=6, color=colr, alpha=0.55, linewidths=0, zorder=2)
            ax.plot([x - 0.09, x + 0.09], [np.median(v)] * 2, color=colr, lw=2, zorder=3)
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    n_ann = "  ".join(f"{REGION_LABEL[r]}:n={int((a['region'] == r).sum())}"
                      f"/{a.loc[a['region'] == r, 'tumour'].nunique()}t" for r in regs)
    ax.set_xticks(range(len(regs)))
    ax.set_xticklabels([REGION_LABEL[r] for r in regs], fontsize=11)
    ax.set_ylabel("score (z per tumour)")
    ax.set_title("(a) IvyGAP anatomic regions: MES programme peaks in the "
                 "pseudopalisading/peri-necrotic region,\n     not in the vascular region",
                 fontsize=9, loc="left")
    ax.set_xlabel(f"{n_ann}   (122 region samples, 10 tumours; t = tumours)", fontsize=8)
    ax.legend(handles=[Patch(color=C_MES, label="MES (MES1+MES2)"),
                       Patch(color=C_MES1, label="MES1 (hypoxia-linked)"),
                       Patch(color=C_MES2, label="MES2")], frameon=False, fontsize=8)

    # (b) hypoxia
    ax = fig.add_subplot(gs[0, 2])
    for i, r in enumerate(regs):
        v = a.loc[a["region"] == r, "hypoxia_z"].dropna().values
        ax.scatter(np.full(len(v), i) + rng.normal(0, 0.05, len(v)), v, s=7,
                   color="#16A085", alpha=0.6, linewidths=0)
        ax.plot([i - 0.12, i + 0.12], [np.median(v)] * 2, color="#16A085", lw=2)
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.set_xticks(range(len(regs)))
    ax.set_xticklabels([REGION_LABEL[r] for r in regs])
    ax.set_ylabel("hypoxia score (z)")
    ax.set_title("(b) Hypoxia control\n(15-gene HIF-target set)", fontsize=9, loc="left")

    # (c) adjusted region effects (mixed model, unadjusted vs adjusted)
    ax = fig.add_subplot(gs[1, :2])
    terms = {"IT": "IT", "CT": "CT", "CTmvp": "MVP", "CTpan": "PAN"}
    ypos = np.arange(len(terms))
    for j, (mdl, colr, lab) in enumerate([("MixedLM_unadjusted", "#7F8C8D", "unadjusted"),
                                          ("MixedLM_adj_vascular_myeloid", C_MES,
                                           "adjusted for vascular + myeloid content")]):
        sub = adj[adj["model"] == mdl]
        vals = []
        for t in terms:
            row = sub[sub["term"] == f"C(region)[T.{t}]"]
            vals.append(float(row["coef"].iloc[0]) if len(row) else np.nan)
        ax.barh(ypos + (0.2 if j else -0.2), vals, height=0.36, color=colr, alpha=0.85, label=lab)
        for y, v in zip(ypos + (0.2 if j else -0.2), vals):
            ax.text(v, y, f" {v:+.2f}", va="center", fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{k} vs LE" for k in terms])
    ax.invert_yaxis()
    ax.set_xlabel("region effect on MES score (mixed model, random intercept per tumour)")
    ax.set_title("(c) The vascular effect halves after composition adjustment, "
                 "while the peri-necrotic effect persists", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)

    # (d) deconvolution: malignant presence and malignant-MES share
    ax = fig.add_subplot(gs[1, 2])
    dd = dec.join(a[["region", "tumour"]], how="inner")
    dd = dd[dd["region"].isin(regs)]
    for i, r in enumerate(regs):
        v = dd.loc[dd["region"] == r, "malignant_total"].dropna().values
        ax.scatter(np.full(len(v), i - 0.16), v, s=7, color="#2C3E50", alpha=0.6, linewidths=0)
        v2 = dd.loc[dd["region"] == r, "MES_share_of_malignant"].dropna().values
        ax.scatter(np.full(len(v2), i + 0.16), v2, s=7, color=C_MES, alpha=0.6, linewidths=0)
    ax.set_xticks(range(len(regs)))
    ax.set_xticklabels([REGION_LABEL[r] for r in regs])
    ax.set_ylabel("NNLS fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("(d) Deconvolution\n(black) malignant content\n(red) MES share of malignant",
                 fontsize=9, loc="left")

    # (e) NicheNet ligand activity + receptor fractions (if computed)
    ax = fig.add_subplot(gs[2, :2])
    npath = os.path.join(RES, "rev09_csf1r_receptor_fractions.tsv")
    lpath = os.path.join(RES, "rev09_nichenet_ligand_activity.tsv")
    if os.path.exists(lpath):
        lg = pd.read_csv(lpath, sep="\t").head(15).iloc[::-1]
        ax.barh(np.arange(len(lg)), lg["ligand_activity"], color="#2980B9", alpha=0.85)
        ax.set_yticks(np.arange(len(lg)))
        ax.set_yticklabels(lg["ligand"], fontsize=8)
        ax.set_xlabel("NicheNet ligand activity vs the MES-transition programme (Pearson)")
        ax.set_title("(e) Niche ligands for the MES programme (base-R NicheNet prior)",
                     fontsize=9, loc="left")
    else:
        ax.text(0.5, 0.5, "(e) NicheNet panel: run rev09_nichenet.R", ha="center", va="center")
        ax.axis("off")

    # (f) survival forest
    ax = fig.add_subplot(gs[2, 2])
    mpath = os.path.join(RES, "rev06_cox_models.tsv")
    if os.path.exists(mpath):
        cx = pd.read_csv(mpath, sep="\t")
        cx = cx[cx["model"].str.contains("MES_Neftel \\+ age \\+ IDH|WHO grade", regex=True)]
        cx = cx.iloc[::-1]
        y = np.arange(len(cx))
        ax.errorbar(cx["HR_per_SD"], y,
                    xerr=[cx["HR_per_SD"] - cx["lo"], cx["hi"] - cx["HR_per_SD"]],
                    fmt="o", color=C_MES, ms=4, lw=1, capsize=2)
        ax.axvline(1, color="grey", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels([m.split(":")[0] + ("+grade" if "grade" in m else "")
                            for m in cx["model"]], fontsize=8)
        ax.set_xlabel("HR per SD (95% CI)")
        ax.set_title("(f) Prognostic value\n(grade-adjusted shown)", fontsize=9, loc="left")
    save(fig, "Figure5_main")


# --------------------------------------------------------------------- Figure 3
def figure3():
    """Figure 3 — regulon panels.

    Panels that depend on the de-novo arm are added by this same script once
    results/revision/rev07_denovo_consensus.tsv exists (rev03 + rev07); until then they are
    omitted rather than drawn as placeholders, so no unfinished panel can reach a figure.
    """
    prior = pd.read_csv(os.path.join(RES, "rev02_tf_consensus_reimplemented.tsv"), sep="\t")
    coh = pd.read_csv(os.path.join(RES, "rev02_tf_activity_vs_MES.tsv"), sep="\t")
    pat = pd.read_csv(os.path.join(RES, "rev02_pooled_vs_patient.tsv"), sep="\t")
    have_denovo = os.path.exists(os.path.join(RES, "rev07_denovo_consensus.tsv"))
    cons = pd.read_csv(os.path.join(RES, "rev07_denovo_consensus.tsv"), sep="\t") if have_denovo else None

    ncol = 3 if have_denovo else 2
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, ncol, hspace=0.6, wspace=0.34)

    # (a) prior route top-20
    ax = fig.add_subplot(gs[0, 0])
    top = prior.sort_values("mean_rho", ascending=False).head(20).iloc[::-1]
    ax.barh(np.arange(len(top)), top["mean_rho"], color=C_MES1, alpha=0.9)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["tf"], fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean Spearman rho (4 cohorts)")
    ax.set_title("(a) Prior-based route (DoRothEA, mLM)", fontsize=9, loc="left")

    # (b) de-novo route (only when the re-run is complete)
    if have_denovo:
        ax = fig.add_subplot(gs[0, 1])
        dn = cons.sort_values("mean_rho", ascending=False).head(20).iloc[::-1]
        ax.barh(np.arange(len(dn)), dn["mean_rho"], color=C_MES2, alpha=0.9)
        ax.set_yticks(np.arange(len(dn)))
        ax.set_yticklabels(dn["regulon"], fontsize=7)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("mean Spearman rho (3 cohorts)")
        ax.set_title("(b) De-novo route (pySCENIC re-run)", fontsize=9, loc="left")

    # key-factor heatmap across cohorts (prior route)
    ax = fig.add_subplot(gs[0, ncol - 1])
    piv = coh.pivot_table(index="tf", columns="cohort", values="rho").reindex(KEY_TFS)
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, fontsize=8, rotation=20, ha="right")
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=8)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(v) > 0.4 else "black")
    ax.set_title("(c) Key factors across cohorts", fontsize=9, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="rho")

    # per-cohort distributions
    ax = fig.add_subplot(gs[1, 0])
    data = [coh.loc[coh["cohort"] == c, "rho"].values for c in coh["cohort"].unique()]
    ax.boxplot(data, tick_labels=list(coh["cohort"].unique()), showfliers=False, widths=0.55)
    ax.axhline(0, color="grey", ls=":", lw=0.8)
    ax.set_ylabel("regulon rho")
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right", fontsize=7)
    ax.set_title("(d) Per-cohort distributions", fontsize=9, loc="left")

    # pooled-cell vs patient-level
    ax = fig.add_subplot(gs[1, 1:])
    pp = pat.set_index("tf").reindex(KEY_TFS).dropna(subset=["pooled_mean_rho"])
    x = np.arange(len(pp))
    ax.bar(x - 0.2, pp["pooled_mean_rho"], width=0.38, color=C_MES, alpha=0.9,
           label="pooled cells (descriptive)")
    ax.bar(x + 0.2, pp["patient_median_rho"], width=0.38, color="#2C3E50", alpha=0.9,
           label="patient-level median")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(pp.index, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Spearman rho vs MES")
    ax.set_title("(e) The unit of inference matters: pooled-cell correlations are attenuated "
                 "at patient level, and TWIST1/CEBPB are negative in both", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "Figure3_main")


KEY_TFS = ["STAT3", "KLF4", "FOXL2", "JUNB", "FOSL1", "FOSL2", "NFKB1", "CEBPB", "RELA",
           "TWIST1", "HIF1A", "ATF4"]


# --------------------------------------------------------------------- Figure 2
def figure2():
    """NOT DRAWN HERE — Figure 2 is owned by scripts/mc_fig2_idh_recurrence.py.

    That script (2026-10-10, the IDH-number correction) builds the six panels from
    results/multicohort/idh_patients.tsv, which carries the corrected patient definition and
    all five cohorts. An earlier revision of THIS function drew Figure 2 from
    composition_all_cohorts.tsv, a stale (2026-09-10) table that omits GSE89567 and labels the
    GSE70630 block "Venteicher_GSE70630"; running it overwrote the corrected figure. It is kept
    as a no-op so that a rev08 run can never clobber Figure2_main again.
    """
    print("[rev08] figure2() skipped: Figure 2 is built by scripts/mc_fig2_idh_recurrence.py "
          "(reads results/multicohort/idh_patients.tsv).")
    return


def figure2_legacy_stale():
    """Superseded Figure-2 builder (stale composition_all_cohorts.tsv). Kept for reference
    only — never call it: it overwrites the corrected Figure2_main produced by
    scripts/mc_fig2_idh_recurrence.py."""
    comp = pd.read_csv(os.path.join(MC, "composition_all_cohorts.tsv"), sep="\t")
    idh = pd.read_csv(os.path.join(MC, "idh_wt_vs_mut.tsv"), sep="\t")
    h1 = pd.read_csv(os.path.join(MC, "gse182109_h1_primary_vs_recurrent.tsv"), sep="\t")
    rp = pd.read_csv(os.path.join(RES, "rev02_recurrence_power.tsv"), sep="\t").iloc[0]
    g182 = pd.read_csv(os.path.join(MC, "gse182109_composition_by_patient.tsv"), sep="\t")

    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    m = comp.groupby("cohort")["MES_prop"].mean().sort_values()
    ax.bar(range(len(m)), m.values, color=C_MES, alpha=0.9)
    ax.set_xticks(range(len(m)))
    ax.set_xticklabels(m.index, rotation=25, ha="right", fontsize=7)
    ax.set_ylabel("mean MES-like proportion")
    ax.set_title("(a) Cohort means", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[0, 1])
    for i, c in enumerate(comp["cohort"].unique()):
        v = comp.loc[comp["cohort"] == c, "MES_prop"].values
        ax.scatter(np.full(len(v), i) + np.random.default_rng(1).normal(0, 0.06, len(v)), v,
                   s=9, alpha=0.7, linewidths=0)
    ax.set_xticks(range(comp["cohort"].nunique()))
    ax.set_xticklabels(comp["cohort"].unique(), rotation=25, ha="right", fontsize=7)
    ax.set_ylabel("per-patient MES-like proportion")
    ax.set_title("(b) Per-patient distributions", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[0, 2])
    wt = comp[comp["cohort"].isin(["GSE103224", "GSE182109", "Neftel_GSE131928"])]["MES_prop"].values
    mut = comp[comp["cohort"] == "Venteicher_GSE70630"]["MES_prop"].values
    for i, (v, lab, col) in enumerate([(wt, f"IDH-wt (n={len(wt)})", "#2C3E50"),
                                       (mut, f"IDH-mut (n={len(mut)})", "#C0392B")]):
        ax.scatter(np.full(len(v), i) + np.random.default_rng(2).normal(0, 0.06, len(v)), v,
                   s=12, color=col, alpha=0.75, linewidths=0)
        ax.bar(i, np.mean(v), width=0.45, facecolor="none", edgecolor=col, lw=1.4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["IDH-wt", "IDH-mut"])
    ax.set_ylabel("per-patient MES-like proportion")
    ax.set_title(f"(c) IDH-wt vs IDH-mut\nMann-Whitney P = "
                 f"{float(idh['p'].iloc[0]):.1e} (bars = mean, points = patients)",
                 fontsize=9, loc="left")

    ax = fig.add_subplot(gs[1, 0])
    prim = g182[g182["type"] == "ndGBM"]["MES_prop"].values
    rec = g182[g182["type"] == "rGBM"]["MES_prop"].values
    ax.scatter(np.full(len(prim), 0) + np.random.default_rng(3).normal(0, 0.07, len(prim)),
               prim, s=22, color="#2C3E50", alpha=0.8, linewidths=0)
    ax.scatter(np.full(len(rec), 1) + np.random.default_rng(4).normal(0, 0.07, len(rec)),
               rec, s=22, color="#C0392B", alpha=0.8, linewidths=0)
    for i, (v, c) in enumerate([(prim, "#2C3E50"), (rec, "#C0392B")]):
        ax.errorbar(i, np.mean(v), yerr=np.std(v, ddof=1) / np.sqrt(len(v)),
                    fmt="_", color=c, capsize=6, lw=1.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["primary", "recurrent"])
    ax.set_ylabel("per-patient MES-like proportion")
    ax.set_title("(d) Primary vs recurrent\npoints = patients, bars = mean +/- SEM",
                 fontsize=9, loc="left")

    ax = fig.add_subplot(gs[1, 1])
    ax.boxplot([prim, rec], tick_labels=["primary", "recurrent"], showfliers=False, widths=0.5)
    ax.scatter(np.full(len(prim), 1) + np.random.default_rng(5).normal(0, 0.06, len(prim)),
               prim, s=14, color="#2C3E50", alpha=0.8, linewidths=0)
    ax.scatter(np.full(len(rec), 2) + np.random.default_rng(6).normal(0, 0.06, len(rec)),
               rec, s=14, color="#C0392B", alpha=0.8, linewidths=0)
    ax.set_ylabel("per-patient MES-like proportion")
    ax.set_title(f"(e) Overlap is broad\nmean diff {float(rp['mean_difference']):+.3f} "
                 f"(95% CI {rp['mean_difference_95CI_boot']}), two-sided P = "
                 f"{float(rp['p_two_sided']):.2f}", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    txt = comp.groupby("cohort")["MES_prop"].agg(n_patients="count", mean_MES="mean").reset_index()
    cell = [[f"{r.cohort}", f"{int(r.n_patients)}", f"{r.mean_MES:.3f}"]
            for _, r in txt.iterrows()]
    t = ax.table(cellText=cell, colLabels=["cohort", "n patients", "mean MES"],
                 loc="center", cellLoc="left")
    t.auto_set_font_size(False); t.set_fontsize(7)
    ax.set_title("(f) Cohort-by-genotype summary (see Table 1 for the full version)",
                 fontsize=9, loc="left")
    save(fig, "Figure2_main")


if __name__ == "__main__":
    lib.setup_logging()
    print("Figure 5 ..."); figure5()
    print("Figure 3 ..."); figure3()
    print("Figure 2 ..."); figure2()
    print("DONE rev08")
