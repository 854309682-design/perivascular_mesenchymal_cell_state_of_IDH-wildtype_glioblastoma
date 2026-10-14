#!/usr/bin/env python3
"""rev19 — Figure 4 rebuilt: per-cohort DPT axis + PAGA (reviewer items 4.2, 4.6, 4.7).

Why this script exists
----------------------
Reviewer item 4.2 (MAJOR) is about Fig. 4b, whose legend gave the reader no usable colour ->
cohort mapping (the swatches were drawn ~1.4 pt, so they read as grey smudges); 4.6: the panel
labels mix `gse103224` and `GSE182109`; 4.7: the bar order in Fig. 4e differs from the legend
order in 4b. Figure 4 was produced on 2026-09-13 by a script that no longer exists anywhere on
this volume (searched 2026-10-14, including the `Datasets/Paper_3/backups/` trees), so the
figure is rebuilt from data rather than patched.

Panels and their source
-----------------------
(a) GSE103224 pooled axis          <- `Datasets/processed/gse103224_traj.h5ad`, the shipped
(c) PAGA connectivity, GSE103224      discovery-cohort axis (`dpt_pseudotime`, `score_MES`,
(f) state composition along dpt       `state`, stored kNN `connectivities`) — the same axis that
                                      `results/step04_traj/dpt_vs_MES_corr.tsv` reports, so
                                      panel (a)'s number stays reproducible from the shipped code
(b) multi-cohort axis              <- per-cohort DPT recomputed by this script for GSE182109 and
(d) per-sample rho(DPT, MES)          Neftel with the *same* procedure as the discovery cohort
(e) per-cohort rho(DPT, MES)          (`release/scripts/04_trajectory.py`, which does exist):
                                      log-normalised malignant cells -> variance-selected HVGs
                                      (2,000) -> scale -> PCA(50) -> kNN(15) -> diffusion
                                      map(15) -> DPT rooted at the cell with the highest
                                      NPC+OPC score -> oriented so pseudotime increases with
                                      the MES score. One cohort order (GSE103224, GSE182109,
                                      Neftel) and one legend per colour are used in every panel.

Neftel platform scale (data-hygiene note; declared in the caption and Methods)
-----------------------------------------------------------------------------
The two Neftel matrices in GEO are delivered on *different* scales: the Smart-seq2 file is
log2(TPM+1) (nonzero values 0.02-10) while the 10x file is raw TPM (86-691). Concatenating
them — which the submitted analysis did — makes this cohort's MES score and pseudotime axis
mostly a PLATFORM axis (measured: rho(platform, MES) = +0.82, rho(platform, dpt) = -0.90).
The cohort is therefore restricted to its Smart-seq2 block (`NEFTEL_SCALE=ss2`, the default and
the only variant used in the revised manuscript): internally consistent, and the platform that
carries that study's per-patient design. With both platforms log2-transformed
(`NEFTEL_SCALE=log10x`) a platform association remains (rho(platform, MES) = +0.43). See
scripts/rev23_neftel_ss2.py for the full before/after table.

Outputs
-------
figures/Figure4_main.pdf, figures/Figure4_main.jpeg   (JPEG + vector PDF, no grid lines)
results/revision/rev19_traj_axis_<cohort>.tsv         per-cell dpt + MES + state + sample
results/revision/rev19_traj_summary.tsv               per-cohort rho summary
results/revision/rev19_traj_per_sample.tsv            per-sample rho
results/revision/rev19_paga_gse103224.tsv             PAGA connectivity matrix

Run:
  GLIOMA_PROJ=/run/media/meng/Seqdata/glioma/Paper_2 \
  GLIOMA_DATA=/run/media/meng/Seqdata/glioma/Datasets \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev19_figure4_trajectory.py
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad_mod
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules, four_state_tables
import rev01_export_scenic as r1

sc.settings.verbosity = 0
RES = os.path.join(lib.PROJ, "results", "revision")
FIG = os.path.join(lib.PROJ, "figures")
TRAJ = os.path.join(lib.PROC_DIR, "gse103224_traj.h5ad")
os.makedirs(RES, exist_ok=True)

SEED = 42
NHVG = 2000
NEW_COHORTS = [("GSE182109", "GSE182109"), ("Neftel_GSE131928", "Neftel")]
COHORT_ORDER = ["GSE103224", "GSE182109", "Neftel"]
COHORT_COLOR = {"GSE103224": "#C0392B", "GSE182109": "#2980B9", "Neftel": "#27AE60"}
STATE_COLOR = {"AC": "#7F8C8D", "NPC": "#2980B9", "OPC": "#27AE60", "MES": "#C0392B"}
PAGA_ORDER = ["AC", "MES", "NPC", "OPC"]
NEFTEL_SCALE = os.environ.get("NEFTEL_SCALE", "ss2")         # "ss2" (used) | "log10x"

plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": False, "font.family": "DejaVu Sans"})


# ------------------------------------------------------------------ cohort matrices
def neftel_matrix():
    """Both Neftel platforms on the log2(TPM+1) scale (or the Smart-seq2 block only)."""
    a1 = r1.load_log_matrix(lib.NEFTEL_TPM_SS2, "Neftel_GSE131928",
                            lambda c: c.split("-")[0], "Smartseq2")
    if NEFTEL_SCALE == "ss2":
        lib.log.info("[rev19] Neftel: Smart-seq2 block only")
        return a1
    a2 = r1.load_log_matrix(lib.NEFTEL_TPM_10X, "Neftel_GSE131928",
                            lambda c: c.split("_")[0], "10x")
    if sp.issparse(a2.X):                        # log2(TPM+1) in place; stored zeros stay 0
        a2.X = a2.X.copy()
        a2.X.data = np.log2(a2.X.data + 1.0).astype(np.float32)
    else:
        a2.X = np.log2(np.asarray(a2.X, dtype=np.float32) + 1.0)
    a = ad_mod.concat([a1, a2], join="outer")
    a.var_names_make_unique()
    lib.log.info("[rev19] Neftel: SS2 + 10x, both log2(TPM+1) (scale-corrected)")
    return a


def malignant_cells(name, sub):
    """Malignant log-normalised cells for one cohort (Neftel via the corrected matrix)."""
    if name == "Neftel_GSE131928":
        a = r1.score_states(neftel_matrix(), sub)
        mal = a[a.obs["nonmalignant"] == 0].copy()
        lib.log.info(f"[rev19] [{name}] {a.n_obs} cells -> {mal.n_obs} malignant")
        return mal
    return r1.load_cohort(name, sub)


def build_axis(mal, label):
    """DPT axis, mirroring release/scripts/04_trajectory.py. Returns (AnnData, raw rho)."""
    X = mal.X
    if not sp.issparse(X):
        X = sp.csr_matrix(np.asarray(X, dtype=np.float32))
        mal.X = X
    mu = np.asarray(X.mean(axis=0)).ravel()
    var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mu ** 2
    keep = np.sort(np.argsort(var)[::-1][:NHVG])
    mal = mal[:, keep].copy()
    lib.log.info(f"[rev19] [{label}] {mal.n_obs} cells x {mal.n_vars} variance-HVGs")
    sc.pp.scale(mal, max_value=10)
    sc.tl.pca(mal, n_comps=50, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(mal, n_neighbors=15, random_state=SEED)
    sc.tl.diffmap(mal, n_comps=15, random_state=SEED)
    prog = mal.obs["score_NPC"].values + mal.obs["score_OPC"].values
    mal.uns["iroot"] = int(np.argmax(prog))
    lib.log.info(f"[rev19] [{label}] root (max NPC+OPC) = {mal.obs_names[mal.uns['iroot']]}")
    sc.tl.dpt(mal, n_dcs=15)
    raw = float(spearmanr(mal.obs["dpt_pseudotime"].values,
                          mal.obs["score_MES"].values).statistic)
    if raw < 0:                                  # orient so pseudotime increases with MES
        mal.obs["dpt_pseudotime"] = -mal.obs["dpt_pseudotime"].values
    lib.log.info(f"[rev19] [{label}] rho(dpt,MES) raw {raw:+.3f} -> |rho| {abs(raw):.3f}")
    return mal, raw


def as_frame(mal, label):
    return pd.DataFrame({"cell": mal.obs_names,
                         "sample": mal.obs["patient"].astype(str).values
                         if "patient" in mal.obs else mal.obs["sample"].astype(str).values,
                         "state": mal.obs["state"].astype(str).values,
                         "mes": np.asarray(mal.obs["score_MES"], dtype=float),
                         "dpt": np.asarray(mal.obs["dpt_pseudotime"], dtype=float)})


def per_sample_rho(df, cohort, min_cells=30):
    out = []
    for s, g in df.groupby("sample", observed=True):
        if len(g) >= min_cells and g["dpt"].nunique() > 5:
            out.append(dict(cohort=cohort, sample=s, n_cells=len(g),
                            rho=float(spearmanr(g["dpt"], g["mes"]).statistic)))
    return pd.DataFrame(out)


# ------------------------------------------------------------------ figure
def draw(frames, conn, summ, per_sample):
    rho = dict(zip(summ["cohort"], summ["rho_pooled"]))
    med = dict(zip(summ["cohort"], summ["rho_median_per_sample"]))
    fig = plt.figure(figsize=(13, 9.2))
    gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.30)

    # (a) discovery-cohort axis
    d = frames["GSE103224"]
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(d["dpt"], d["mes"], s=3, alpha=0.25, c=d["mes"], cmap="magma", linewidths=0)
    ax.set_xlabel("dpt pseudotime (MES-oriented)"); ax.set_ylabel("MES module score")
    ax.set_title(f"(a) GSE103224 pooled dpt vs MES\nSpearman \u03c1 = {rho['GSE103224']:.2f}",
                 fontsize=9, loc="left")

    # (b) multi-cohort axis, z-scored within cohort
    ax = fig.add_subplot(gs[0, 1])
    handles = []
    for c in COHORT_ORDER:
        dd = frames[c]
        x = (dd["dpt"] - dd["dpt"].mean()) / dd["dpt"].std()
        y = (dd["mes"] - dd["mes"].mean()) / dd["mes"].std()
        ax.scatter(x, y, s=2.5, alpha=0.22, color=COHORT_COLOR[c], linewidths=0)
        handles.append(Line2D([], [], marker="o", ls="", ms=7, color=COHORT_COLOR[c],
                              label=f"{c}: pooled \u03c1 = {rho[c]:.2f} "
                                    f"(median per sample {med[c]:.2f})"))
    ax.set_xlabel("oriented pseudotime (z per cohort)")
    ax.set_ylabel("MES score (z per cohort)")
    ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="upper left")
    ax.set_title("(b) Per-cohort axis, z-scored per cohort", fontsize=9, loc="left")

    # (c) PAGA heatmap
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(conn.values, cmap="viridis", vmin=0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(PAGA_ORDER)))
    ax.set_xticklabels(PAGA_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(PAGA_ORDER))); ax.set_yticklabels(PAGA_ORDER)
    for i in range(len(PAGA_ORDER)):
        for j in range(len(PAGA_ORDER)):
            ax.text(j, i, f"{conn.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if conn.values[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="PAGA connectivity")
    ax.set_title("(c) PAGA state connectivity (GSE103224)", fontsize=9, loc="left")

    # (d) per-sample rho
    ax = fig.add_subplot(gs[1, 0])
    data = [per_sample.loc[per_sample["cohort"] == c, "rho"].values for c in COHORT_ORDER]
    ax.boxplot(data, tick_labels=COHORT_ORDER, widths=0.6, showfliers=False)
    for i, v in enumerate(data, start=1):
        ax.scatter(np.full(len(v), i) + np.random.default_rng(i).normal(0, 0.05, len(v)), v,
                   s=10, color=COHORT_COLOR[COHORT_ORDER[i - 1]], linewidths=0, zorder=3)
    ax.axhline(0, color="grey", ls=":", lw=0.8)
    ax.set_ylabel("\u03c1(DPT, MES) per sample")
    ax.set_title("(d) Per-sample correlation by cohort", fontsize=9, loc="left")

    # (e) per-cohort rho
    ax = fig.add_subplot(gs[1, 1])
    ax.bar(range(len(COHORT_ORDER)), [rho[c] for c in COHORT_ORDER],
           color=[COHORT_COLOR[c] for c in COHORT_ORDER], width=0.62)
    for i, c in enumerate(COHORT_ORDER):
        ax.text(i, rho[c] + 0.012, f"{rho[c]:.2f}", ha="center", fontsize=8)
    for i, c in enumerate(COHORT_ORDER):
        ax.plot([i - 0.31, i + 0.31], [med[c], med[c]], color="black", lw=1.6, zorder=4)
    ax.plot([], [], color="black", lw=1.6, label="median per sample")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_ylim(0, max(rho[c] for c in COHORT_ORDER) * 1.30)
    ax.set_xticks(range(len(COHORT_ORDER))); ax.set_xticklabels(COHORT_ORDER, fontsize=8)
    ax.set_ylabel("\u03c1(oriented DPT, MES)")
    ax.set_title("(e) Pooled per-cohort correlation\n(black line = median per sample)",
                 fontsize=9, loc="left")

    # (f) state composition along pseudotime
    ax = fig.add_subplot(gs[1, 2])
    d = frames["GSE103224"]
    # quantile (equal-count) bins: the stored discovery-cohort pseudotime is heavily skewed
    # (median 0.057, 95th percentile 0.22, max 1.0), so equal-width bins would leave >95% of
    # cells in the first three bins and put a handful of cells in the rest
    idx = np.clip(np.digitize(d["dpt"], np.quantile(d["dpt"], np.linspace(0, 1, 11)[1:-1])),
                  0, 9)
    rows = []
    for st in ["AC", "NPC", "OPC", "MES"]:
        frac = [float((d.loc[idx == b, "state"] == st).mean()) for b in range(10)]
        for b, f in enumerate(frac):
            rows.append(dict(bin=b, state=st, fraction=f))
        ax.plot(range(10), frac, "-o", ms=3, lw=1.3, color=STATE_COLOR[st], label=st)
    bdf = pd.DataFrame(rows)
    bdf.to_csv(os.path.join(RES, "rev19_bins_gse103224.tsv"), sep="\t", index=False)
    piv = bdf.pivot(index="bin", columns="state", values="fraction")
    print(piv.round(3).to_string())
    lib.log.info(f"[rev19] MES fraction peaks in bin {int(piv['MES'].idxmax())} "
                 f"({piv['MES'].max():.3f}); AC peaks in bin {int(piv['AC'].idxmax())} "
                 f"({piv['AC'].max():.3f}); NPC peaks in bin {int(piv['NPC'].idxmax())} "
                 f"({piv['NPC'].max():.3f})")
    ax.set_xlabel("DPT decile"); ax.set_ylabel("state fraction")
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    ax.set_title("(f) State composition along pseudotime (GSE103224)", fontsize=9, loc="left")

    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"Figure4_main.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote Figure4_main.pdf/.jpeg")


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev19.log"))
    sub = four_state_tables(load_neftel_modules())[1]

    frames, summary, persample = {}, [], []

    # (1) discovery cohort — the shipped axis, used unchanged
    a = sc.read_h5ad(TRAJ)
    d = pd.DataFrame({"cell": a.obs_names, "sample": a.obs["sample"].astype(str).values,
                      "state": a.obs["state"].astype(str).values,
                      "mes": np.asarray(a.obs["score_MES"], dtype=float),
                      "dpt": np.asarray(a.obs["dpt_pseudotime"], dtype=float)})
    rho = float(spearmanr(d["dpt"], d["mes"]).statistic)
    if rho < 0:
        d["dpt"] = -d["dpt"]
        rho = -rho
    lib.log.info(f"[rev19] [GSE103224] shipped axis rho(dpt,MES) = {rho:+.3f} (n={len(d)})")
    frames["GSE103224"] = d
    d.to_csv(os.path.join(RES, "rev19_traj_axis_GSE103224.tsv"), sep="\t", index=False)
    ps = per_sample_rho(d, "GSE103224")
    persample.append(ps)
    summary.append(dict(cohort="GSE103224", n_cells=len(d), n_samples=len(ps), rho_pooled=rho,
                        rho_median_per_sample=float(ps["rho"].median()),
                        n_samples_rho_pos=int((ps["rho"] > 0).sum()),
                        n_samples_rho_gt_0_3=int((ps["rho"] > 0.3).sum()),
                        source="shipped step04 axis (gse103224_traj.h5ad)"))

    # (2) PAGA on the discovery cohort (stored kNN graph + four state labels)
    agg = a.copy()
    sc.tl.paga(agg, groups="state")
    pc = agg.uns["paga"]["connectivities"]
    # uns["paga"]["groups"] is the *key name* ("state"); the matrix is ordered by the
    # categorical's categories
    pg = list(agg.obs["state"].cat.categories)
    if sp.issparse(pc):
        pc = np.asarray(pc.todense())
    conn = pd.DataFrame(np.asarray(pc), index=pg, columns=pg)
    conn = conn.reindex(index=PAGA_ORDER, columns=PAGA_ORDER)
    conn.to_csv(os.path.join(RES, "rev19_paga_gse103224.tsv"), sep="\t")
    lib.log.info("[rev19] PAGA connectivity:\n" + conn.round(3).to_string())

    # (3) the two cohorts whose axes are recomputed here
    for name, label in NEW_COHORTS:
        mal, raw = build_axis(malignant_cells(name, sub), label)
        d = as_frame(mal, label)
        frames[label] = d
        d.to_csv(os.path.join(RES, f"rev19_traj_axis_{label}.tsv"), sep="\t", index=False)
        ps = per_sample_rho(d, label)
        persample.append(ps)
        summary.append(dict(cohort=label, n_cells=len(d), n_samples=len(ps),
                            rho_pooled=float(spearmanr(d["dpt"], d["mes"]).statistic),
                            rho_median_per_sample=float(ps["rho"].median()),
                            n_samples_rho_pos=int((ps["rho"] > 0).sum()),
                            n_samples_rho_gt_0_3=int((ps["rho"] > 0.3).sum()),
                            source=f"rev19 ({'Neftel ' + NEFTEL_SCALE})"))

    summ = pd.DataFrame(summary)
    ps_all = pd.concat(persample, ignore_index=True)
    summ.to_csv(os.path.join(RES, "rev19_traj_summary.tsv"), sep="\t", index=False)
    ps_all.to_csv(os.path.join(RES, "rev19_traj_per_sample.tsv"), sep="\t", index=False)
    print(summ.to_string(index=False))
    draw(frames, conn, summ, ps_all)
    lib.log.info("DONE rev19")


if __name__ == "__main__":
    main()
