#!/usr/bin/env python3
"""rev07 — De-novo (pySCENIC) arm: consensus ranking and reconciliation with the prior route.

Reviewer items M2 / DA2 / DA3 (all CRITICAL):
  (a) the de-novo top list printed in Results §3.3 did not match Fig. 3b (NPAS2, CD59,
      NFATC4, FOXC1, TPI1 appeared in the figure; BHLHE40 and JUNB appeared in the text);
  (b) CD59 and TPI1 are not transcription factors, so their presence as "regulons" was an
      artefact that had to be explained or removed;
  (c) the claimed convergence of the two inference routes was never quantified (a 4/30
      overlap was described as "share exactly these factors");
  (d) TWIST1 rested on one route only;
  (e) CEBPB flipped sign between the two routes and that was not addressed.

This script recomputes everything from the *reproducible* pySCENIC run (rev01/rev03) and
prints the reconciliation that the revised manuscript reports:

  * per-cohort regulon-vs-MES Spearman rho from AUCell activity;
  * the three-cohort consensus ranking, with the number of cohorts reaching FDR < 0.05;
  * a TF-anchoring audit of every ranked regulon against the cisTarget TF list
    (allTFs_hg38.txt) — this is what decides (b) objectively;
  * the top-30 overlap with the prior-based route as a Fisher exact test against the
    scored-TF universe — (c);
  * TWIST1's rank in each route — (d);
  * CEBPB's rho in each route per cohort — (e).

Outputs (results/revision/):
    rev07_denovo_activity_<cohort>.tsv   regulon x cell activity summary (optional full matrix)
    rev07_denovo_vs_MES.tsv              per-cohort per-regulon rho + FDR
    rev07_denovo_consensus.tsv           mean rho across cohorts, n cohorts FDR<0.05
    rev07_tf_anchoring_audit.tsv         every ranked "regulon" vs the cisTarget TF list
    rev07_route_overlap.tsv              top-30 overlap + Fisher exact test
    rev07_key_factors.tsv                TWIST1 / CEBPB / AP-1 / C/EBP ranks in both routes
    rev07_summary.txt                    human-readable reconciliation summary

Run:
  GLIOMA_PROJ=... GLIOMA_DATA=... \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev07_denovo_reconcile.py
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(__file__))
import lib

OUT = os.path.join(lib.PROJ, "results", "revision")
SCENIC = os.path.join(lib.PROJ, "results", "scenic")
os.makedirs(OUT, exist_ok=True)
KEY = ["TWIST1", "CEBPB", "CEBPD", "CEBPA", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
       "STAT3", "KLF4", "RELA", "RELB", "NFKB1", "BHLHE40", "ATF3", "ATF4", "HIF1A",
       "TEAD1", "TEAD4", "CD59", "TPI1", "NPAS2", "FOXC1", "NFATC4"]


def tf_universe():
    p = os.path.join(lib.DATASETS, "scenic_resources", "allTFs_hg38.txt")
    with open(p) as f:
        return {l.strip() for l in f if l.strip()}


def load_aucell(cohort):
    p = os.path.join(SCENIC, f"rev03_{cohort}_aucell.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p, index_col=0)
    df.columns = [c.replace("(+)", "").strip() for c in df.columns]
    return df


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev07_denovo.log"))
    tfs = tf_universe()
    lib.log.info(f"cisTarget TF universe: {len(tfs)} TFs")

    rho_rows, audit_rows, per_cohort = [], [], {}
    for cohort in ["GSE103224", "GSE182109", "Neftel_GSE131928", "GSE70630"]:
        act = load_aucell(cohort)
        meta_p = os.path.join(SCENIC, f"rev01_{cohort}_cellmeta.tsv")
        if act is None or not os.path.exists(meta_p):
            lib.log.info(f"[{cohort}] not run / no output yet — skipped")
            continue
        meta = pd.read_csv(meta_p, sep="\t", index_col=0)
        common = act.index.intersection(meta.index)
        act = act.loc[common]
        mes = meta.loc[common, "score_MES"].values
        lib.log.info(f"[{cohort}] {act.shape[0]} cells x {act.shape[1]} regulons; "
                     f"regulons that are not in the cisTarget TF list: "
                     f"{sorted(set(act.columns) - tfs)}")
        rows = []
        for r in act.columns:
            v = act[r].values
            if np.std(v) < 1e-12:
                continue
            rho, p = spearmanr(v, mes)
            rows.append({"cohort": cohort, "regulon": r, "rho": rho, "p": p,
                         "cell_level_p_is_descriptive": True})
        d = pd.DataFrame(rows)
        d["fdr"] = multipletests(d["p"].values, method="fdr_bh")[1]
        d["is_cisTarget_TF"] = d["regulon"].isin(tfs)
        d.to_csv(os.path.join(OUT, f"rev07_denovo_vs_MES_{cohort}.tsv"), sep="\t", index=False)
        rho_rows.append(d)
        per_cohort[cohort] = act
        for r in act.columns:
            audit_rows.append({"cohort": cohort, "regulon": r, "is_cisTarget_TF": r in tfs})

    if not rho_rows:
        lib.log.warning("no AUCell output yet — rerun after rev03 finishes")
        return
    allrho = pd.concat(rho_rows, ignore_index=True)
    allrho.to_csv(os.path.join(OUT, "rev07_denovo_vs_MES.tsv"), sep="\t", index=False)
    pd.DataFrame(audit_rows).drop_duplicates().to_csv(
        os.path.join(OUT, "rev07_tf_anchoring_audit.tsv"), sep="\t", index=False)

    cons = allrho.groupby("regulon").agg(cohorts=("cohort", "nunique"),
                                         mean_rho=("rho", "mean"),
                                         n_fdr005=("fdr", lambda s: int((s < 0.05).sum())),
                                         is_cisTarget_TF=("is_cisTarget_TF", "max")).reset_index()
    cons = cons.sort_values("mean_rho", ascending=False)
    cons.to_csv(os.path.join(OUT, "rev07_denovo_consensus.tsv"), sep="\t", index=False)

    # ---- top-30 overlap with the prior-based route + Fisher exact test ----
    prior = pd.read_csv(os.path.join(OUT, "rev02_tf_consensus_reimplemented.tsv"), sep="\t")
    prior_top = set(prior.sort_values("mean_rho", ascending=False).head(30)["tf"])
    dno_top = set(cons.head(30)["regulon"])
    a = len(prior_top & dno_top)
    b = len(prior_top - dno_top)
    c = len(dno_top - prior_top)
    universe = set(prior["tf"]) | set(cons["regulon"])
    d = len(universe - prior_top - dno_top)
    odds, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")
    shared = sorted(prior_top & dno_top)
    ov = pd.DataFrame([{"prior_top30": len(prior_top), "denovo_top30": len(dno_top),
                        "shared": a, "shared_factors": ", ".join(shared),
                        "prior_only": b, "denovo_only": c, "universe": len(universe),
                        "odds_ratio": round(float(odds), 2), "fisher_p_greater": float(p_fisher)}])
    ov.to_csv(os.path.join(OUT, "rev07_route_overlap.tsv"), sep="\t", index=False)

    # ---- key factors in both routes ----
    key_rows = []
    for tf in KEY:
        pr = prior[prior["tf"] == tf]
        dn = cons[cons["regulon"] == tf]
        key_rows.append({
            "factor": tf,
            "prior_mean_rho": round(float(pr["mean_rho"].iloc[0]), 3) if len(pr) else np.nan,
            "prior_n_fdr005_cohorts": int(pr["n_fdr005"].iloc[0]) if len(pr) else np.nan,
            "denovo_mean_rho": round(float(dn["mean_rho"].iloc[0]), 3) if len(dn) else np.nan,
            "denovo_n_fdr005_cohorts": int(dn["n_fdr005"].iloc[0]) if len(dn) else np.nan,
            "denovo_is_cisTarget_TF": bool(dn["is_cisTarget_TF"].iloc[0]) if len(dn) else np.nan,
            "sign_agreement": ("agree" if (len(pr) and len(dn) and
                                           np.sign(pr["mean_rho"].iloc[0]) == np.sign(dn["mean_rho"].iloc[0]))
                               else ("DISAGREE" if (len(pr) and len(dn)) else "na")),
        })
    kt = pd.DataFrame(key_rows)
    kt["denovo_rank"] = kt["factor"].map(
        {r: i + 1 for i, r in enumerate(cons["regulon"])} if "regulon" in cons else {})
    kt.to_csv(os.path.join(OUT, "rev07_key_factors.tsv"), sep="\t", index=False)

    lines = []
    lines.append(f"de-novo arm: cohorts analysed = {sorted(allrho['cohort'].unique())}")
    lines.append(f"regulons per cohort: " + ", ".join(
        f"{c}={per_cohort[c].shape[1]}" for c in per_cohort))
    bad = sorted(set(audit_rows and [a['regulon'] for a in audit_rows if not a['is_cisTarget_TF']]))
    lines.append(f"regulons NOT in the cisTarget TF list (would be artefacts): {bad if bad else 'none'}")
    lines.append(f"CD59/TPI1 among ranked regulons: CD59={('CD59' in set(cons['regulon']))}, "
                 f"TPI1={('TPI1' in set(cons['regulon']))}")
    lines.append(f"top-30 overlap with the prior route: {a}/{len(prior_top)} shared "
                 f"({', '.join(shared)}); Fisher exact P = {p_fisher:.3g}")
    lines.append("top-10 de-novo regulons:")
    for i, r in enumerate(cons.head(10).itertuples(), 1):
        lines.append(f"   {i:2d}. {r.regulon:10s} mean rho={r.mean_rho:+.3f}  "
                     f"FDR<0.05 in {r.n_fdr005}/{r.cohorts} cohorts  TF={r.is_cisTarget_TF}")
    for tf in ["TWIST1", "CEBPB", "CEBPD", "STAT3"]:
        rr = kt[kt["factor"] == tf]
        if len(rr):
            r = rr.iloc[0]
            lines.append(f"{tf:7s} prior rho={r.prior_mean_rho} ({r.prior_n_fdr005_cohorts} cohorts FDR<0.05)"
                         f" | de-novo rho={r.denovo_mean_rho} ({r.denovo_n_fdr005_cohorts} cohorts)"
                         f" | {r.sign_agreement}")
    txt = "\n".join(lines)
    open(os.path.join(OUT, "rev07_summary.txt"), "w").write(txt + "\n")
    lib.log.info("RECONCILIATION:\n" + txt)
    lib.log.info("DONE rev07")


if __name__ == "__main__":
    main()
