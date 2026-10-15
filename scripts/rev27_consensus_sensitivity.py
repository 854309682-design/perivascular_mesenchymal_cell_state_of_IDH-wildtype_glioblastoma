#!/usr/bin/env python3
"""rev27 — Prior-route consensus: leave-one-cohort-out sensitivity (reviewer minor 5).

The reviewer's point: GSE89567 (Venteicher 2017) has a mean MES-like proportion of ~0.01, so
correlating per-cell TF activity against an almost constant score gives unstable ρ estimates,
and including such a cohort in a "4/4 cohorts" claim dilutes what that claim means. They asked
us either to justify the inclusion (e.g. with a variance floor) or to report the consensus with
and without it.

This script does the second, from the re-implemented prior-route tables:

  * the consensus (mean ρ over cohorts, and the number of cohorts with FDR < 0.05) computed on
    all four cohorts, on the three IDH-wildtype cohorts only, and leave-one-cohort-out;
  * the rank of every pre-specified key factor in each of those variants, so a reader can see
    exactly which conclusions move;
  * the spread of the MES score per cohort (SD and IQR), which is the quantity the review
    worried about.

Outputs (results/revision/):
    rev27_consensus_sensitivity.tsv      consensus rank + mean ρ of the key factors per variant
    rev27_key_factor_ranks.tsv           full key-factor table per variant
    rev27_mes_score_spread.tsv           per-cohort MES-score SD/IQR (variance-floor evidence)

Run:
  GLIOMA_PROJ=/run/media/meng/Seqdata/glioma/Paper_2 \
  GLIOMA_DATA=/run/media/meng/Seqdata/glioma/Datasets \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev27_consensus_sensitivity.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

RES = os.path.join(lib.PROJ, "results", "revision")
OUT = RES
KEY = ["STAT3", "KLF4", "SP1", "FOSL1", "FOSL2", "JUNB", "CEBPD", "CEBPB", "RELA", "TWIST1",
       "TEAD1", "HIF1A"]
IDH_MUT = ["Venteicher", "GSE70630", "GSE89567"]
MES_CELLS = {  # per-cell MES scores, one file per cohort (whichever the cohort's run wrote)
    "GSE103224": os.path.join(lib.PROJ, "results", "step03_state", "state_scores_per_cell.tsv"),
    "GSE182109": os.path.join(lib.PROJ, "results", "multicohort", "gse182109_state_scores.tsv"),
    "Neftel": os.path.join(lib.PROJ, "results", "multicohort",
                           "Neftel_GSE131928_state_scores.tsv"),
    "GSE70630": os.path.join(lib.PROJ, "results", "multicohort",
                             "Venteicher_GSE70630_state_scores.tsv"),
    "GSE89567": os.path.join(lib.PROJ, "results", "multicohort", "gse89567_state_scores.tsv"),
}


def consensus(rho, cohorts):
    """Mean ρ and the number of cohorts with FDR < 0.05, over the given cohorts."""
    sub = rho[rho["cohort"].isin(cohorts)]
    g = sub.groupby("tf").agg(mean_rho=("rho", "mean"),
                              cohorts=("cohort", "nunique"),
                              n_fdr005=("fdr", lambda x: int((x < 0.05).sum())))
    return g.sort_values("mean_rho", ascending=False)


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev27.log"))
    rho = pd.read_csv(os.path.join(RES, "rev02_tf_activity_vs_MES.tsv"), sep="\t")
    rho = rho[["cohort", "tf", "rho", "fdr"]].dropna()

    variants = {"all four cohorts": sorted(rho["cohort"].unique())}
    variants["IDH-wildtype cohorts only"] = [c for c in sorted(rho["cohort"].unique())
                                             if c not in IDH_MUT]
    for c in sorted(rho["cohort"].unique()):
        variants[f"without {c}"] = [x for x in sorted(rho["cohort"].unique()) if x != c]

    rows, ranks = [], {}
    for name, cohorts in variants.items():
        cons = consensus(rho, cohorts)
        rk = {tf: int(np.where(cons.index == tf)[0][0]) + 1 for tf in KEY if tf in cons.index}
        ranks[name] = rk
        rows.append(dict(variant=name, n_cohorts=len(cohorts), cohorts=", ".join(cohorts),
                         top5=", ".join(cons.index[:5]),
                         **{f"rank_{tf}": rk.get(tf, np.nan) for tf in KEY[:6]}))
        lib.log.info(f"[rev27] {name}: top-5 = {list(cons.index[:5])}; "
                     f"STAT3 rank {rk.get('STAT3')}, FOSL1 rank {rk.get('FOSL1')}, "
                     f"JUNB rank {rk.get('JUNB')}, CEBPB rank {rk.get('CEBPB')}, "
                     f"TWIST1 rank {rk.get('TWIST1')}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "rev27_consensus_sensitivity.tsv"),
                              sep="\t", index=False)

    cons = consensus(rho, variants["all four cohorts"])
    kt = pd.DataFrame([{"factor": tf,
                        "mean_rho_all4": float(cons.loc[tf, "mean_rho"]),
                        "n_fdr005_all4": int(cons.loc[tf, "n_fdr005"]),
                        **{f"rank_{v}": ranks[v].get(tf, np.nan) for v in variants},
                        "sign_stable": len({np.sign(cons.loc[tf, "mean_rho"])} |
                                           {np.sign(consensus(rho, c).loc[tf, "mean_rho"])
                                            for c in variants.values() if tf in
                                            consensus(rho, c).index}) == 1}
                       for tf in KEY if tf in cons.index])
    kt.to_csv(os.path.join(OUT, "rev27_key_factor_ranks.tsv"), sep="\t", index=False)

    spread = []
    for c, p in MES_CELLS.items():
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p, sep="\t")
        col = "score_MES" if "score_MES" in d.columns else "mes"
        v = pd.to_numeric(d[col], errors="coerce").dropna()
        spread.append(dict(cohort=c, n_cells=len(v), sd=round(float(v.std()), 3),
                           iqr=round(float(v.quantile(.75) - v.quantile(.25)), 3),
                           mean=round(float(v.mean()), 3)))
    sp = pd.DataFrame(spread)
    sp.to_csv(os.path.join(OUT, "rev27_mes_score_spread.tsv"), sep="\t", index=False)

    print("\nconsensus sensitivity:")
    print(pd.DataFrame(rows)[["variant", "n_cohorts", "top5"]].to_string(index=False))
    print("\nkey-factor ranks (1 = top of the consensus):")
    print(kt[["factor", "mean_rho_all4", "n_fdr005_all4"] +
             [c for c in kt.columns if c.startswith("rank_")] + ["sign_stable"]].to_string(index=False))
    print("\nMES-score spread per cohort:")
    print(sp.to_string(index=False))
    lib.log.info("DONE rev27")


if __name__ == "__main__":
    main()
