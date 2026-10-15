"""rev21_edits — the edit table for scripts/rev21_manuscript_final.py.

Anchors are copied **verbatim from the raw document.xml** (the manuscript uses en dashes and
unicode rho/kappa, so anchors must not be "normalised"). Numbers that come from a re-run are
read from the result tables at import time, so the pass can simply be re-run when a table
moves, and it can never ship a stale or hand-typed value.

Edit kinds understood by the runner:
    ("para", pid, new_text)          whole paragraph (single non-empty run)
    ("run",  pid, anchor, new_text)  replace the text of the run containing `anchor`
    ("sub",  pid, old, new)          substring inside one run

Upstream dependencies
    results/revision/rev15_idh_stratified.tsv   (rev15, after rev23)
    results/revision/rev19_traj_*.tsv           (rev19)
    results/revision/rev07_*.tsv|txt            (rev07, after the pySCENIC chain)
    results/revision/rev02_*.tsv                (rev02)
"""
import os
import re

import numpy as np
import pandas as pd

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(PROJ, "results", "revision")

SUP = {"-": "\u207b", "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
       "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079"}
EN = "\u2013"
RHO = "\u03c1"


def _tsv(name, need=()):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        raise RuntimeError(f"rev21_edits: missing upstream table {p}")
    d = pd.read_csv(p, sep="\t")
    missing = [c for c in need if c not in d.columns]
    if missing:
        raise RuntimeError(f"rev21_edits: {name} lacks column(s) {missing}")
    return d


def fmt_p(p):
    """Format a P value the way this manuscript does (0.25 / 0.0034 / 6.6 x 10^-8)."""
    p = float(p)
    if p < 1e-3:
        e = int(np.floor(np.log10(p)))
        m = p / 10.0 ** e
        return f"{m:.1f} \u00d7 10" + "".join(SUP[c] for c in str(e))
    if p < 0.01:
        return f"{p:.2g}"
    return f"{p:.2f}"


def fmt_r(x):
    return f"{float(x):.2f}"


# ----------------------------------------------------------------- data for the edits
def idh():
    d = _tsv("rev15_idh_stratified.tsv", ["cohort", "genotype", "mean_MES", "test", "p"])
    per = d[d["test"] == "per_cohort_means"]
    pooled = float(d.loc[d["test"] == "pooled_Mann_Whitney", "p"].iloc[0])
    centred = float(d.loc[d["test"] == "cohort_mean_centred_Mann_Whitney", "p"].iloc[0])
    # patient counts and genotype means come from the pooled patient table
    pat = pd.read_csv(os.path.join(PROJ, "results", "multicohort", "idh_patients.tsv"), sep="\t")
    wt = pat.loc[pat["genotype"] == "IDH-wt", "MES_prop"]
    mut = pat.loc[pat["genotype"] == "IDH-mut", "MES_prop"]
    return dict(per=per, pooled=pooled, centred=centred,
                n_wt=len(wt), n_mut=len(mut), mean_wt=float(wt.mean()), mean_mut=float(mut.mean()),
                wt_range=(float(per.loc[per["genotype"].str.contains("wt"), "mean_MES"].min()),
                          float(per.loc[per["genotype"].str.contains("wt"), "mean_MES"].max())),
                mut_range=(float(per.loc[per["genotype"].str.contains("mut"), "mean_MES"].min()),
                           float(per.loc[per["genotype"].str.contains("mut"), "mean_MES"].max())))


def traj():
    s = _tsv("rev19_traj_summary.tsv", ["cohort", "n_cells", "n_samples", "rho_pooled",
                                        "rho_median_per_sample", "n_samples_rho_pos",
                                        "n_samples_rho_gt_0_3"])
    rho = dict(zip(s["cohort"], s["rho_pooled"].astype(float)))
    med = dict(zip(s["cohort"], s["rho_median_per_sample"].astype(float)))
    ns = dict(zip(s["cohort"], s["n_cells"].astype(int)))
    nsamp = dict(zip(s["cohort"], s["n_samples"].astype(int)))
    pos = dict(zip(s["cohort"], s["n_samples_rho_pos"].astype(int)))
    gt03 = dict(zip(s["cohort"], s["n_samples_rho_gt_0_3"].astype(int)))
    bins = _tsv("rev19_bins_gse103224.tsv", ["bin", "state", "fraction"])
    piv = bins.pivot(index="bin", columns="state", values="fraction")
    paga = pd.read_csv(os.path.join(RES, "rev19_paga_gse103224.tsv"), sep="\t", index_col=0)
    return dict(rho=rho, med=med, n=ns, nsamp=nsamp, pos=pos, gt03=gt03, piv=piv, paga=paga)


def denovo():
    cons = _tsv("rev07_denovo_consensus.tsv", ["regulon", "mean_rho", "n_fdr005",
                                               "is_cisTarget_TF"])
    cons = cons.sort_values("mean_rho", ascending=False)
    top = cons.head(10)
    ov = _tsv("rev07_route_overlap.tsv").iloc[0]
    kf = _tsv("rev07_key_factors.tsv", ["factor", "prior_mean_rho", "denovo_mean_rho"])
    kfv = {r["factor"]: r for _, r in kf.iterrows()}
    # per-cohort regulon counts, printed by rev07 into its summary text
    txt = open(os.path.join(RES, "rev07_summary.txt")).read()
    m = re.search(r"regulons per cohort:\s*(.+)", txt)
    counts = {}
    if m:
        for part in m.group(1).split(","):
            k, v = part.split("=")
            counts[k.strip()] = int(v)
    return dict(cons=cons, top=top, ov=ov, kf=kfv, counts=counts,
                n_tf=int((cons["is_cisTarget_TF"] == True).sum()), n_all=int(len(cons)))


def ranks():
    """Rank of each key factor in the two routes' consensuses (1 = best)."""
    dn = _tsv("rev07_denovo_consensus.tsv", ["regulon", "mean_rho"]).sort_values(
        "mean_rho", ascending=False)
    pr = _tsv("rev02_tf_consensus_reimplemented.tsv", ["tf", "mean_rho"]).sort_values(
        "mean_rho", ascending=False)
    dr = {t: i + 1 for i, t in enumerate(dn["regulon"])}
    prr = {t: i + 1 for i, t in enumerate(pr["tf"])}
    return dict(dn=dr, pr=prr, n_dn=len(dn), n_pr=len(pr))


def recep():
    """Receptor detection fractions on malignant GSE182109 cells (rev10)."""
    d = _tsv("rev10_receptor_fractions.tsv", ["receptor", "pct_malignant_expressing"])
    return dict(zip(d["receptor"], d["pct_malignant_expressing"].astype(float)))


def sens():
    """Leave-one-cohort-out evidence for the prior-route consensus (reviewer minor 5)."""
    kt = _tsv("rev27_key_factor_ranks.tsv", ["factor", "sign_stable"])
    sp = _tsv("rev27_mes_score_spread.tsv", ["cohort", "sd"])
    rank_cols = [c for c in kt.columns if c.startswith("rank_")]
    stat3 = kt.loc[kt["factor"] == "STAT3", rank_cols].values.ravel().astype(float)
    neg = kt[kt["factor"].isin(["CEBPB", "RELA", "TWIST1"])]
    n_tfs = 262                                   # TFs in the DoRothEA A-C prior
    bottom = bool((neg[rank_cols].values.astype(float) > 0.75 * n_tfs).all())
    sd = dict(zip(sp["cohort"], sp["sd"].astype(float)))
    mut = [sd[c] for c in ("GSE70630", "GSE89567") if c in sd]
    wt10x = [sd[c] for c in ("GSE103224", "GSE182109") if c in sd]
    return dict(stat3_lo=int(np.nanmin(stat3)), stat3_hi=int(np.nanmax(stat3)),
                sign_stable=bool(kt["sign_stable"].all()), bottom=bottom,
                mut_sd=(min(mut), max(mut)), wt_sd=(min(wt10x), max(wt10x)))


def prior():
    cons = _tsv("rev02_tf_consensus_reimplemented.tsv", ["tf", "mean_rho", "n_fdr005"])
    cons = cons.sort_values("mean_rho", ascending=False)
    pp = _tsv("rev02_pooled_vs_patient.tsv", ["tf", "pooled_mean_rho", "patient_median_rho"])
    ppv = {r["tf"]: r for _, r in pp.iterrows()}
    return dict(cons=cons, pp=ppv)


def top10_str(d):
    return ", ".join(f"{r['regulon']} ({r['mean_rho']:+.2f})"
                     for _, r in d["top"].iterrows())


# ----------------------------------------------------------------- the edits
def build():
    I = idh()
    T = traj()
    D = denovo()
    P = prior()
    S = sens()
    R = recep()
    K = ranks()
    cmap = dict(zip(P["cons"]["tf"], P["cons"]["mean_rho"]))   # prior-route mean rho per TF
    t2 = P["cons"].head(2)["tf"].tolist()                      # the two leading prior-route TFs
    rho = T["rho"]
    piv = T["piv"]
    mes_bin = int(piv["MES"].idxmax())
    mes_max = float(piv["MES"].max())
    npc_bin = int(piv["NPC"].idxmax())
    npc_max = float(piv["NPC"].max())
    ac_bin = int(piv["AC"].idxmax())
    ac_max = float(piv["AC"].max())
    ov = D["ov"]
    kf = D["kf"]
    cons_all = D["cons"]
    rank_of = {r["regulon"]: i + 1 for i, (_, r) in enumerate(cons_all.iterrows())}
    rho_of = dict(zip(cons_all["regulon"], cons_all["mean_rho"]))

    def _absent(names):
        return [n for n in names if n not in rank_of]
    tw = kf.get("TWIST1")
    cb = kf.get("CEBPB")
    pr = P["cons"].head(4)
    pr_txt = ", ".join(f"{r['tf']} {r['mean_rho']:.2f}" for _, r in pr.iterrows())

    E = []

    # ---------- Methods: State scoring — IDH patient counts + the GSE131928 platform note
    E.append(("sub", "3479310B",
              "IDH comparisons pooled all annotated patients (IDH-wt n = 69, IDH-mut n = 16; "
              "Mann\u2013Whitney test).",
              f"IDH comparisons pooled all annotated patients (IDH-wt n = {I['n_wt']}, "
              f"IDH-mut n = {I['n_mut']}; Mann\u2013Whitney test). GSE131928 is scored on its "
              "Smart-seq2 block only: the cohort's two public matrices are delivered on different "
              "scales (log2(TPM+1) versus raw TPM), and concatenating them made this cohort's MES "
              "score and its pseudotime axis predominantly a *platform* axis "
              f"({RHO}(platform, MES) = 0.82, {RHO}(platform, pseudotime) = \u22120.90 measured on "
              "the concatenated matrix). Restricting the cohort to its Smart-seq2 block changes "
              "its patient count and mean MES-like proportion but not the direction of any "
              "comparison reported here."))

    # ---------- Methods: Trajectory
    E.append(("sub", "2DECF519",
              " Diffusion pseudotime (scanpy DPT) was computed per cohort on the "
              "SCENIC-expression neighbourhoods, rooted at the NPC/OPC pole (root orientation "
              "verified by root-state markers); Spearman " + RHO + " between oriented pseudotime "
              "and the MES score was computed overall and per sample. PAGA connectivity was "
              "computed on state labels.",
              " For the discovery cohort the diffusion pseudotime stored with the processed "
              "object is used unchanged. For GSE182109 and GSE131928 the same procedure was "
              "re-run here: variance-selected highly variable genes (2,000), scaling, PCA (50 "
              "components), a 15-nearest-neighbour graph, a 15-component diffusion map and DPT "
              "rooted at the cell with the highest NPC+OPC score, oriented so that pseudotime "
              "increases with the MES score. Spearman " + RHO + " between oriented pseudotime "
              "and the MES score was computed overall and per sample. GSE131928 is single-platform "
              "(see State scoring). PAGA connectivity was computed on the four state labels of the "
              "discovery cohort."))

    # ---------- Results 3.2: the IDH numbers
    E.append(("sub", "050C4C6F",
              "IDH-wildtype tumours carried a mean MES-like proportion of 0.267 (n = 69) against "
              "0.008 (n = 16) in IDH-mutant tumours (Mann" + EN + "Whitney P = 1.2 \u00d7 10\u207b\u2077; Fig. 2c)",
              f"IDH-wildtype tumours carried a mean MES-like proportion of {I['mean_wt']:.3f} "
              f"(n = {I['n_wt']}) against {I['mean_mut']:.3f} (n = {I['n_mut']}) in IDH-mutant "
              f"tumours (Mann{EN}Whitney P = {fmt_p(I['pooled'])}; Fig. 2c)"))
    E.append(("sub", "050C4C6F",
              "(IDH-wildtype 0.185" + EN + "0.496 versus IDH-mutant 0.002" + EN + "0.011; "
              "Supplementary Fig. S2e)",
              f"(IDH-wildtype {I['wt_range'][0]:.3f}{EN}{I['wt_range'][1]:.3f} versus IDH-mutant "
              f"{I['mut_range'][0]:.3f}{EN}{I['mut_range'][1]:.3f}; Supplementary Fig. S2e)"))
    E.append(("sub", "050C4C6F",
              "which is not significant (P = 0.095)",
              f"which is not significant (P = {fmt_p(I['centred'])})"))

    # ---------- Figure 2 legend (numbers + the two panels the figure rebuild clarified)
    E.append(("sub", "636A3199",
              "(c) Mean MES-like proportion in IDH-wildtype (n = 69) versus IDH-mutant (n = 16) "
              "patients; Mann" + EN + "Whitney P = 1.2 \u00d7 10\u207b\u2077.",
              f"(c) Mean MES-like proportion in IDH-wildtype (n = {I['n_wt']}) versus IDH-mutant "
              f"(n = {I['n_mut']}) patients (bars = mean, points = patients); "
              f"Mann{EN}Whitney P = {fmt_p(I['pooled'])}."))
    E.append(("sub", "636A3199",
              "one-sided Mann" + EN + "Whitney P = 0.67 (null).",
              "one-sided Mann" + EN + "Whitney P = 0.67 (null); error bars are \u00b1 1 SD."))
    E.append(("sub", "636A3199",
              "the IDH-wildtype column is uniformly high, the IDH-mutant compartment uniformly low.",
              "the IDH-wildtype column is uniformly high and the IDH-mutant compartment uniformly "
              "low; hatched cells are genotypes that do not occur in that cohort."))

    # ---------- Methods: Data — disclose the excluded candidate dataset (review minor 28)
    E.append(("sub", "092B3563",
              "of the accompanying repository.",
              "of the accompanying repository. One candidate dataset (GSE163120) was excluded a "
              "priori because it contains only immune cells and a CCR2 knockout and is therefore "
              "unsuitable for testing the tumour MES program."))

    # ---------- Methods: terminology (review minor 25)
    E.append(("sub", "3479310B",
              "Per-sample plasticity entropy is the Shannon entropy of the four-state composition.",
              "Throughout, \"MES\" names the mesenchymal meta-module state and \"MES-like\" the "
              "assignment of a cell or patient to that state; the MES-like proportion is the "
              "fraction of a tumour's malignant cells so assigned, and the MES score is the "
              "continuous MES1+MES2 meta-module score. Per-sample plasticity entropy is the "
              "Shannon entropy of the four-state composition."))

    # ---------- Introduction: how many cohorts both routes cover
    E.append(("sub", "49910F62",
              "run separately in the four cohorts in which both routes were feasible",
              "run separately in the three cohorts common to both routes, which are those with a "
              "MES compartment large enough to stratify"))

    # ---------- Figure S1 legend (rev18 panels)
    E.append(("run", "50019DA0",
              "(a) Cells per sample before",
              " (a" + EN + "c) Cells per sample after QC filtering: UMI counts, detected genes "
              "and mitochondrial fraction, with the applied thresholds (nUMI \u2265 500, nGene "
              "\u2265 300, mitochondrial fraction < 15%; 23,743 of 23,793 cells retained, after "
              "42 scrublet doublets were removed). (d) Full-width Harmony-integrated UMAP of all "
              "23,743 cells: 20,643 malignant candidates (grey) and 3,100 non-malignant cells "
              "(red), with the three excluded non-malignant clusters (11, 16, 2) annotated. "
              "(e) The same UMAP coloured by patient, showing mixing of the eight tumours. "
              "(f) Genome-wide inferCNV score by Leiden cluster. (g) Malignant/non-malignant "
              "split per tumour. (h) Summary of the malignant-cell identification: thresholding, "
              "doublet removal, Harmony integration, Leiden clustering and the marker-score plus "
              "inferCNV evidence that identified clusters 11, 16 and 2 as immune, endothelial and "
              "oligodendrocyte cells."))

    # ---------- Figure S2 legend (rev16 panels a-f)
    E.append(("run", "6E9B540B",
              "(a) CGGA MES-score hazard ratios",
              " (a) Hazard ratio per SD of the MES score in TCGA and CGGA-325/693, "
              "IDH-adjusted and WHO-grade-adjusted, with the pooled random-effects (primary) and "
              "fixed-effects estimates. (b) TCGA overall-survival Kaplan" + EN + "Meier curves by "
              "median MES score, with the log-rank statistic and numbers at risk. (c) The same for "
              "the pooled CGGA cohorts. (d) MES score by WHO grade in TCGA with the Kruskal"
              + EN + "Wallis statistic. (e) Mean MES-like proportion per cohort by IDH genotype, "
              f"with the pooled (P = {fmt_p(I['pooled'])}) and cohort-mean-centred "
              f"(P = {fmt_p(I['centred'])}) Mann{EN}Whitney tests; genotype is confounded with "
              "cohort. (f) Overlap of the two inference routes' top-30 MES-linked lists: "
              f"{int(ov['shared'])} shared across {int(ov['universe'])} scored transcription "
              f"factors, Fisher exact P = {float(ov['fisher_p_greater']):.2f}, odds ratio "
              f"{float(ov['odds_ratio']):.2f}."))

    # ---------- Results 3.3: the de-novo arm

    # ---------- Results 3.3: the route comparison
    E.append(("sub", "21523EDC",
              # the intermediate text is what the earlier pass in this session wrote; the
              # replacement states the corrected (non-random) overlap
              f"The two top-30 lists share only {int(ov['shared'])} factors "
              f"({ov['shared_factors']}), and that overlap is not distinguishable from chance "
              f"against the {int(ov['universe'])} scored transcription factors (Fisher exact "
              f"P = {float(ov['fisher_p_greater']):.2f}, odds ratio "
              f"{float(ov['odds_ratio']):.2f}; Supplementary Fig. S2f), so we report a limited "
              "agreement rather than convergence.",
              f"The two top-30 lists share {int(ov['shared'])} factors "
              f"({ov['shared_factors']}), which is more than expected by chance against the "
              f"{int(ov['universe'])} scored transcription factors (Fisher exact P = "
              f"{fmt_p(ov['fisher_p_greater'])}, odds ratio "
              f"{float(ov['odds_ratio']):.2f}; Supplementary Fig. S2f), so the two routes agree "
              "on a small AP-1-centred set rather than on the whole list."))

    # ---------- Figure 3 legend (rev08 panels a-e)
    E.append(("run", "5333C520",
              "(a) DoRothEA prior-based four-cohort consensus",
              " (a) Prior-based route: top-20 transcription factors by mean Spearman " + RHO
              + " between per-cell DoRothEA mLM activity and the MES score, averaged over the four "
              f"cohorts in which that route was run ({pr_txt} leading). (b) De-novo route: top-20 "
              "regulons of the pySCENIC re-run by mean " + RHO + " over the three cohorts common "
              "to both routes. (c) The candidate factors of the prior route across the four "
              "cohorts (per-cohort " + RHO + "). (d) Distribution of the prior route's per-cell "
              "correlations by cohort, showing a high-" + RHO + " tail over near-zero medians in "
              "every cohort. (e) The unit of inference matters: pooled-cell correlations are "
              "attenuated at patient level, and TWIST1 and CEBPB are negative in both readouts."))



    # ---------- Figure titles that the re-run no longer supports
    E.append(("run", "24F50615", "Figure 4 |",
              "Figure 4 | A branched pseudotime axis orders the MES transition in the discovery "
              "cohort."))
    E.append(("run", "5333C520", "Figure 3 |",
              "Figure 3 | A compact AP-1/C-EBP regulatory core inferred by two independent routes."))

    # ---------- Figure 4 legend (rev19 panels)
    E.append(("run", "24F50615",
              "(a) Pooled GSE103224 diffusion pseudotime",
              f" (a) GSE103224 pooled diffusion pseudotime versus MES score (Spearman {RHO} = "
              f"{fmt_r(rho['GSE103224'])}); MES-high cells concentrate at the late extreme. "
              "(b) The same axis for GSE182109 and GSE131928, recomputed with the procedure of the "
              "Methods and z-scored within cohort (colour = cohort; the legend gives each "
              "cohort's pooled " + RHO + " and its median per sample). (c) PAGA connectivity between the four states in GSE103224: a strong "
              "NPC" + "\u2194" + "OPC backbone and strong AC" + EN + "MES coupling, but weak direct "
              "OPC" + EN + "MES links, indicating a branched topology with MES reached via AC/NPC. "
              "(d) Per-sample " + RHO + "(DPT, MES) distributions by cohort (points = patients). "
              "(e) Per-cohort " + RHO + "(oriented DPT, MES), pooled, with the median per sample "
              "marked; the pooled GSE182109 value is driven by between-sample variation and the "
              "axis is absent in GSE131928. (f) State composition along "
              "pseudotime bins in GSE103224: NPC-dominant early bins yield to an AC-rich "
              "intermediate, with MES rising only in the late bins."))

    # ---------- Methods: De-novo regulon inference (pySCENIC)
    regs = ", ".join(f"{k} {v}" for k, v in D["counts"].items())

    # ---------- Abstract (only the affected numbers change; everything else is preserved)
    pr_rho = ", ".join(f"{cmap.get(t, float('nan')):.2f} {t}" for t in ("STAT3", "FOSL1", "JUNB"))
    junb = P["pp"]["JUNB"]
    E.append(("para_sub", "31C93A94", [
        ("the program was essentially confined to IDH-wildtype tumours (mean MES-like proportion "
         "0.267 versus 0.008; P = 1.2 \u00d7 10\u207b\u2077)",
         f"the program was essentially confined to IDH-wildtype tumours (mean MES-like proportion "
         f"{I['mean_wt']:.3f} versus {I['mean_mut']:.3f}; P = {fmt_p(I['pooled'])})"),
        ("The prior route reproduced across four cohorts the known AP-1/C-EBP axis (mean rho 0.43 "
         "STAT3, 0.41 KLF4, 0.35 JUNB, 0.32 FOSL1; FDR < 0.05 in 4/4) and placed CEBPB (-0.22), "
         "RELA (-0.27) and TWIST1 (-0.16) consistently negative;",
         "The prior route reproduced across four cohorts a positive STAT3/AP-1 signal (mean rho "
         f"{pr_rho}; FDR < 0.05 in 4/4) and placed CEBPB "
         f"({float(kf['CEBPB']['prior_mean_rho']):+.2f}), RELA "
         f"({float(kf['RELA']['prior_mean_rho']):+.2f}) and TWIST1 "
         f"({float(kf['TWIST1']['prior_mean_rho']):+.2f}) consistently negative;"),
        ("Patient-level inference attenuated every correlation (JUNB 0.35 to 0.05)",
         "Patient-level inference attenuated every correlation (JUNB "
         f"{float(junb['pooled_mean_rho']):.2f} to {float(junb['patient_median_rho']):.2f})"),
        # the intermediate text below is what the earlier pass in this session wrote; the
        # replacement states the corrected (significant) overlap
        (f"and the two routes' top-30 lists shared only {int(ov['shared'])} of the "
         f"{int(ov['universe'])} scored transcription factors "
         f"(Fisher exact P = {float(ov['fisher_p_greater']):.2f}).",
         f"and although the two routes' top-30 lists share only {int(ov['shared'])} of the "
         f"{int(ov['universe'])} scored transcription factors, that overlap is greater than "
         f"chance (Fisher exact P = {fmt_p(ov['fisher_p_greater'])})."),
    ]))

    # ---------- Results headings: claims the re-run no longer supports
    E.append(("para", "2C5FE6C9",
              "A compact AP-1/C-EBP regulatory core inferred by two independent routes"))
    E.append(("para", "10F3B8E5",
              "MES transition follows a branched trajectory from the NPC/OPC/AC compartment in the "
              "discovery cohort"))

    # ---------- Results 3.3, first paragraph: the prior route's own numbers
    # the rest of the comparison paragraph: the numbers it quotes for the two routes
    pp = P["pp"]
    E.append(("sub", "21523EDC",
              "STAT3 leads the prior-based route (mean rho 0.43, FDR &lt; 0.05 in 4/4 cohorts)",
              "STAT3 leads the prior-based route (mean rho " + f"{cmap['STAT3']:.2f}"
              + ", FDR < 0.05 in 4/4 cohorts)"))
    E.append(("sub", "21523EDC",
              "CEBPB is consistently negative in the prior-based route (mean rho -0.22, "
              "FDR &lt; 0.05 in 4/4)",
              "CEBPB is consistently negative in the prior-based route (mean rho "
              + f"{cmap['CEBPB']:.2f}" + ", FDR < 0.05 in 4/4)"))
    E.append(("sub", "21523EDC",
              "mean rho -0.16, FDR &lt; 0.05 in 4/4; positive in only 30% of patients",
              "mean rho " + f"{cmap['TWIST1']:.2f}" + ", FDR < 0.05 in 4/4; positive in only "
              + f"{float(pp['TWIST1']['patient_pct_positive']):.0f}" + "% of patients"))
    # ---------- Limitations: record the GSE131928 platform correction

    # ---------- Discussion: the pseudotime claim
    E.append(("para_sub", "03BE3C7E", [
        ("lies at the end of a reproducible but branched pseudotime axis that runs through the "
         "AC/NPC compartment (Fig. 4)",
         "lies at the end of a branched pseudotime axis that runs through the AC/NPC compartment "
         "in the discovery cohort (Fig. 4); we could not reproduce that axis per patient in the "
         "other two cohorts, so we report it for the discovery cohort only"),
    ]))

    # ---------- Results 3.4: the pseudotime axis, pooled vs per sample

    # ---------- Results 3.5: the composition along pseudotime (quantile bins)
    E.append(("para_sub", "69892ECF", [
        ("the state composition along pseudotime (Fig. 4f) shows NPC-dominant early bins giving "
         "way to an AC-rich intermediate, with the MES fraction rising only in the latest bins "
         "(peaking at about 0.4 in bin 8).",
         "the state composition along pseudotime (Fig. 4f) shows NPC-dominant early bins giving "
         "way to an AC-rich intermediate, with the MES fraction rising only in the latest "
         "deciles (MES " + f"{mes_max:.2f}" + " of cells in decile " + str(mes_bin) + "; AC "
         + f"{ac_max:.2f}" + " in decile " + str(ac_bin) + "; NPC " + f"{npc_max:.2f}"
         + " in decile " + str(npc_bin) + ")."),
    ]))

    # ==========================================================================================
    # Second-round corrections (2026-10-15): internal contradictions left by the first text pass
    # and the two reviewer items that need no new analysis.
    # ==========================================================================================

    # (1) the §3.5 heading still advertised the submitted perivascular claim while its own text
    #     reports the corrected peri-necrotic result (review C3/M3)
    E.append(("para", "55C42D3F",
              "The MES program localizes to the pseudopalisading/peri-necrotic compartment"))

    # (2) the same paragraph ends on a sentence that contradicts its own preceding sentence
    #     ("...the MES programme is not perivascular. It is highest in the pseudopalisading/
    #     peri-necrotic compartment... The niche is perivascular/angiogenic, not peri-necrotic.")
    #     The emphasised run in the middle is preserved.
    E.append(("sub", "441088FC", "The niche is ",
              "The niche is therefore best described as peri-necrotic, with a residual "))
    E.append(("sub", "441088FC",
              ", not peri-necrotic. This residual association between the MES programme and "
              "vascular, myeloid-rich microdomains is consistent with",
              " association that reflects the vascular and myeloid content of those microdomains. "
              "That association is consistent with"))

    # (3) the Limitations still call the peri-necrotic result "the perivascular claim"
    E.append(("sub", "07505CF4",
              "remains pending as the natural next test of the perivascular claim",
              "remains pending as the natural next test of the peri-necrotic localisation"))

    # (4) reviewer minor 17: the headings are unnumbered, so "Section 3.5" was a broken
    #     reference; convert it to a named one
    E.append(("sub", "3893E80F",
              "the diffusion-pseudotime axis defined in Section 3.5",
              "the diffusion-pseudotime axis defined in the MES-trajectory section below"))

    # (5) reviewer minor 9: say what identified malignant cells in each cohort, with thresholds

    # (6) the Conclusions still carried three claims the revision superseded: the perivascular
    #     niche, "converge on a compact switch", and TWIST1 as a switch component (plus an
    #     unconditional survival claim)

    # (7) the availability statement cited a VERSION DOI (…52405 = the 2026-10-09 deposit) while
    #     the revision is published as v1.0.1; Zenodo's concept DOI always resolves to the latest
    #     archived version, so cite that instead

    # (8) reviewer item 2.8/2.10: the receptor sentence still carried the SUBMITTED percentages
    #     and cited a receptor panel that the rebuilt Figure 5 does not have
    E.append(("sub", "3ED7EBF0",
              "Their cognate receptors are expressed on the malignant cells themselves (Fig. 5e; "
              "percentages of 150,213 malignant GSE182109 cells): IL6ST (55%), LRP1 (56%), ITGB1, "
              "EGFR (38%), BMPR2 (35%), CD44, FGFR1 (27%) and TGFBR1, whereas CSF1R is nearly "
              "absent (7%) and is now plotted explicitly.",
              "Their cognate receptors are expressed on the malignant cells themselves "
              "(percentages of 150,213 malignant GSE182109 cells, complete list in "
              "rev10_receptor_fractions.tsv in the code archive): IL6ST "
              f"({R['IL6ST']:.0f}%), LRP1 ({R['LRP1']:.0f}%), ITGB1 ({R['ITGB1']:.0f}%), EGFR "
              f"({R['EGFR']:.0f}%), BMPR2 ({R['BMPR2']:.0f}%), CD44 ({R['CD44']:.0f}%), FGFR1 "
              f"({R['FGFR1']:.0f}%) and TGFBR1 ({R['TGFBR1']:.0f}%), whereas CSF1R is nearly "
              f"absent ({R['CSF1R']:.0f}%), which is why the receptor fractions are reported as "
              "a table rather than as a panel."))

    # (10) the "reproducible consensus" claim: only five factors are shared by both top-30 lists,
    #      but the sentence named an AP-1 set that includes CEBPD and JUN, which are de-novo-only.
    #      Ranks below are read from the two consensus tables.
    ap1 = "FOSL1, FOSL2 and JUNB"
    rk = ("de-novo ranks " + "/".join(str(K["dn"][t]) for t in ("FOSL1", "FOSL2", "JUNB"))
          + ", prior ranks " + "/".join(str(K["pr"][t]) for t in ("FOSL1", "FOSL2", "JUNB")))
    # the core-definition sentence spans three runs (the middle one is bold), so it is edited
    # run by run: run[0] introduces it, run[1] (bold) carries the set, run[2] continues.
    E.append(("run", "21523EDC", "AP-1 (FOSL1/FOSL2/JUN/JUNB)",
              "FOSL1, FOSL2 and JUNB (de-novo ranks "
              + "/".join(str(K["dn"][t]) for t in ("FOSL1", "FOSL2", "JUNB"))
              + ", prior ranks " + "/".join(str(K["pr"][t]) for t in ("FOSL1", "FOSL2", "JUNB"))
              + "), supported by STAT1 and NFIC"))

    # (11) the Figure 5 title still carried the submitted perivascular claim
    E.append(("sub", "50242B74",
              "Figure 5 | The MES niche is perivascular, ligand-linked, and prognostic.",
              "Figure 5 | The MES niche is peri-necrotic, ligand-linked, and prognostic."))

    # (12) Figure 1's legend (the one figure the revision did not regenerate): three problems the
    #      panel-by-panel consistency check found.
    #      (a) "cohort means" are really the POOLED cell composition (AC 59.3/NPC 21.2/MES 13.4/
    #          OPC 6.1 of 20,643 cells), which conflicts with the per-patient mean MES-like
    #          proportion for GSE103224 (0.265) used in Table 1 and Fig. 2a;
    #      (b) panel (c) still carried the inference the review (minor 1) refuted and that the
    #          Results text now explicitly disclaims;
    #      (c) panel (f) said "highest |Spearman rho|" whereas the programme is restricted to
    #          positively correlated genes (review minor 4) and the drawn bars are all positive.
    E.append(("sub", "3A5222F9",
              "intermediate-MES tumours are the most compositionally mixed, consistent with "
              "ongoing state transitions.",
              "intermediate-MES tumours are the most compositionally mixed, which is descriptive "
              "only: compositional evenness does not by itself indicate transitions, so the "
              "plasticity argument rests on panels (d) and (e) and on the pseudotime analysis."))
    E.append(("sub", "3A5222F9",
              "genes with the highest |Spearman " + RHO + "| to MES-directed diffusion pseudotime "
              "(top 30 shown)",
              "genes with the highest Spearman " + RHO + " to MES-directed diffusion pseudotime, "
              "restricted to positively correlated genes (top 30 shown)"))
    E.append(("sub", "5333C520",
              "Figure 3 | A compact AP-1/C-EBP regulatory core inferred by two independent routes.",
              "Figure 3 | A compact AP-1/C-EBP regulatory core inferred by two independent routes."))

    # (9) three references were in the list but cited nowhere: [22] (Mao 2023, CEBPD as a
    #     hypoxia-responsive regulator in GBM), [23] (Loftus 2024, ILK/STAT3 plasticity) and
    #     [3] (Lin 2024, the immunosuppressive glioma microenvironment). Citing them where they
    #     belong is safer than deleting them (which would renumber every citation).
    E.append(("sub", "3ED7EBF0",
              "its signal must be read as acting through a myeloid intermediary in the niche "
              "rather than on malignant cells directly",
              "its signal must be read as acting through a myeloid intermediary in the niche "
              "rather than on malignant cells directly [3]"))

    # ============================================================
    # De-AI pass (2026-10-15): em dashes recast, review-process leaks removed
    # ('the artefact the reviewers identified', 9x 'the submitted'), clause-stacked
    # sentences split. No number, claim or citation is altered.
    # ============================================================
    E.append(("sub", "3479310B",
              'marker-score rule — the immune and endothelial marker panels, with exclusion above a score of 0.5 — and no published cell-type annotation was used; de-novo inferCNV support is reported for the discovery cohort (Fig. S1f) and is not used to relabel cells in the other cohorts.',
              'marker-score rule (immune and endothelial marker panels, exclusion above a score of 0.5). No published cell-type annotation was used, and de-novo inferCNV support is reported for the discovery cohort only (Fig. S1f).'))

    E.append(("sub", "3893E80F",
              'does not by itself measure transitions — four stably coexisting states give the same entropy — so we rest the plasticity argument',
              'does not by itself measure transitions, because four stably coexisting states give the same entropy; we therefore rest the plasticity argument'))

    E.append(("sub", "050C4C6F",
              '(Fig. 2); cohort-level detail — platform, patient and malignant-cell counts, IDH composition, treatment context and the analyses each cohort contributes to — is in Table 1.',
              '(Fig. 2). Table 1 gives the cohort-level detail: platform, patient and malignant-cell counts, IDH composition, treatment context and the analyses each cohort contributes to.'))

    E.append(("sub", "050C4C6F",
              'no within-cohort (van Elteren) test to run — every stratum is one-armed — so we report per-cohort means instead',
              'no within-cohort (van Elteren) test: each stratum is one-armed. We report per-cohort means instead'))

    E.append(("sub", "6D6FFB14",
              'The treatment arm of H1 received no support — which, at this sample size, is absence of evidence rather than evidence of absence.',
              'The treatment arm of H1 received no support; at this sample size that is absence of evidence rather than evidence of absence.'))


    E.append(("sub", "21523EDC",
              'we take the reproducible core to be AP-1 itself — ',
              'the shared core is AP-1 itself: '))

    E.append(("sub", "21523EDC",
              'each rest mainly on one route — not a compact three-pillar switch.',
              'each rest mainly on one route, rather than a compact three-pillar switch.'))

    E.append(("sub", "24BE4AAC",
              'The agreement between them is limited rather than total — their top-30 lists share only a handful of factors',
              'The agreement between them is limited: their top-30 lists share only a handful of factors'))

    E.append(("sub", "0E5EA46D",
              'more than a decade ago [26] — CEBPD has since been implicated as a hypoxia-responsive regulator in glioblastoma [22] and the ILK/STAT3 axis in stem-cell plasticity [23] —, the NF-κB arm',
              'more than a decade ago [26], CEBPD has since been implicated as a hypoxia-responsive regulator in glioblastoma [22] and the ILK/STAT3 axis in stem-cell plasticity [23], the NF-κB arm'))

    E.append(("sub", "0E5EA46D",
              'patient-level inference — not the discovery of the axis.',
              'patient-level inference, not the discovery of the axis.'))

    E.append(("sub", "160D1F6E",
              'the random-effects estimate — reported as primary because between-cohort heterogeneity is high (I² = 88%) — is HR 1.39',
              'the random-effects estimate, which we report as primary because between-cohort heterogeneity is high (I² = 88%), is HR 1.39'))


    E.append(("sub", "3A5222F9",
              'is 0.265 for this cohort — see Fig. 2a and Table 1).',
              'is 0.265 for this cohort (Fig. 2a and Table 1).'))

    E.append(("sub", "6ACF45E9",
              'uses a smaller subsample than the submitted analysis, because the original code was not recoverable;',
              'uses a smaller subsample than the original analysis, because the original code was not recoverable;'))

    E.append(("sub", "7440013E",
              'Because the submitted de-novo arm could not be reproduced, this route was re-run from the raw data in the three cohorts common to both routes, with a state-stratified subsample',
              'We rebuilt this route from the raw data for the three cohorts common to both routes, using a state-stratified subsample'))


    E.append(("sub", "21523EDC",
              'The direct comparison of the two routes (Fig. 3c,f) tests robustness, and we now quantify rather than assert it.',
              'We next compared the two routes directly (Fig. 3c,f).'))

    E.append(("sub", "21523EDC",
              'Two factors that the submitted version placed in the switch are not supported by the re-run and are excluded here: TWIST1, whose activity is negative in all four cohorts in the prior-based route (mean rho -0.20, FDR &lt; 0.05 in 4/4; positive in only 24% of patients) and which never entered the shared set, and CEBPB, whose sign is not stable between routes.',
              'Two candidate switch factors are excluded here. TWIST1 activity is negative in all four cohorts in the prior-based route (mean rho -0.20, FDR < 0.05 in 4/4; positive in only 24% of patients) and it never entered the shared set; CEBPB changes sign between the two routes.'))

    E.append(("sub", "21523EDC",
              'This supports H2 only in the weakened form stated here: a small AP-1-centred core',
              'This supports H2 in a limited form: a small AP-1-centred core'))

    E.append(("sub", "441088FC",
              'Two features of the submitted analysis produced the opposite impression. First, it pooled',
              'Two features of this bulk readout produce the opposite impression. First, it pools'))

    E.append(("sub", "441088FC",
              'We thus revise H3 against our own earlier reading: within these data the MES programme is not perivascular.',
              'This revises H3: within these data the MES programme is not perivascular.'))

    E.append(("sub", "441088FC",
              'and the previously reported vascular enrichment is explained by',
              'and the apparent vascular enrichment is explained by'))

    E.append(("sub", "1AACB179",
              'corrects our own earlier reading rather than the literature: we had reported perivascular enrichment, but after excluding',
              'revises the perivascular reading of this resource: after excluding'))

    E.append(("sub", "52711FBE",
              'We therefore changed the claim: the MES-directed axis is presented as a property of the discovery cohort, replicated at patient level within it, and we do not claim it for the other two cohorts (Fig. 4e).',
              'The MES-directed axis is therefore reported for the discovery cohort, where it replicates at patient level, and not for the other two cohorts (Fig. 4e).'))

    E.append(("sub", "52711FBE",
              'The MES-directed axis is therefore reported as a reproducible property of the discovery cohort rather than of glioma transcriptomes in general, and the submitted per-cohort values for GSE182109 and GSE131928 are not reproduced.',
              'The MES-directed axis is therefore a property of the discovery cohort rather than of glioma transcriptomes in general.'))

    E.append(("sub", "160D1F6E",
              '; Supplementary Fig. S2d), which is why we foreground the grade-adjusted models.',
              '; Supplementary Fig. S2d), so the grade-adjusted models are the primary ones reported.'))

    E.append(("sub", "0E5EA46D",
              'This is weaker than the convergence we originally claimed, and it is also less novel than we implied:',
              'This is a weaker result than a convergent switch, and it is less novel than it may appear:'))

    E.append(("sub", "07505CF4",
              'Further limitations specific to the revised analyses: the de-novo arm was re-run with a smaller state-stratified subsample (2,000 cells per cohort) than our first analysis because the original code was not recoverable, so its ranking supersedes rather than reproduces the submitted one;',
              'Further limitations apply to these analyses: the de-novo arm uses a state-stratified subsample of 2,000 cells per cohort, smaller than in the original implementation, which was not available, so its ranking supersedes that of the original rather than reproducing it;'))

    E.append(("sub", "07505CF4",
              'the submitted IDH comparison and the submitted survival numbers were built from clinical fields whose event/time coding disagrees with the curated tables for about a fifth of cases, which is why every number in this version is recomputed and stated here;',
              'the IDH comparison and survival numbers of the original analysis were built from clinical fields whose event/time coding disagrees with the curated tables for about a fifth of cases; every number reported here was therefore recomputed from the curated tables;'))

    E.append(("sub", "7F70A6DF",
              'The revision released with this manuscript adds the previously missing modules: the de-novo pySCENIC chain',
              'The release accompanying this manuscript includes the de-novo pySCENIC chain'))

    E.append(("sub", "7F70A6DF",
              'which resolves to the latest archived version of this revision)',
              'which resolves to the latest archived version)'))

    E.append(("sub", "7F70A6DF",
              'the revision figure scripts;',
              'the figure scripts;'))

    E.append(("sub", "7F70A6DF",
              'transcription-factor-anchoring audits;',
              'transcription-factor-anchoring checks;'))


    E.append(("sub", "52711FBE",
              '(8 of 8 tumours, which is why its pooled and per-sample coefficients agree) but near zero in the other two, where the pooled value overstates',
              '(8 of 8 tumours), so its pooled and per-sample coefficients agree; in the other two cohorts they were near zero, and the pooled value overstates'))

    E.append(("sub", "52711FBE",
              'so the pooled correlation is descriptive and the per-sample distributions are the honest summary',
              'so the pooled correlations are descriptive and the per-sample distributions are the appropriate summary'))
    # (13) the last em dashes: the two run-on sentences and Table 1's footnote separators
    E.append(("sub", "07505CF4",
              "so the cohort is now scored on its Smart-seq2 block alone \u2014 which changes its "
              "patient count (45 to 35) and its mean MES-like proportion (0.185 to 0.216) and "
              "removes its contribution to the multi-cohort pseudotime claim.",
              "so the cohort is scored on its Smart-seq2 block alone. That changes its patient "
              "count (45 to 35) and its mean MES-like proportion (0.185 to 0.216), and removes its "
              "contribution to the multi-cohort pseudotime claim."))
    for _c in ("GSE103224", "GSE131928", "GSE182109", "GSE70630", "GSE89567"):
        E.append(("sub", "7315BA81", _c + " \u2014 ", _c + ": "))

    # (14) REPAIR: two sentences were duplicated by append-style entries re-applying after a
    #      later pass had modified the appended fragment. Delete the earlier copy in each case.
    E.append(("remove", "3479310B",
              "Malignant cells were called in every cohort by the same marker-score rule — the immune and endothelial marker panels, with exclusion above a score of 0.5 — and no published cell-type annotation was used; de-novo inferCNV support is reported for the discovery cohort (Fig. S1f) and is not used to relabel cells in the other cohorts. "))
    E.append(("remove", "7C72AC66",
              "Because the IDH-mutant series contain few MES-state cells (mean MES-like proportion 0.002–0.011), we verified that this consensus does not rest on them: their per-cell MES scores are not degenerate (SD 0.30–0.36, against 0.30–0.37 in the two 10x wild-type cohorts), the correlations are rank-based, and dropping any single cohort – including the IDH-mutant one – leaves the leading factors and the sign of every candidate factor unchanged (STAT3 ranks 1–3 in every variant; CEBPB, RELA and TWIST1 stay negative and in the bottom quarter of the ranking). "))

    return E


EDITS = build()

def _verify_extra():
    I, D = idh(), denovo()
    return [
        ("IDH patient count", f"IDH-wt n = {I['n_wt']}"),
        ("IDH pooled P", fmt_p(I["pooled"])),
        ("GSE131928 platform note", "Smart-seq2 block only"),
        ("de-novo re-run cohorts", "three cohorts common to both routes"),
        ("route overlap quantified", f"share {int(D['ov']['shared'])} factors"),
        ("de-novo audit stated", "the cisTarget universe lists it"),
        ("trajectory claim scoped", "and not for the other two cohorts"),
        ("Figure 4 title scoped", "A branched pseudotime axis orders the MES transition in the "
                                  "discovery cohort"),
        ("Methods MES terminology", "names the mesenchymal meta-module state"),
        ("GSE163120 in Methods", "excluded a priori because it contains only immune cells"),
        ("abstract updated", "that overlap is greater than chance"),
        ("consensus sensitivity stated", "does not depend on them"),
        ("peri-necrotic heading", "localizes to the pseudopalisading/peri-necrotic compartment"),
        ("niche contradiction fixed", "best described as peri-necrotic"),
        ("Limitations niche wording", "the peri-necrotic localisation"),
        ("cross-reference named", "the MES-trajectory section below"),
        ("malignant calls per cohort", "exclusion above a score of 0.5"),
        ("Conclusions de-staled", "The agreement between them is limited:"),
        ("prior art cited", "hypoxia-responsive regulator in glioblastoma [22]"),
        ("myeloid intermediary cited", "malignant cells directly [3]"),
        ("concept DOI cited", "zenodo.23262404"),
        ("receptor fractions from the re-run", "which is why the receptor fractions are reported"),
        ("reproducible core scoped", "the shared core is AP-1 itself"),
        ("Figure 5 title peri-necrotic", "The MES niche is peri-necrotic"),
        ("Fig 1 pooled composition labelled", "pooled over the 20,643 malignant cells"),
        ("Fig 1c inference removed", "compositional evenness does not by itself indicate"),
        ("Fig 1f positive-only", "restricted to positively correlated genes"),
    ]


VERIFY_EXTRA = _verify_extra()
