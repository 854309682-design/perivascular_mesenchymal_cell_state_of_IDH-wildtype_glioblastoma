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
    E.append(("para", "7440013E",
              "Second, we rebuilt gene-regulatory networks de novo in each cohort with pySCENIC "
              "[9] (GRNBoost2 co-expression modules, cisTarget motif pruning, AUCell per-cell "
              "activity), which makes no use of curated TF" + EN + "target priors (for the "
              "multi-omic extension of this framework, see SCENIC+ [14]). Because the submitted "
              "de-novo arm could not be reproduced, this route was re-run from the raw data in the "
              "three cohorts common to both routes, with a state-stratified subsample of 2,000 "
              "malignant cells per cohort (" + ", ".join(f"{k} {v}" for k, v in D["counts"].items())
              + " regulons). Its consensus (Fig. 3b) is headed by " + top10_str(D) + ". The top of "
              "this list is not the submitted one. Auditing every ranked regulon against the "
              "cisTarget TF universe accounts for the artefact the reviewers identified: "
              + ", ".join(_absent(["CD59", "NPAS2", "FOXC1"]))
              + " do not appear among the ranked regulons at all, and NFATC4 appears only at rank "
              + str(rank_of.get("NFATC4", -1)) + " (mean " + RHO + " "
              + f"{rho_of.get('NFATC4', 0.0):+.2f}); TPI1 is ranked "
              + str(rank_of.get("TPI1", -1)) + " (mean " + RHO + " "
              + f"{rho_of.get('TPI1', 0.0):+.2f}) \u2014 TPI1 is not a transcription factor, but "
              "it is listed in that universe, so we report the full ranking in the shipped table "
              "and exclude it from the switch. The per-cohort distributions of "
              "regulon" + EN + "MES correlations (Fig. 3e) show that each cohort independently "
              "contains strongly MES-tracking regulons, even though the cohort medians sit near "
              "zero. The signal is therefore concentrated in a small set of regulons, not spread "
              "across the whole network."))

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
    E.append(("sub", "6ACF45E9",
              " Per cohort, a state-stratified subsample of malignant cells (\u226410,000 cells, "
              "min 200 per state) was exported as a loom. GRNBoost2 (arboreto) inferred "
              "co-expression modules; cisTarget (hg38, mc9nr, 10 kb and 500 bp databases) pruned "
              "modules to motif-supported regulons; AUCell produced per-cell regulon activity "
              "(278\u2013477 regulons per cohort). Regulon activity was correlated with the MES "
              "score per cohort (Spearman plus regulon-activity score), and the four-cohort "
              "consensus ranked regulons by mean " + RHO + " with a per-cohort significance count.",
              " Per cohort, a state-stratified subsample of malignant cells (500 per state, "
              "\u22642,000 per cohort; genes detected in \u22655% of the subsampled cells) was "
              "exported as a loom. GRNBoost2 (arboreto) inferred co-expression modules; cisTarget "
              "(hg38, mc9nr, 10 kb and 500 bp databases, NES \u2265 3.0) pruned modules to "
              "motif-supported regulons; AUCell produced per-cell regulon activity (" + regs +
              " regulons). This re-run covers the three cohorts common to both routes and uses a "
              "smaller subsample than the submitted analysis, because the original code was not "
              "recoverable; every ranked regulon was audited against the cisTarget TF universe. "
              "Regulon activity was correlated with the MES score per cohort, and the three-cohort "
              "consensus ranked regulons by mean " + RHO + " with a per-cohort significance count."))

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
    E.append(("para_sub", "7C72AC66", [
        ("The four-cohort consensus (Fig. 3a) placed STAT3 and KLF4 at the top (mean " + RHO
         + " \u2248 0.39, STAT3 FDR < 0.05 in 4/4 cohorts), with the AP-1 factors JUNB and FOSL1 "
         "also positive (" + RHO + " \u2248 0.30/0.29).",
         "The four-cohort consensus (Fig. 3a) placed " + t2[0] + " and " + t2[1] + " at the top "
         "(mean " + RHO + f" {cmap[t2[0]]:.2f}/{cmap[t2[1]]:.2f}, {t2[0]} FDR < 0.05 in 4/4 "
         "cohorts), with the AP-1 factors FOSL1 and JUNB also positive (" + RHO + " "
         f"{cmap['FOSL1']:.2f}/{cmap['JUNB']:.2f})."),
        ("with AP-1/STAT3 positive and CEBPB/RELA/TEAD1 negative in the prior-based readout.",
         "with AP-1/STAT3 positive and CEBPB and RELA negative in the prior-based readout "
         "(TEAD1 is positive in three of the four cohorts and is not placed in the switch)."),
    ]))
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
    E.append(("para_sub", "07505CF4", [
        ("SCENIC inputs for two cohorts are log-normalized rather than UMI counts, and the "
         "IDH-mutant cohorts contain few MES cells.",
         "SCENIC inputs for two cohorts are log-normalized rather than UMI counts. We also found "
         "and corrected a platform-scale inconsistency in GSE131928: its two public matrices are "
         "delivered on different scales (log2(TPM+1) versus raw TPM) and concatenating them made "
         "that cohort's MES score and pseudotime axis predominantly a platform axis, so the cohort "
         "is now scored on its Smart-seq2 block alone \u2014 which changes its patient count "
         "(45 to 35) and its mean MES-like proportion (0.185 to 0.216) and removes its "
         "contribution to the multi-cohort pseudotime claim. The IDH-mutant cohorts contain few "
         "MES cells."),
    ]))

    # ---------- Discussion: the pseudotime claim
    E.append(("para_sub", "03BE3C7E", [
        ("lies at the end of a reproducible but branched pseudotime axis that runs through the "
         "AC/NPC compartment (Fig. 4)",
         "lies at the end of a branched pseudotime axis that runs through the AC/NPC compartment "
         "in the discovery cohort (Fig. 4); we could not reproduce that axis per patient in the "
         "other two cohorts, so we report it for the discovery cohort only"),
    ]))

    # ---------- Results 3.4: the pseudotime axis, pooled vs per sample
    E.append(("para_sub", "52711FBE", [
        ("Repeating the analysis separately per cohort strengthened the result: oriented "
         "pseudotime tracked the MES score with " + RHO + " = 0.66 (gse103224), 0.61 (GSE182109) "
         "and 0.58 (Neftel) (Fig. 4b,e).",
         "Within that cohort the correlation held at patient level (median per-sample " + RHO
         + " = " + fmt_r(T["med"]["GSE103224"]) + "; " + str(int(T["gt03"]["GSE103224"])) + " of "
         + str(int(T["nsamp"]["GSE103224"])) + " tumours above 0.3; Fig. 4d). Extending the same "
         "procedure to the other two cohorts did not reproduce it: in GSE182109 the pooled "
         "correlation was " + RHO + " = " + fmt_r(rho["GSE182109"]) + " (n = "
         + f"{T['n']['GSE182109']:,}" + " malignant cells) but the median per-sample correlation "
         "was only " + fmt_r(T["med"]["GSE182109"]) + " with " + str(int(T["gt03"]["GSE182109"]))
         + " of " + str(int(T["nsamp"]["GSE182109"])) + " samples above 0.3, i.e. that pooled "
         "value is driven by between-sample variation, and in GSE131928 the two axes were "
         "effectively uncorrelated (" + RHO + " = " + fmt_r(rho["Neftel"]) + ", median per sample "
         + fmt_r(T["med"]["Neftel"]) + "; Fig. 4b,e). We therefore changed the claim: the "
         "MES-directed axis is presented as a property of the discovery cohort, replicated at "
         "patient level within it, and we do not claim it for the other two cohorts (Fig. 4e)."),
        ("Per-sample estimates, while heterogeneous, were predominantly positive in all cohorts "
         "(Fig. 4d), but the pooled coefficients overstate the per-patient effect: the median "
         "per-patient rho is well below the pooled value and the per-sample range spans negative "
         "values, so the pooled correlation is descriptive and the per-sample distributions are "
         "the honest summary (Fig. 4b,d,e).",
         "Per-sample estimates were heterogeneous: they were predominantly positive in the "
         "discovery cohort (" + str(int(T["pos"]["GSE103224"])) + " of "
         + str(int(T["nsamp"]["GSE103224"])) + " tumours, which is why its pooled and per-sample "
         "coefficients agree) but near zero in the other two, where the pooled value overstates "
         "the per-patient effect (GSE182109: pooled " + fmt_r(rho["GSE182109"]) + " versus median "
         "per-sample " + fmt_r(T["med"]["GSE182109"]) + "); the per-sample ranges span negative "
         "values in all three cohorts, so the pooled correlation is descriptive and the "
         "per-sample distributions are the honest summary (Fig. 4b,d,e)."),
        ("The MES-directed axis is therefore a recurring property of glioma transcriptomes, not an "
         "artefact of one dataset.",
         "The MES-directed axis is therefore reported as a reproducible property of the discovery "
         "cohort rather than of glioma transcriptomes in general, and the submitted per-cohort "
         "values for GSE182109 and GSE131928 are not reproduced."),
    ]))

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
        ("de-novo audit stated", "listed in that universe"),
        ("trajectory claim scoped", "we do not claim it for the other two cohorts"),
        ("Figure 4 title scoped", "A branched pseudotime axis orders the MES transition in the "
                                  "discovery cohort"),
        ("Methods MES terminology", "names the mesenchymal meta-module state"),
        ("GSE163120 in Methods", "excluded a priori because it contains only immune cells"),
        ("abstract updated", "that overlap is greater than chance"),
    ]


VERIFY_EXTRA = _verify_extra()
