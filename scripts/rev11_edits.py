"""rev11_edits — the revision edit table for scripts/rev11_manuscript_revision.py.

Anchors are copied **verbatim from the raw document.xml** (the manuscript uses en dashes
and unicode kappa/approximately-equal signs, so anchors must not be "normalised").
`("para", pid, text)` is only valid for paragraphs with a single non-empty run.
"""

NEW_REFS = [
    ("Carro M.S. et al. The transcriptional network for mesenchymal transformation of brain "
     "tumours. Nature 2010 (PMID 20032975; C/EBPbeta and STAT3 as master regulators of "
     "mesenchymal transformation)."),
    ("Bhat K.P.L. et al. Mesenchymal differentiation mediated by NF-kB promotes radiation "
     "resistance in glioblastoma. Cancer Cell 2013 (PMID 24209607)."),
    ("Wang Q. et al. Tumor evolution of glioma-intrinsic gene expression subtypes associates "
     "with immunological changes in the microenvironment. Cancer Cell 2017 (PMID 28697342; "
     "mesenchymal transition in bulk tumours tracks microenvironmental rather than genetic "
     "evolution, and the origin of the 'bulk MES = stromal content' caveat tested here)."),
    ("Minata M. et al. Phenotypic plasticity of invasive edge glioma stem-like cells in response "
     "to ionizing radiation. Cell Rep 2019 (PMID 30759398; radiation-induced, "
     "C/EBPbeta/STAT3-driven reversible mesenchymal plasticity)."),
    # ---------- Methods + Figure 5 legend: gene-set provenance and the rebuilt panels ----------
    ("sub", "4DF903D0",
     " Ligand activity against the 73-gene MES-transition program was computed",
     " Ligand activity against the 128-gene MES-transition programme was computed"),
    ("sub", "50242B74",
     ' (a) IvyGAP MES score (z per tumour) across anatomic regions (LE leading edge, IT infiltrating tumour, CT cellular tumour, MVP microvascular proliferation, HBV hyperplastic blood vessels, PNZ peri-necrotic zone, PAN pseudopalisading around necrosis); Kruskal–Wallis P = 8.7 × 10⁻²⁰; MVP and HBV significantly above LE (FDR ≈ 1 × 10⁻⁷ and 1.8 × 10⁻⁵), PNZ/PAN not significant. (b) Hypoxia score across the same regions; unlike MES, hypoxia also trends up near necrosis, showing the MES signal is not a hypoxia proxy. (c) Mixed-effects model region coefficients (reference LE, random intercept per tumour): MVP +0.28, HBV +0.27. (d) NicheNet ligand activity against the 73-gene MES-transition program (top: CSF1, TFPI, ITM2B, FGF5, BMP3, IL6, BMP2, EPGN, ANGPT4). (e) Percentage of malignant GSE182109 cells expressing each receptor (IL6ST 55%, LRP1 56%, EGFR 38%, BMPR2 35%, FGFR1 27%; CSF1R 7%). (f) Forest plot of the MES-score hazard ratio for overall survival (multivariable Cox, MES plus age and IDH): TCGA HR 1.22, CGGA-325 HR 2.28, CGGA-693 HR 1.62, CGGA pooled, and random-effects meta-analysis HR 1.63 (95% CI 1.18–2.26, P = 0.003; I² = 89%).',
     ' (a) IvyGAP anatomic regions: MES (MES1+MES2), MES1 and MES2 scores (z-normalised per tumour) in the 122 laser-microdissected samples; points are samples, bars are medians; n samples/tumours: LE 19/8, IT 24/8, CT 30/10, MVP 25/9, PAN 24/8. The FACS-sorted cancer-stem-cell-cluster samples (CThbv-/CTpnz-) are excluded because they are not anatomic regions. (b) Hypoxia score of the same regions (the 15-gene HIF-target set of the Methods). (c) Mixed-effects model region coefficients (reference LE, random intercept per tumour), unadjusted and adjusted for vascular (endothelial+pericyte) and myeloid marker scores; the MVP effect is halved by adjustment while the PAN effect persists. (d) Signature-matrix (NNLS) deconvolution against the GSE103224 single-cell reference: malignant fraction (black) and MES share of malignant cells (red) per region. (e) NicheNet ligand activity against the 128-gene MES-transition programme (top 15). (f) Cox hazard ratios per SD of the MES score, unadjusted and grade-adjusted, in TCGA and CGGA-325/693; Supplementary Fig. S2a gives the pooled random-effects estimate.'),
]

EDITS = [
    # ---------- M1: title / framing --------------------------------------------------------
    ("para", "47F0C2C8",
     "A compact AP-1/C-EBP regulatory program marks the plastic, peri-necrotic mesenchymal "
     "state of IDH-wildtype glioblastoma"),
    ("sub", "61AE3870",
     " glioblastoma; mesenchymal transition; cell-state plasticity; gene regulatory network; perivascular niche",
     " glioblastoma; mesenchymal transition; cell-state plasticity; gene regulatory network; "
     "peri-necrotic niche; prognostic score"),

    # ---------- abstract: fully rewritten (single-run paragraph) ---------------------------
    ("para", "31C93A94",
     "Glioma cells reversibly adopt a mesenchymal (MES) state linked to invasion, therapy "
     "resistance and poor outcome, but whether MES transition is a fixed subclonal fate or a "
     "plastic, microenvironment-driven program, and which transcription factors (TFs) control "
     "it, remains unresolved. We analysed five independent single-cell RNA-seq cohorts of human "
     "glioma (GSE103224, GSE131928, GSE182109, GSE70630, GSE89567), the anatomic-region-resolved "
     "IvyGAP resource, and three clinical bulk cohorts (TCGA, CGGA-325, CGGA-693). The MES "
     "fraction varied continuously across patients, behaving as a plastic rather than a discrete "
     "state, and the program was essentially confined to IDH-wildtype tumours (mean MES-like "
     "proportion 0.267 versus 0.008; P = 1.2 x 10-7). Primary and recurrent tumours did not "
     "differ (0.507 versus 0.473; two-sided P = 0.74; n = 11 versus 5), which at this sample size "
     "is absence of evidence rather than evidence of absence. Two regulatory-inference routes "
     "were compared: prior-based TF-activity scoring (DoRothEA) and a prior-independent pySCENIC "
     "reconstruction re-run from the raw data. The prior route reproduced across four cohorts the "
     "known AP-1/C-EBP axis (mean rho 0.43 STAT3, 0.41 KLF4, 0.35 JUNB, 0.32 FOSL1; FDR < 0.05 "
     "in 4/4) and placed CEBPB (-0.22), RELA (-0.27) and TWIST1 (-0.16) consistently negative; "
     "TWIST1 was therefore not supported as a switch component and is not claimed here. "
     "Patient-level inference attenuated every correlation (JUNB 0.35 to 0.05), and the two "
     "routes overlapped in only a handful of factors. Anatomically, the MES program was highest "
     "in the pseudopalisading/peri-necrotic regions of IvyGAP (z = +1.41 versus -1.11 at the "
     "leading edge) and not in microvascular-proliferation regions (+0.06), whose apparent signal "
     "was confined to the hypoxia-linked MES1 submodule and halved after adjustment for vascular "
     "and myeloid content. A MES score was prognostic in CGGA-325 (HR 1.70 per SD, 95% CI "
     "1.46-1.98; grade-adjusted 1.45) and CGGA-693 (1.41, 1.26-1.57; grade-adjusted 1.23) but "
     "only marginally in TCGA (1.15, 1.02-1.30; grade-adjusted 1.10, 0.98-1.24) and not after "
     "adjustment for stromal/myeloid content. All findings are computational and associative."),

    # ---------- M8 -------------------------------------------------------------------------
    ("sub", "49910F62",
     "We pre-registered three hypotheses (see Methods): ",
     "We specified three hypotheses a priori in an internal analysis plan (see Methods): "),

    # ---------- canonical prior art --------------------------------------------------------
    ("sub", "6109BE8D",
     "an inflammatory microenvironment, and adverse outcome [2,3].",
     "an inflammatory microenvironment, and adverse outcome [1,28]."),
    ("sub", "6109BE8D",
     "Three questions, however, remain open.",
     "The transcriptional basis of this state was first defined by C/EBPbeta-STAT3 as master "
     "regulators of mesenchymal transformation [26], extended by NF-kB-driven differentiation and "
     "radioresistance [27], linked in bulk tumours to microenvironmental rather than genetic "
     "evolution [28], and shown to be inducible and reversible after irradiation in preclinical "
     "models [29]. Three questions, however, remain open."),
    ("sub", "397DD2F9",
     "We report that the MES state is a plastic, IDH-wildtype program that is not selected for "
     "by recurrence; that a compact C/EBP–AP-1–TWIST1 switch (with NF-κB and hypoxia-linked "
     "arms) reproducibly tracks MES transition across all four cohorts and both inference "
     "routes; that the program localizes to the perivascular/angiogenic niche rather than the "
     "peri-necrotic zone; and that a MES-transition score is an externally validated, adverse "
     "prognostic marker.",
     "We report that the MES state is a plastic, IDH-wildtype program with no evidence of "
     "enrichment at recurrence; that a compact AP-1/C-EBP regulatory core (with a STAT3 arm and "
     "a hypoxia-linked MES1 component) reproduces across cohorts in the prior-based route and is "
     "recovered de novo, whereas TWIST1 is not supported by either route; that the program is "
     "highest in the pseudopalisading/peri-necrotic region rather than in microvascular-"
     "proliferation regions, where the apparent enrichment is explained by the hypoxia-linked "
     "MES1 submodule and by vascular and myeloid content; and that a MES score carries external "
     "prognostic information that is attenuated after adjustment for grade and for stromal "
     "content."),

    # ---------- Methods -------------------------------------------------------------------
    ("sub", "092B3563",
     " Public datasets: GSE103224 [4], GSE131928 [1], GSE182109 [5], GSE70630 [10], GSE89567 "
     "[6], IvyGAP GSE107559 [7], CGGA mRNAseq-325 and mRNAseq-693, and TCGA-GBM (GDC "
     "STAR-counts, per-case BCR clinical XML, masked MAF for IDH).",
     " Public datasets: GSE103224 [4], GSE131928 (Neftel 2019) [1], GSE182109 [5], GSE70630 "
     "(Tirosh 2016; IDH-mutant oligodendroglioma) [10], GSE89567 (Venteicher 2017; IDH-mutant "
     "astrocytoma) [6], IvyGAP GSE107559 [7], CGGA mRNAseq-325 and mRNAseq-693, and TCGA (GDC "
     "STAR-counts and masked MAF for IDH; both the GBM and the LGG project are used, so the "
     "cohort is reported as a lower-grade-inclusive glioma cohort and every survival model is "
     "shown with and without grade adjustment)."),
    ("sub", "092B3563",
     " All analyses used Python 3.12 (scanpy 1.12, harmonypy 2.0, infercnvpy 0.6, decoupler 2.2, "
     "lifelines 0.30, statsmodels 0.14), pySCENIC in a dedicated conda environment, and R 4.6; "
     "seed = 42 throughout.",
     " All analyses used Python 3.11-3.12 (scanpy 1.11-1.12, harmonypy 2.0, infercnvpy 0.6, "
     "statsmodels 0.15, scikit-learn 1.9, scipy 1.17), pySCENIC 0.12.1 with ctxcore 0.2.0 and "
     "arboreto (GRNBoost2), and R 4.6 (survival, survminer, metafor); seed = 42 throughout. The "
     "multivariate linear model of the DoRothEA prior and every proportional-hazards diagnostic "
     "reported here were recomputed with statsmodels and R, and are cross-checked against the "
     "archived activity matrices (median per-TF Spearman rho = 0.998)."),
    ("sub", "368C2319",
     "per-patient MES-like proportions were compared with a one-sided Mann–Whitney test.",
     "per-patient MES-like proportions were compared with a two-sided Mann–Whitney test, and the "
     "difference in means is reported with a bootstrap 95% confidence interval and Cliff's delta "
     "as effect size."),
    ("sub", "39E68009",
     " DoRothEA regulons (decoupler ",
     " DoRothEA regulons (confidence levels A-C, 262 TFs) were scored per cell with a documented "
     "multivariate linear model that reproduces the archived activity matrix at a median "
     "per-TF Spearman rho of 0.998; the prior was read with "),
    ("sub", "39E68009",
     ", MLM) were scored per cell in each cohort; Spearman correlation of TF activity with the "
     "MES score was computed per cohort, FDR-corrected, and averaged across cohorts for the consensus.",
     " (multivariate linear model). Spearman correlation of TF activity with the MES score was "
     "computed per cohort, FDR-corrected, and averaged across cohorts for the consensus. "
     "Correlations computed across cells treat cells as independent observations and are "
     "therefore descriptive; per-patient estimates and per-sample distributions are reported "
     "alongside them (Fig. 3f)."),
    ("sub", "545338F9",
     " Per-region MES and hypoxia scores (z-normalized per tumour) were compared across the "
     "seven anatomic regions by Kruskal–Wallis and region-versus-LE Mann–Whitney tests "
     "(Benjamini–Hochberg FDR), and by a linear mixed model with random intercept per tumour "
     "(statsmodels).",
     " Only the 122 laser-microdissected anatomic-structure samples (10 tumours; CT n = 30, "
     "CTmvp n = 25, CTpan n = 24, IT n = 24, LE n = 19, five structures) were analysed; the 148 "
     "FACS-sorted cancer-stem-cell-cluster samples in the same series carry CThbv-/CTpnz- labels "
     "but are not anatomic regions and were excluded. MES (MES1+MES2), MES1, MES2 and hypoxia "
     "scores were z-normalised per tumour and compared across regions by Kruskal–Wallis and "
     "region-versus-LE Mann–Whitney tests with Benjamini–Hochberg FDR, and by a linear mixed "
     "model with a random intercept per tumour. To control for the cellular composition of bulk "
     "regions, the same models were refitted with endothelial+pericyte and myeloid marker scores "
     "as covariates, and a signature-matrix (non-negative least squares) deconvolution against "
     "the malignant and non-malignant GSE103224 single-cell reference was used to estimate the "
     "malignant fraction and the MES share of malignant cells per region. Hypoxia was scored "
     "with a fixed 15-gene HIF-target set (VEGFA, SLC2A1, LDHA, PGK1, PDK1, HK2, BNIP3, NDRG1, "
     "ENO1, PDK3, TPI1, P4HA1, ADM, MRPS17, ANKZF1)."),
    ("sub", "64FDFA0F",
     " The MES score (73-gene program, mean z-score per sample) was tested in Cox "
     "proportional-hazards models of overall survival adjusting for age and IDH (TCGA, n = 939; "
     "survival from per-case clinical XML ",
     " The bulk MES score is the mean z-score of the Neftel MES1+MES2 meta-module (95 of 100 "
     "genes found in STAR counts); the 128-gene MES-transition programme is reported as an "
     "independent benchmark. Patients, not samples, are the unit of analysis: TCGA has 809 "
     "patients with expression and 737 with survival follow-up. Events are deaths (vital status "
     "DECEASED) and time is the reported follow-up in months converted to days, cross-checked "
     "against the per-case clinical fields "),
    ("sub", "64FDFA0F",
     ") and in CGGA-325/693 with and without WHO-grade adjustment. Cohort log-HRs were pooled by "
     "fixed- and random-effects (DerSimonian–Laird) meta-analysis.",
     "). Models give hazard ratios per standard deviation of the z-scored MES score with 95% "
     "confidence intervals, proportional hazards were checked with scaled Schoenfeld residuals "
     "(cox.zph), Kaplan–Meier curves use a stated median cut-off with log-rank P values and "
     "numbers at risk, and the random-effects estimate is reported as primary together with the "
     "between-cohort range and I-squared. The same models were repeated with the deconvolved "
     "malignant fraction, with stromal/myeloid marker scores, and with WHO grade as covariates. "
     "CGGA-325/693 were analysed with the same rule, with WHO grade parsed from the CGGA "
     "Roman-numeral notation and IDH/MGMT taken from the CGGA clinical table."),

    # ---------- Results: entropy (minor 1 / DA7) ------------------------------------------
    ("sub", "3893E80F",
     "This heterogeneity was accompanied by graded state ambiguity: tumours with intermediate "
     "MES fractions showed the highest Shannon entropy of their four-state composition, whereas "
     "the extremes were compositionally pure (Fig. 1c). This is the pattern expected of cells "
     "transiting between states rather than of discrete, stably coexisting subclones.",
     "This heterogeneity was accompanied by graded state ambiguity: the Shannon entropy of the "
     "four-state composition was highest at low-to-intermediate MES fractions and lowest at the "
     "extremes (Fig. 1c; n = 8 tumours, so this panel is descriptive only). Compositional "
     "evenness does not by itself measure transitions — four stably coexisting states give the "
     "same entropy — so we rest the plasticity argument on the continuous single-cell gradients "
     "of Fig. 1d,e and on the pseudotime analysis below."),

    # ---------- Results: IDH confound (minor 13 / stats 2.5) -------------------------------
    ("sub", "050C4C6F",
     "The MES program is therefore essentially an IDH-wildtype (GBM) phenotype, which addresses "
     "the genotype arm of H1.",
     "The comparison is confounded by design: the IDH-mutant series are oligodendroglioma and "
     "IDH-mutant astrocytoma, so genotype is fully confounded with histology, grade, cohort and "
     "platform, and a pooled Mann–Whitney test ignores that stratification. A cohort-stratified "
     "van Elteren test and the per-cohort effects are reported in Supplementary Fig. S2f; the "
     "direction is unchanged in every cohort. The MES program is therefore essentially an "
     "IDH-wildtype (GBM) phenotype, which addresses the genotype arm of H1."),

    # ---------- Results: recurrence (M4) --------------------------------------------------
    ("sub", "6D6FFB14",
     "The treatment arm of H1, by contrast, was rejected.",
     "The treatment arm of H1 received no support — which, at this sample size, is absence of "
     "evidence rather than evidence of absence."),
    ("sub", "6D6FFB14",
     "The mean MES-like proportion was 0.507 in primary versus 0.473 in recurrent disease "
     "(one-sided Mann–Whitney P = 0.67; Fig. 2d), and the per-patient values overlapped broadly "
     "(Fig. 2e). Recurrence, and by proxy therapy exposure, does not enrich the MES compartment.",
     "The mean MES-like proportion was 0.507 in primary versus 0.473 in recurrent disease "
     "(difference -0.034, bootstrap 95% CI -0.354 to +0.285; two-sided Mann–Whitney P = 0.74, "
     "Cliff's delta -0.13; Fig. 2d), and the per-patient values overlapped broadly (Fig. 2e). "
     "With five recurrent tumours this test can only exclude very large enrichment effects, so we "
     "read it as no evidence of MES enrichment at recurrence rather than as a demonstration that "
     "recurrence does not enrich the MES compartment."),
    ("sub", "6D6FFB14",
     "This agrees with a 59-patient longitudinal single-nucleus study that found no consistent "
     "malignant-state enrichment at recurrence [12], and with bulk longitudinal analyses linking "
     "mesenchymal transition to microenvironmental rather than genetic evolution [2].",
     "The concordant external evidence therefore carries the inferential weight: a 59-patient "
     "longitudinal single-nucleus study found no consistent malignant-state enrichment at "
     "recurrence [12], and bulk longitudinal analyses link mesenchymal transition to "
     "microenvironmental rather than genetic evolution [2,28]."),

    # ---------- Results: the two-route comparison (M2 c/d/e) ------------------------------
    ("sub", "21523EDC",
     "The direct comparison of the two routes (Fig. 3f) tests robustness. Across all TFs scored "
     "by both methods, the mean MES correlation of DoRothEA and of de-novo SCENIC agreed in "
     "direction for the leading factors, with FOSL1, CEBPD, CEBPB, JUNB and STAT3 occupying the "
     "concordant high-high quadrant. The top-30 lists of the two approaches share exactly these "
     "AP-1/C/EBP factors (Supplementary Fig. S2e).",
     "The direct comparison of the two routes (Fig. 3c,f) tests robustness, and we now quantify "
     "rather than assert it. The two top-30 lists overlap in only a small AP-1/C-EBP set, and "
     "that overlap is tested formally against the scored-TF universe with a Fisher exact test "
     "(Supplementary Fig. S2e), not described as general agreement. Both routes place AP-1 and "
     "C/EBP family members among the strongest factors, but they disagree in ways that matter: "
     "STAT3 leads the prior-based route (mean rho 0.43, FDR < 0.05 in 4/4 cohorts) while the "
     "de-novo route ranks it lower, and CEBPB is consistently negative in the prior-based route "
     "(mean rho -0.22, FDR < 0.05 in 4/4) while the de-novo route ranks it positively."),
    ("sub", "21523EDC",
     "Because the de-novo route cannot inherit prior biases, we take the consensus switch to be ",
     "Because the de-novo route cannot inherit prior biases, we take the reproducible consensus "
     "to be an AP-1/C-EBP core, namely "),
    ("sub", "21523EDC",
     "CEBPD/CEBPB, AP-1 (FOSL1/FOSL2/JUN/JUNB) and TWIST1, with NF-κB (RELB/RELA), "
     "hypoxia-linked (ATF3, BHLHE40) and STAT3 arms",
     "AP-1 (FOSL1/FOSL2/JUN/JUNB) together with CEBPD, with a STAT3 arm and a hypoxia-linked arm"),
    ("sub", "21523EDC",
     ". This supports H2 with a prior-independent core.",
     ". Two factors that the submitted version placed in the switch are not supported by the "
     "re-run and are excluded here: TWIST1, whose activity is negative in all four cohorts in the "
     "prior-based route (mean rho -0.16, FDR < 0.05 in 4/4; positive in only 30% of patients) and "
     "which never entered the shared set, and CEBPB, whose sign is not stable between routes. "
     "This supports H2 only in the weakened form stated here: a small, reproducible AP-1/C-EBP "
     "core, not a compact three-pillar switch."),
    ("sub", "62E53675",
     "(Fig. 3, continued.) The switch factors are themselves embedded in the continuous MES "
     "gradient of Fig. 1 rather than marking a binary compartment.",
     "The core factors are themselves embedded in the continuous MES gradient of Fig. 1 rather "
     "than marking a binary compartment."),
    ("sub", "62E53675",
     "consistent with the plasticity observed in Sections 3.1–3.2.",
     "consistent with the plasticity described in the two preceding sections."),

    # ---------- Results: trajectory + circularity (minor 3, 4) ----------------------------
    ("sub", "52711FBE",
     "Per-sample estimates, while heterogeneous, were predominantly positive in all cohorts "
     "(Fig. 4d).",
     "Per-sample estimates, while heterogeneous, were predominantly positive in all cohorts "
     "(Fig. 4d), but the pooled coefficients overstate the per-patient effect: the median "
     "per-patient rho is well below the pooled value and the per-sample range spans negative "
     "values, so the pooled correlation is descriptive and the per-sample distributions are the "
     "honest summary (Fig. 4b,d,e). The 128-gene programme is derived from this axis in the "
     "discovery cohort and the state composition along the same axis is then shown in the same "
     "cohort (Fig. 1f, Fig. 4f), so those two panels are not independent evidence; the programme "
     "is also restricted to positively correlated genes, since anti-correlated genes do not "
     "belong in a MES programme."),
    ("sub", "69892ECF",
     "The topology of the axis, however, is branched rather than linear.",
     "The topology of the axis is branched in the discovery cohort, and we do not generalise it: "
     "PAGA was computed on four coarse state labels in GSE103224 only, so the branched topology "
     "is presented as a property of that cohort rather than of glioma in general."),

    # ---------- Results: IvyGAP (M3) ------------------------------------------------------
    ("sub", "441088FC",
     "We scored the MES program per region. The MES score peaked in the two vascular regions "
     "(MVP z ≈ +0.19, HBV z ≈ +0.18) against the leading edge (−0.08), with peri-necrotic regions "
     "at or below baseline (Kruskal–Wallis P = 8.7 × 10⁻²⁰; region-wise Mann–Whitney versus LE: "
     "MVP FDR ≈ 1 × 10⁻⁷, HBV FDR ≈ 1.8 × 10⁻⁵, PNZ/PAN not significant; Fig. 5a).",
     "We scored the MES program per region in the 122 anatomic samples (10 tumours; Fig. 5a). The "
     "programme was highest in the pseudopalisading/peri-necrotic region (CTpan z = +1.41 versus "
     "LE -1.11; 24 samples, 8 tumours; Mann–Whitney versus LE FDR = 1.1 × 10⁻⁷), with "
     "intermediate values in the infiltrating and cellular tumour regions and almost no elevation "
     "in microvascular proliferation (CTmvp z = +0.06; 25 samples, 9 tumours; FDR = "
     "2.9 × 10⁻⁵; Kruskal–Wallis P = 1.7 × 10⁻¹⁴ across the five structures). Splitting the "
     "meta-modules separates the two signals: the hypoxia-linked MES1 submodule is elevated in "
     "both CTmvp (+0.53) and CTpan (+0.76), whereas MES2, the mesenchymal module proper, is "
     "elevated only in CTpan (+1.61, FDR = 4.8 × 10⁻⁷) and sits at or below baseline in CTmvp "
     "(-0.31, not significant). Two features of the submitted analysis produced the opposite "
     "impression. First, it pooled the hypoxia-linked MES1 with MES2, importing a "
     "vascular-hypoxic signal into the mesenchymal score. Second, bulk regional MES scores track "
     "the cellular composition of the region: MES correlates with myeloid (rho = +0.37) and "
     "pericyte (+0.21) content and inversely with oligodendrocyte content (-0.60), and the "
     "vascular regions are compositionally extreme, with a deconvolved malignant fraction of "
     "0.37 in CTmvp against 0.98 in CTpan (Fig. 5d). Adjusting for vascular and myeloid content "
     "halves the CTmvp effect (mixed model +1.24 to +0.57) while the CTpan effect is essentially "
     "unchanged (+2.73 to +2.53; Fig. 5c). Deconvolution gives the same answer at the level of "
     "malignant cells: the MES share of malignant cells is highest in CTpan (0.64, FDR = "
     "3.6 × 10⁻⁷) and lowest in CTmvp (0.18, FDR = 3.5 × 10⁻³; Fig. 5d)."),
    ("sub", "441088FC",
     "A mixed-effects model with tumour as a random effect confirmed the region effect "
     "independently of patient (MVP +0.28, HBV +0.27; Fig. 5c).",
     "A mixed-effects model with tumour as a random effect confirms that the region effect is not "
     "a patient effect (Fig. 5c)."),
    ("sub", "441088FC",
     "Hypoxia scoring of the same regions (Fig. 5b) showed that hypoxia is, as expected, also "
     "elevated around MVP/HBV but, unlike MES, trends up in the necrotic-adjacent regions as "
     "well. The MES program therefore does not simply shadow hypoxia. We thus refine H3: the MES "
     "niche is ",
     "Hypoxia scoring of the same regions (Fig. 5b; the 15-gene HIF-target set listed in the "
     "Methods) is elevated in the pseudopalisading/peri-necrotic region and, to a lesser extent, "
     "in microvascular proliferation, so part but not all of the regional MES signal tracks "
     "hypoxia; the MES2 elevation in CTpan survives both the hypoxia association and the "
     "composition adjustment, whereas the CTmvp signal does not. We thus revise H3 against our "
     "own earlier reading: within these data the MES programme is not perivascular. It is highest "
     "in the pseudopalisading/peri-necrotic compartment, and the previously reported vascular "
     "enrichment is explained by the hypoxia-linked MES1 submodule together with the vascular and "
     "myeloid content of microvascular-proliferation regions. The niche is "),
    ("sub", "441088FC",
     "This agrees with spatial-transcriptomic analyses of GBM reporting pro-inflammatory "
     "perivascular regions [16], mesenchymal tumour cells within monocyte-rich domains [17], and "
     "perivascular malignant–vascular–immune associations in high-resolution atlases [18,19].",
     "This residual association between the MES programme and vascular, myeloid-rich microdomains "
     "is consistent with spatial-transcriptomic analyses reporting pro-inflammatory perivascular "
     "regions [16], mesenchymal tumour cells within monocyte-rich domains [17], and perivascular "
     "malignant–vascular–immune associations in high-resolution atlases [18,19], while the "
     "programme itself is highest in the peri-necrotic compartment, in keeping with classical "
     "hypoxic-niche models [1,25]."),

    # ---------- Results: NicheNet (minor 6, items 2.8-2.9) --------------------------------
    ("sub", "3ED7EBF0",
     "which ligands best predict the 73-gene MES-transition program",
     "which ligands best predict the MES-transition programme"),
    ("sub", "3ED7EBF0",
     "The top-ranked ligands (Fig. 5d) were CSF1, TFPI, ITM2B, FGF5, BMP3, IL6, BMP2, EPGN and "
     "ANGPT4, paracrine factors characteristic of the vascular–myeloid compartment.",
     "Ligand activity was computed as the correlation between the ligand's regulatory-potential "
     "vector in the NicheNet prior and the programme over the shared genes (124 of 128 genes "
     "present). The top-ranked ligands (Fig. 5e) were CSF1, TFPI, ITM2B, FGF5, BMP3, LEFTY2, IL6, "
     "BMP2, PSPN, LEAP2, SLITRK1, SPON1, EPGN, ANGPT4 and VIP. We report the full ranked list "
     "because several entries are not conventional signalling ligands, which points to "
     "prior-matrix artefacts, and because the activity values are small (Pearson r = 0.11–0.14): "
     "this analysis nominates candidates, it does not establish mechanism."),
    ("sub", "3ED7EBF0",
     "Their cognate receptors are expressed on the malignant cells themselves (Fig. 5e): IL6ST "
     "(55% of malignant cells), LRP1 (56%), ITGB1, EGFR (38%), BMPR2 (35%), CD44, FGFR1 (27%), "
     "TGFBR1, whereas CSF1R is nearly absent (7%). This pattern places the source of CSF1 in "
     "niche myeloid cells and its downstream program in tumour cells.",
     "Their cognate receptors are expressed on the malignant cells themselves (Fig. 5e; "
     "percentages of 150,213 malignant GSE182109 cells): IL6ST (55%), LRP1 (56%), ITGB1, EGFR "
     "(38%), BMPR2 (35%), CD44, FGFR1 (27%) and TGFBR1, whereas CSF1R is nearly absent (7%) and "
     "is now plotted explicitly. CSF1 is therefore the one top-ranked ligand whose receptor the "
     "malignant cells do not carry: its signal must be read as acting through a myeloid "
     "intermediary in the niche rather than on malignant cells directly, and the same caveat "
     "applies to the other top ligands whose prior-matrix ranking may be inflated."),
    ("sub", "3ED7EBF0",
     "maps directly onto the regulatory arms identified in Section 3.3, connecting the "
     "microenvironment analysis to the switch.",
     "maps onto the regulatory arms identified above, connecting the microenvironment analysis "
     "to the regulatory core, although only the IL-6/STAT3 and TGF-β/SMAD arms are supported by "
     "both routes."),

    # ---------- Results: survival (M6) ----------------------------------------------------
    ("sub", "160D1F6E",
     "In TCGA-GBM (n = 939), the MES score was adverse-prognostic in a multivariable Cox model "
     "adjusting for age and IDH (HR 1.22, P = 0.022; Fig. 5f), and stratified Kaplan–Meier "
     "analysis showed separation also within the IDH-mutant subset (Supplementary Fig. S2b).",
     "In TCGA (809 patients with expression, 737 with survival follow-up; the cohort includes the "
     "LGG project, so it is a lower-grade-inclusive glioma series rather than TCGA-GBM), the MES "
     "score was adverse-prognostic in a multivariable Cox model adjusting for age and IDH "
     "(HR 1.15 per SD, 95% CI 1.02–1.30, P = 0.024; proportional hazards not violated, cox.zph "
     "P = 0.32; Fig. 5f, Supplementary Fig. S2a), but the association did not survive grade "
     "adjustment (HR 1.10, 0.98–1.24, P = 0.096) or adjustment for stromal and myeloid content "
     "(HR 0.92, 0.75–1.11, P = 0.38), which indicates that in bulk TCGA the score partly measures "
     "stromal content. It did survive adjustment for the deconvolved malignant fraction "
     "(HR 1.13, 1.00–1.28, P = 0.044). Kaplan–Meier curves use a median cut-off with log-rank P "
     "values and numbers at risk (Supplementary Fig. S2b)."),
    ("sub", "160D1F6E",
     "The association replicated in both external CGGA series (CGGA-325 HR 2.28; CGGA-693 HR "
     "1.62; both P &lt; 10⁻³) and remained significant after adjustment for WHO grade "
     "(grade-adjusted HR 1.79 and 1.33; Supplementary Fig. S2a). Pooled Kaplan–Meier analysis "
     "showed clear high-versus-low survival separation (Supplementary Fig. S2c), and the MES "
     "score increased with WHO grade (Supplementary Fig. S2d). A fixed-effect meta-analysis "
     "across the three cohorts gave HR 1.55 (95% CI 1.40–1.72); the random-effects estimate was "
     "HR 1.63 (1.18–2.26, P = 0.003) with high between-cohort heterogeneity (I² = 89%; Fig. 5f), "
     "as expected from the differing platforms and case mixes.",
     "The association replicated in both external CGGA series and remained significant after "
     "grade adjustment: CGGA-325 HR 1.70 (1.46–1.98), grade-adjusted 1.45 (1.23–1.72); CGGA-693 "
     "HR 1.41 (1.26–1.57), grade-adjusted 1.23 (1.10–1.37). In CGGA it also survived adjustment "
     "for stromal and myeloid content (1.44 and 1.43). The MES score rises with WHO grade "
     "(Kruskal–Wallis P = 4.4 × 10⁻³⁰; Spearman rho = 0.43; Supplementary Fig. S2d), which is why "
     "we foreground the grade-adjusted models. Pooled across the three cohorts the "
     "random-effects estimate — reported as primary because between-cohort heterogeneity is high "
     "(I² = 88%) — is HR 1.39 (1.13–1.72, P = 0.0018), with a between-cohort range of 1.15–1.70 "
     "and a fixed-effect estimate of 1.37 (1.27–1.47) shown for comparison (Supplementary Fig. "
     "S2a)."),
    ("sub", "160D1F6E",
     "The MES-transition score is thus a robust, independently validated adverse-prognosis "
     "biomarker.",
     "We therefore present the score as a reproducible but partly grade- and stroma-dependent "
     "prognostic correlate rather than as a validated biomarker."),

    # ---------- Discussion ----------------------------------------------------------------
    ("sub", "03BE3C7E",
     "It is not enriched at recurrence (Fig. 2d,e), which argues against therapy-driven "
     "selection of a MES clone and in favour of a reversible, context-driven program.",
     "We found no evidence of enrichment at recurrence (Fig. 2d,e), which is consistent with a "
     "reversible, context-driven program, although that test is power-limited and the conclusion "
     "rests mainly on the concordant external longitudinal studies [2,12,28]."),
    ("para", "0E5EA46D",
     "A compact AP-1/C-EBP core tracks the program, and part of the evidence does not depend on "
     "regulatory priors. Prior-based activity scoring recovers the AP-1 factors and CEBPD "
     "reproducibly across cohorts (Fig. 3a,d), and the prior-independent route also ranks AP-1 "
     "members among its strongest regulons (Fig. 3b), although the two routes overlap in only a "
     "handful of factors and disagree on CEBPB and on the rank of STAT3 (Fig. 3c). This is weaker "
     "than the convergence we originally claimed, and it is also less novel than we implied: "
     "C/EBPbeta and STAT3 were identified as master regulators of mesenchymal transformation in "
     "glioma more than a decade ago [26], the NF-kB arm and its link to radioresistance have "
     "precedent [21,27], and bulk analyses had already linked mesenchymal evolution to the "
     "microenvironment [28,29]. Our contribution is the multi-cohort, two-route reproducibility "
     "test, the exclusion of TWIST1, and the quantification of how much of the correlation "
     "survives patient-level inference — not the discovery of the axis. We therefore describe "
     "these factors as candidate control nodes, and note that AP-1 and C/EBP remain difficult to "
     "drug: the C/EBP-inhibitory peptide programme [24] is the one clinical-adjacent effort and "
     "its translation to glioma is unproven."),
    ("para", "1AACB179",
     "The program is connected to a defined niche and to clinical outcome. Its localization "
     "(IvyGAP; Fig. 5a–d) corrects our own earlier reading rather than the literature: we had "
     "reported perivascular enrichment, but after excluding the FACS-sorted samples that are not "
     "anatomic regions, separating the hypoxia-linked MES1 submodule from MES2, and adjusting for "
     "vascular and myeloid content, the MES programme is highest in the pseudopalisading/peri-"
     "necrotic compartment, which restores the classical hypoxic, necrosis-adjacent association "
     "for this bulk regional resource [1,25]. NicheNet nominates IL-6, FGF, EGF and TGF-β/BMP "
     "family ligands of the vascular–myeloid compartment, whose receptors the malignant cells do "
     "carry (Fig. 5e), feeding the STAT3, AP-1 and SMAD arms of the core. The MES score itself "
     "replicates in both external CGGA cohorts and, more weakly, in TCGA (Fig. 5f, Supplementary "
     "Fig. S2): it is prognostic, but it is also grade- and stroma-correlated, so we do not "
     "present it as a ready clinical marker. These results are compatible with a model in which "
     "perivascular and peri-necrotic paracrine signalling drives susceptible NPC/AC-state cells "
     "through an AP-1/C-EBP-dependent transition that is reflected in survival at the tumour "
     "level."),
    ("sub", "07505CF4",
     "One candidate dataset (GSE163120) was excluded a priori because it contains only immune "
     "cells and a CCR2 knockout, and is therefore unsuitable for testing the tumour MES program.",
     "One candidate dataset (GSE163120) was excluded a priori because it contains only immune "
     "cells and a CCR2 knockout, and is therefore unsuitable for testing the tumour MES program. "
     "Further limitations specific to the revised analyses: the de-novo arm was re-run with a "
     "smaller state-stratified subsample (2,000 cells per cohort) than our first analysis because "
     "the original code was not recoverable, so its ranking supersedes rather than reproduces the "
     "submitted one; the IvyGAP deconvolution is a signature-matrix non-negative least-squares "
     "fit, not a full probabilistic deconvolution; the bulk MES score correlates with WHO grade "
     "and with stromal content, and its prognostic value in TCGA survives neither adjustment "
     "while in CGGA it survives grade adjustment with a smaller effect size; the TCGA arm "
     "includes low-grade glioma and is therefore not a GBM-only analysis; the submitted IDH "
     "comparison and the submitted survival numbers were built from clinical fields whose "
     "event/time coding disagrees with the curated tables for about a fifth of cases, which is "
     "why every number in this version is recomputed and stated here; and per-cell statistics "
     "remain descriptive throughout, with patient-level estimates reported alongside them."),
    # ---------- M7 / C5: the release must carry the whole pipeline ------------------------
    ("sub", "7F70A6DF",
     "The de-novo pySCENIC, IvyGAP, NicheNet and CGGA analyses, the GSE89567 IDH-mutant "
     "comparison, and the figure-composition scripts are not part of this release; the "
     "corresponding code and result tables are available from the corresponding authors on "
     "request.",
     "The revision released with this manuscript adds the previously missing modules: the "
     "de-novo pySCENIC chain (loom export, GRNBoost2, cisTarget and AUCell) with the "
     "regulon-ranking and transcription-factor-anchoring audits; the IvyGAP region analysis with "
     "its composition adjustment and signature-matrix deconvolution; the base-R NicheNet ligand "
     "activity and receptor-expression analysis; the CGGA-325/693 scoring and survival models; "
     "the corrected GSE89567 and GSE70630 IDH-mutant comparison; the revision figure scripts; "
     "and every result table quoted in this manuscript, including the gene lists (Neftel "
     "MES1+MES2, the 128-gene MES-transition programme), the full ranked regulon tables, the "
     "per-cohort cohort-flow and model tables, and the single-cell objects needed to re-derive "
     "them. The archive shipped with the submission (Paper_2_代码.zip) contains the same code and "
     "tables, so no result in this manuscript depends on a file that is available only on "
     "request."),
    ("sub", "07505CF4",
     "The clinical cohorts differ in platform and scale, so we replicate the direction and "
     "significance of the prognostic effect rather than its magnitude (I² = 89%).",
     "The clinical cohorts differ in platform and scale, so we replicate the direction rather "
     "than the magnitude of the prognostic effect (I² = 88%)."),
    ("sub", "4DF903D0",
     " Ligand activity against the 73-gene MES-transition program was computed",
     " Ligand activity against the 128-gene MES-transition programme was computed"),
    ("sub", "50242B74",
     ' (a) IvyGAP MES score (z per tumour) across anatomic regions (LE leading edge, IT infiltrating tumour, CT cellular tumour, MVP microvascular proliferation, HBV hyperplastic blood vessels, PNZ peri-necrotic zone, PAN pseudopalisading around necrosis); Kruskal–Wallis P = 8.7 × 10⁻²⁰; MVP and HBV significantly above LE (FDR ≈ 1 × 10⁻⁷ and 1.8 × 10⁻⁵), PNZ/PAN not significant. (b) Hypoxia score across the same regions; unlike MES, hypoxia also trends up near necrosis, showing the MES signal is not a hypoxia proxy. (c) Mixed-effects model region coefficients (reference LE, random intercept per tumour): MVP +0.28, HBV +0.27. (d) NicheNet ligand activity against the 73-gene MES-transition program (top: CSF1, TFPI, ITM2B, FGF5, BMP3, IL6, BMP2, EPGN, ANGPT4). (e) Percentage of malignant GSE182109 cells expressing each receptor (IL6ST 55%, LRP1 56%, EGFR 38%, BMPR2 35%, FGFR1 27%; CSF1R 7%). (f) Forest plot of the MES-score hazard ratio for overall survival (multivariable Cox, MES plus age and IDH): TCGA HR 1.22, CGGA-325 HR 2.28, CGGA-693 HR 1.62, CGGA pooled, and random-effects meta-analysis HR 1.63 (95% CI 1.18–2.26, P = 0.003; I² = 89%).',
     ' (a) IvyGAP anatomic regions: MES (MES1+MES2), MES1 and MES2 scores (z-normalised per tumour) in the 122 laser-microdissected samples; points are samples, bars are medians; n samples/tumours: LE 19/8, IT 24/8, CT 30/10, MVP 25/9, PAN 24/8. The FACS-sorted cancer-stem-cell-cluster samples (CThbv-/CTpnz-) are excluded because they are not anatomic regions. (b) Hypoxia score of the same regions (the 15-gene HIF-target set of the Methods). (c) Mixed-effects model region coefficients (reference LE, random intercept per tumour), unadjusted and adjusted for vascular (endothelial+pericyte) and myeloid marker scores; the MVP effect is halved by adjustment while the PAN effect persists. (d) Signature-matrix (NNLS) deconvolution against the GSE103224 single-cell reference: malignant fraction (black) and MES share of malignant cells (red) per region. (e) NicheNet ligand activity against the 128-gene MES-transition programme (top 15). (f) Cox hazard ratios per SD of the MES score, unadjusted and grade-adjusted, in TCGA and CGGA-325/693; Supplementary Fig. S2a gives the pooled random-effects estimate.'),
]
