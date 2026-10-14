#!/usr/bin/env python3
"""rev04 — IvyGAP niche analysis, re-run and composition-controlled.

Reviewer items M3 / DA4 / minor 7 / minor 22-23: the submitted perivascular-niche claim
was (i) confounded by the cellular composition of bulk, laser-microdissected IvyGAP
regions, (ii) reported without sample or tumour counts, (iii) pooling the hypoxia-linked
MES1 submodule with MES2, and (iv) produced by code that no longer exists.

What this script does, in order:

  0. Loads GSE107559 (IvyGAP) and documents the two *different* kinds of sample in it:
     - `anatomic structures`  = 122 laser-microdissected region samples, 5 structures
       (CT 30, CTmvp 25, CTpan 24, IT 24, LE 19) from 10 tumours  <-- the niche analysis
     - `cancer stem cell clusters` = 148 FACS-sorted CSC-cluster samples, which are NOT
       anatomic regions; they carry CThbv-*/CTpnz-* labels.
  1. Scores the Neftel MES1 / MES2 / pooled MES modules, a documented 15-gene HIF-target
     hypoxia set, the 73-gene MES-transition program, and per-compartment marker scores
     (endothelial, pericyte, myeloid, T cell, oligodendrocyte, astrocyte).
  2. z-normalises per tumour (as the submitted Methods state) and tests region effects
     (Kruskal-Wallis + region-vs-LE Mann-Whitney with BH FDR), reporting n per region.
  3. Splits MES1 (hypoxia-linked) from MES2.
  4. Controls the region effect for composition: correlation with compartment markers,
     OLS adjustment, a linear mixed model with a random intercept per tumour, and a
     signature-matrix (NNLS) deconvolution against an in-hand scRNA reference, testing
     the *malignant MES* fraction per region.

Outputs (results/revision/):
    rev04_ivygap_samples.tsv        per-sample scores + region/tumour
    rev04_ivygap_region_stats.tsv   KW + region-vs-LE MWU/FDR, anatomic and pooled sets
    rev04_ivygap_mes12.tsv          MES1 vs MES2 region means + tests
    rev04_ivygap_composition.tsv    MES vs compartment-marker correlations
    rev04_ivygap_adjusted.tsv       OLS- and mixed-model region effects, unadjusted vs adjusted
    rev04_ivygap_deconv.tsv         NNLS cell-type fractions + malignant-MES fraction per sample
    rev04_ivygap_deconv_stats.tsv   region tests on the deconvolved malignant-MES fraction

Run:
  GLIOMA_PROJ=... GLIOMA_DATA=... \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev04_ivygap_composition.py
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from scipy import stats
from scipy.optimize import nnls
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "revision")
os.makedirs(OUT, exist_ok=True)
IVY = os.path.join(lib.DATASETS, "IvyGAP")
SEED = 42

REGIONS = ["LE", "IT", "CT", "CTmvp", "CTpan"]
REGIONS_POOLED = ["LE", "IT", "CT", "CTmvp", "CTpan", "CThbv", "CTpnz"]
REGION_LABEL = {"LE": "Leading edge", "IT": "Infiltrating tumour", "CT": "Cellular tumour",
                "CTmvp": "Microvascular proliferation", "CTpan": "Pseudopalisading around necrosis",
                "CThbv": "CSC-cluster: hyperplastic blood vessels*", "CTpnz": "CSC-cluster: peri-necrotic*"}

# 15-gene HIF-1alpha target set used here (reported verbatim in the revised Methods so the
# reader can reproduce the hypoxia control; the submitted manuscript did not name its set).
HYPOXIA = ["VEGFA", "SLC2A1", "LDHA", "PGK1", "PDK1", "HK2", "BNIP3", "NDRG1", "ENO1",
           "PDK3", "TPI1", "P4HA1", "ADM", "MRPS17", "ANKZF1"]
COMPARTMENTS = {
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "ENG", "FLT1", "KDR", "ESAM"],
    "pericyte":    ["RGS5", "PDGFRB", "ACTA2", "NOTCH3", "CSPG4", "MCAM", "TAGLN"],
    "myeloid":     ["CD68", "AIF1", "C1QA", "C1QB", "C1QC", "CSF1R", "ITGAM", "PTPRC", "LYZ", "CD14"],
    "tcell":       ["CD3D", "CD3E", "CD8A", "CD2", "IL7R"],
    "oligodendrocyte": ["PLP1", "MBP", "MOG", "MAG", "CLDN11"],
    "astrocyte":   ["GFAP", "AQP4", "S100B", "ALDH1L1"],
}


# ------------------------------------------------------------------ loading
def load_ivygap():
    expr = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_fpkm_table.csv.gz"),
                       index_col=0, compression="gzip")
    rg = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_rows-genes.csv.gz"), compression="gzip")
    id2s = dict(zip(rg["gene_id"].astype(str), rg["gene_symbol"].astype(str)))
    expr.index = [id2s.get(str(g), str(g)) for g in expr.index]
    expr = expr[expr.index.notna()]
    expr = expr[~expr.index.duplicated(keep="first")]
    log = np.log2(expr + 1)
    ann = pd.read_excel(os.path.join(IVY, "GSE107559_ivygap_columns-samples.xlsx"))
    ann["rna_well_id"] = ann["rna_well_id"].astype(str)
    ann = ann.set_index("rna_well_id")
    lib.log.info(f"IvyGAP: {log.shape[0]} genes x {log.shape[1]} samples; "
                 f"annotation rows {ann.shape[0]}")
    return log, ann


def region_code(s):
    s = str(s)
    for c in ["CTmvp", "CTpan", "CThbv", "CTpnz", "CT", "IT", "LE"]:
        if s.startswith(c):
            return c
    return np.nan


def mean_z(expr, genes):
    g = [x for x in genes if x in expr.index]
    if len(g) < 3:
        return pd.Series(np.nan, index=expr.columns)
    s = expr.loc[g]
    z = s.sub(s.mean(axis=1), axis=0).div(s.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def score_samples(log, ann):
    mods = load_neftel_modules()
    prog = pd.read_csv(os.path.join(lib.PROJ, "results", "step04_traj",
                                    "mes_transition_program_core.tsv"), sep="\t")
    prog_up = list(prog.loc[prog["rho"] > 0, "gene"])
    s = pd.DataFrame(index=log.columns)
    s["MES"] = mean_z(log, mods.get("MES1", []) + mods.get("MES2", []))
    s["MES1"] = mean_z(log, mods.get("MES1", []))
    s["MES2"] = mean_z(log, mods.get("MES2", []))
    s["PROG73"] = mean_z(log, prog_up)                 # the MES-transition program
    s["hypoxia"] = mean_z(log, HYPOXIA)
    for name, genes in COMPARTMENTS.items():
        s[f"c_{name}"] = mean_z(log, genes)
    s["structure"] = [ann["structure_abbreviation"].get(c, np.nan) for c in s.index]
    s["tumour"] = [ann["tumor_name"].get(c, np.nan) for c in s.index]
    s["study"] = [ann["study"].get(c, np.nan) for c in s.index]
    s["region"] = s["structure"].map(region_code)
    s = s[s["region"].notna() & s["tumour"].notna()].copy()
    # per-tumour z-normalisation of the score columns (as in the submitted Methods)
    for col in ["MES", "MES1", "MES2", "PROG73", "hypoxia"] + [f"c_{k}" for k in COMPARTMENTS]:
        s[col + "_z"] = s.groupby("tumour")[col].transform(
            lambda v: (v - v.mean()) / (v.std(ddof=0) if v.std(ddof=0) > 0 else np.nan))
    s["vascular_z"] = s[["c_endothelial_z", "c_pericyte_z"]].mean(axis=1)
    s["myeloid_z"] = s["c_myeloid_z"]
    return s


# ------------------------------------------------------------------ tests
def region_tests(df, col, regions):
    """Kruskal-Wallis across regions + each region vs the leading edge (BH-FDR)."""
    present = [r for r in regions if (df["region"] == r).sum() >= 3]
    groups = [df.loc[df["region"] == r, col].dropna().values for r in present]
    H, p = stats.kruskal(*groups) if len(groups) > 1 else (np.nan, np.nan)
    le = df.loc[df["region"] == "LE", col].dropna().values
    others = [r for r in present if r != "LE"]
    raw = {}
    for r in others:
        v = df.loc[df["region"] == r, col].dropna().values
        raw[r] = (stats.mannwhitneyu(v, le, alternative="two-sided").pvalue
                  if len(v) >= 3 and len(le) >= 3 else np.nan)
    fdr = {}
    if raw:
        keys = [k for k, v in raw.items() if not np.isnan(v)]
        if keys:
            adj = multipletests([raw[k] for k in keys], method="fdr_bh")[1]
            fdr = dict(zip(keys, adj))
    rows = []
    for r in present:
        v = df.loc[df["region"] == r, col].dropna()
        rows.append({"score": col, "region": r,
                     "n_samples": int(len(v)),
                     "n_tumours": int(df.loc[df["region"] == r, "tumour"].nunique()),
                     "mean_z": round(float(v.mean()), 3),
                     "sd": round(float(v.std(ddof=1)), 3),
                     "raw_MWU_p": raw.get(r, np.nan),
                     "MWU_vs_LE_fdr": fdr.get(r, np.nan),
                     "kruskal_H": round(float(H), 2) if not np.isnan(H) else np.nan,
                     "kruskal_p": float(p) if not np.isnan(p) else np.nan})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ deconvolution
def build_reference():
    """Cell-type signature matrix from the in-hand GSE103224 scRNA objects."""
    mal = sc.read_h5ad(os.path.join(lib.PROC_DIR, "gse103224_regulon.h5ad"))
    refs = {}
    mes = mal.obs["score_MES"]
    hi, lo = mes.quantile(0.9), mes.quantile(0.1)
    refs["malignant_MES"] = np.asarray(mal[mes >= hi].X.mean(axis=0)).ravel()
    refs["malignant_nonMES"] = np.asarray(mal[mes <= lo].X.mean(axis=0)).ravel()
    refs["malignant_all"] = np.asarray(mal.X.mean(axis=0)).ravel()
    genes = mal.var_names.to_numpy()
    # non-malignant compartments from the QC object (log-normalised here)
    q = sc.read_h5ad(os.path.join(lib.PROC_DIR, "gse103224_qc.h5ad"))
    if not sp.issparse(q.X):
        q.X = sp.csr_matrix(q.X.astype(np.float32))
    sc.pp.filter_genes(q, min_cells=3)
    sc.pp.normalize_total(q, target_sum=1e4)
    sc.pp.log1p(q)
    q.var_names_make_unique()
    for name, geneset in [("immune", COMPARTMENTS["myeloid"] + COMPARTMENTS["tcell"]),
                          ("endothelial", COMPARTMENTS["endothelial"]),
                          ("pericyte", COMPARTMENTS["pericyte"]),
                          ("oligodendrocyte", COMPARTMENTS["oligodendrocyte"])]:
        present = [g for g in geneset if g in q.var_names]
        sc.tl.score_genes(q, present, score_name=f"sc_{name}", ctrl_size=50)
        sel = q.obs[f"sc_{name}"] > q.obs[f"sc_{name}"].quantile(0.95)
        refs[name] = np.asarray(q[sel].X.mean(axis=0)).ravel() if sel.sum() > 20 else None
    ref_genes = q.var_names.to_numpy()
    lib.log.info(f"reference: {len(refs)} cell types; malignant genes {len(genes)}, qc genes {len(ref_genes)}")
    return refs, genes, ref_genes


def deconvolve(log, refs, genes, ref_genes):
    common = [g for g in genes if g in ref_genes and g in log.index]
    ref = {}
    for k, v in refs.items():
        if v is None:
            continue
        ref[k] = pd.Series(v, index=(genes if len(v) == len(genes) else ref_genes)).reindex(common)
    R = pd.DataFrame(ref).fillna(0.0)
    M = log.loc[common]
    # centre both matrices (bulk and reference) so NNLS fits expression *shapes*
    Rc = R.sub(R.mean(axis=0), axis=1)
    Mc = M.sub(M.mean(axis=0), axis=1)
    A = Rc.to_numpy()
    frac = np.zeros((M.shape[1], A.shape[1]))
    for j in range(M.shape[1]):
        y = Mc.iloc[:, j].to_numpy()
        y = np.nan_to_num(y)
        w, _ = nnls(A, y)
        frac[j] = w
    frac = pd.DataFrame(frac, index=M.columns, columns=R.columns)
    tot = frac.sum(axis=1).replace(0, np.nan)
    frac_n = frac.div(tot, axis=0)
    out = frac_n.copy()
    out["malignant_total"] = frac_n[[c for c in frac_n.columns if c.startswith("malignant")]].sum(axis=1)
    out["MES_share_of_malignant"] = frac_n["malignant_MES"] / out["malignant_total"].replace(0, np.nan)
    return out


# ------------------------------------------------------------------ main
def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev04_ivygap.log"))
    log, ann = load_ivygap()
    s = score_samples(log, ann)
    s.to_csv(os.path.join(OUT, "rev04_ivygap_samples.tsv"), sep="\t")

    anatomic = s[s["study"] == "anatomic structures"].copy()
    pooled = s.copy()
    lib.log.info(f"anatomic-structure samples: {len(anatomic)} from {anatomic['tumour'].nunique()} tumours; "
                 f"pooled set (incl. FACS-sorted CSC clusters): {len(pooled)}")

    # ---- region statistics, both sample sets (A = correct, B = what the P-value implies) ----
    stats_rows = []
    for col in ["MES_z", "MES1_z", "MES2_z", "PROG73_z", "hypoxia_z"]:
        for tag, df, regs in [("anatomic", anatomic, REGIONS), ("pooled_incl_CSC", pooled, REGIONS_POOLED)]:
            if (df["region"] == "LE").sum() < 3:
                continue
            t = region_tests(df, col, regs)
            t["sample_set"] = tag
            stats_rows.append(t)
    st = pd.concat(stats_rows, ignore_index=True)
    st.to_csv(os.path.join(OUT, "rev04_ivygap_region_stats.tsv"), sep="\t", index=False)
    lib.log.info("anatomic-set region means for MES_z:\n" +
                 st[(st["sample_set"] == "anatomic") & (st["score"] == "MES_z")]
                 [["region", "n_samples", "n_tumours", "mean_z", "raw_MWU_p", "MWU_vs_LE_fdr"]]
                 .to_string(index=False))
    lib.log.info("pooled-set region means for MES_z:\n" +
                 st[(st["sample_set"] == "pooled_incl_CSC") & (st["score"] == "MES_z")]
                 [["region", "n_samples", "n_tumours", "mean_z", "raw_MWU_p", "MWU_vs_LE_fdr"]]
                 .to_string(index=False))

    # ---- MES1 vs MES2 ----
    m12 = st[st["score"].isin(["MES1_z", "MES2_z"])][
        ["sample_set", "score", "region", "n_samples", "mean_z", "MWU_vs_LE_fdr", "kruskal_p"]]
    m12.to_csv(os.path.join(OUT, "rev04_ivygap_mes12.tsv"), sep="\t", index=False)
    lib.log.info("MES1 (hypoxia-linked) vs MES2 (mesenchymal), anatomic set:\n" +
                 m12[m12["sample_set"] == "anatomic"].to_string(index=False))

    # ---- composition: correlation of the MES score with compartment content ----
    crow = []
    for k in COMPARTMENTS:
        rho, p = stats.spearmanr(anatomic["MES_z"], anatomic[f"c_{k}_z"], nan_policy="omit")
        rh2, p2 = stats.spearmanr(anatomic["MES1_z"], anatomic[f"c_{k}_z"], nan_policy="omit")
        rh3, p3 = stats.spearmanr(anatomic["MES2_z"], anatomic[f"c_{k}_z"], nan_policy="omit")
        crow.append({"dataset": "anatomic", "compartment": k,
                     "rho_MES": round(rho, 3), "p_MES": p,
                     "rho_MES1": round(rh2, 3), "p_MES1": p2,
                     "rho_MES2": round(rh3, 3), "p_MES2": p3})
    comp_df = pd.DataFrame(crow)
    comp_df.to_csv(os.path.join(OUT, "rev04_ivygap_composition.tsv"), sep="\t", index=False)
    lib.log.info("MES vs compartment content (anatomic set):\n" + comp_df.to_string(index=False))

    # ---- adjusted region effects: OLS + mixed model with tumour random intercept ----
    d = anatomic.dropna(subset=["MES_z", "vascular_z", "myeloid_z"]).copy()
    d["region"] = pd.Categorical(d["region"], categories=REGIONS)
    rows = []
    for label, formula in [("OLS_unadjusted", "MES_z ~ C(region)"),
                           ("OLS_adj_vascular_myeloid", "MES_z ~ C(region) + vascular_z + myeloid_z"),
                           ("OLS_adj_vascular_only", "MES_z ~ C(region) + vascular_z")]:
        m = smf.ols(formula, data=d).fit()
        for term, coef in m.params.items():
            if "region" not in term:
                continue
            ci = m.conf_int().loc[term]
            rows.append({"model": label, "term": term, "coef": round(float(coef), 3),
                         "ci_lo": round(float(ci[0]), 3), "ci_hi": round(float(ci[1]), 3),
                         "p": float(m.pvalues[term]), "n": int(m.nobs)})
        lib.log.info(f"{label}: n={int(m.nobs)}, R2={m.rsquared:.3f}")
    try:
        mm = smf.mixedlm("MES_z ~ C(region)", d, groups=d["tumour"]).fit(reml=False)
        for term, coef in mm.params.items():
            if "region" not in term:
                continue
            rows.append({"model": "MixedLM_unadjusted", "term": term, "coef": round(float(coef), 3),
                         "ci_lo": np.nan, "ci_hi": np.nan, "p": float(mm.pvalues[term]), "n": int(mm.nobs)})
        lib.log.info(f"MixedLM (random intercept per tumour): {mm.params.to_dict()}")
        mm2 = smf.mixedlm("MES_z ~ C(region) + vascular_z + myeloid_z", d, groups=d["tumour"]).fit(reml=False)
        for term, coef in mm2.params.items():
            if "region" not in term:
                continue
            rows.append({"model": "MixedLM_adj_vascular_myeloid", "term": term, "coef": round(float(coef), 3),
                         "ci_lo": np.nan, "ci_hi": np.nan, "p": float(mm2.pvalues[term]), "n": int(mm2.nobs)})
        lib.log.info(f"MixedLM adjusted: {mm2.params.to_dict()}")
    except Exception as e:
        lib.log.warning(f"mixed model failed: {e}")
    adj = pd.DataFrame(rows)
    adj.to_csv(os.path.join(OUT, "rev04_ivygap_adjusted.tsv"), sep="\t", index=False)

    # ---- deconvolution ----
    try:
        refs, genes, ref_genes = build_reference()
        dec = deconvolve(log, refs, genes, ref_genes)
        dec.to_csv(os.path.join(OUT, "rev04_ivygap_deconv.tsv"), sep="\t")
        dd = anatomic.join(dec, how="inner")
        for col in ["MES_share_of_malignant", "malignant_total", "endothelial", "pericyte", "myeloid", "immune"]:
            if col not in dd.columns:
                continue
            t = region_tests(dd, col, REGIONS)
            t["sample_set"] = "anatomic"
            st2 = pd.concat([st, t], ignore_index=True)
            st = st2
        st.to_csv(os.path.join(OUT, "rev04_ivygap_region_stats.tsv"), sep="\t", index=False)
        mask = ((st["sample_set"] == "anatomic") &
                st["score"].isin(["MES_share_of_malignant", "malignant_total", "endothelial",
                                  "pericyte", "myeloid", "immune"]))
        dstat = st.loc[mask, ["score", "region", "n_samples", "mean_z",
                              "raw_MWU_p", "MWU_vs_LE_fdr"]]
        dstat.to_csv(os.path.join(OUT, "rev04_ivygap_deconv_stats.tsv"), sep="\t", index=False)
        lib.log.info("deconvolution-based region test (malignant MES share + compartments):\n" +
                     dstat.to_string(index=False))
    except Exception as e:
        lib.log.warning(f"deconvolution failed: {e}")
    lib.log.info("DONE rev04")


if __name__ == "__main__":
    main()
