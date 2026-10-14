#!/usr/bin/env python3
"""rev05 — Clinical-cohort re-analysis inputs (TCGA + CGGA-325/693).

Reviewer items M6 / DA8 / minor 12-13 / 2.3-2.6: the submitted survival analysis
reported "TCGA, n = 939" without a sample-flow, without stating its event/time
construction, without HR units or CIs, without proportional-hazards diagnostics,
without Kaplan-Meier cut-off rules / log-rank P / at-risk tables, without any purity
adjustment, and without a benchmark against the published Neftel MES module — and its
CGGA arm used a "73-gene program" that no surviving artefact defines.

This script builds the case-level inputs and resolves the gene-set question:

  * `MES_Neftel`  = z-mean of the Neftel MES1+MES2 meta-module (100 genes) — this is the
    score the surviving bulk pipeline (scripts/tcga_mes_score.py) actually computed.
  * `MES_prog131` = z-mean of the 131-gene core MES-transition program
    (results/step04_traj/mes_transition_program_core.tsv; 128 positive-rho genes),
    reported as a benchmark.
  * `stromal_content` = z-mean of endothelial + pericyte marker genes, and
    `myeloid_content`  = z-mean of myeloid markers, so the prognostic claim can be
    tested against the "the score is just stromal content" alternative.
  * `malignant_frac`, `MES_share_of_malignant` = signature-matrix (NNLS) deconvolution
    against the in-hand GSE103224 scRNA reference, as a bulk-purity surrogate
    (the ESTIMATE package is not available on this machine).

Outputs (results/revision/):
    rev05_tcga_cases.tsv, rev05_cgga325_cases.tsv, rev05_cgga693_cases.tsv
    rev05_cohort_flow.tsv          sample -> case -> survival funnel per cohort
    rev05_reference_matrix.tsv     cached deconvolution reference (cell types x genes)
    rev05_genesets.tsv             the exact gene lists used (MES_Neftel, MES_prog131)

Run:
  GLIOMA_PROJ=... GLIOMA_DATA=... \
  /run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python scripts/rev05_survival_prep.py
"""
import os
import sys
import glob
import zipfile

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
import lib
from genesets import load_neftel_modules
import rev04_ivygap_composition as R4

sc.settings.verbosity = 0
OUT = os.path.join(lib.PROJ, "results", "revision")
os.makedirs(OUT, exist_ok=True)
TCGA_DIR = os.path.join(lib.RAW_DIR, "tcga")
CGGA_DIR = os.path.join(lib.DATASETS, "CGGA")
FLOW = []


def genesets():
    mods = load_neftel_modules()
    mes = sorted(set(mods["MES1"] + mods["MES2"]))
    prog = pd.read_csv(os.path.join(lib.PROJ, "results", "step04_traj",
                                    "mes_transition_program_core.tsv"), sep="\t")
    prog_up = sorted(set(prog.loc[prog["rho"] > 0, "gene"]))
    gs = pd.DataFrame({"MES_Neftel": pd.Series(mes), "MES_prog131": pd.Series(prog_up)})
    gs.to_csv(os.path.join(OUT, "rev05_genesets.tsv"), sep="\t", index=False)
    return mes, prog_up


def reference():
    """Cached cell-type reference (mean log-normalised expression per compartment)."""
    cache = os.path.join(OUT, "rev05_reference_matrix.tsv")
    if os.path.exists(cache):
        r = pd.read_csv(cache, sep="\t", index_col=0)
        lib.log.info(f"reference cache: {r.shape}")
        return r
    refs, genes, ref_genes = R4.build_reference()
    cols = {}
    for k, v in refs.items():
        if v is None:
            continue
        cols[k] = pd.Series(v, index=(genes if len(v) == len(genes) else ref_genes))
    r = pd.DataFrame(cols)
    r.index.name = "gene"
    r.to_csv(cache, sep="\t")
    lib.log.info(f"reference built and cached: {r.shape}")
    return r


def nnls_fractions(logmat, ref):
    """logmat: genes x samples (log2 scale). Returns fraction table per sample."""
    common = [g for g in ref.index if g in logmat.index]
    A = ref.loc[common].to_numpy()
    A = A - A.mean(axis=0)                      # centre each signature
    M = logmat.loc[common]
    M = M.sub(M.mean(axis=0), axis=1)
    out = np.zeros((M.shape[1], A.shape[1]))
    for j in range(M.shape[1]):
        y = np.nan_to_num(M.iloc[:, j].to_numpy())
        w, _ = R4.nnls(A, y)
        out[j] = w
    f = pd.DataFrame(out, index=M.columns, columns=ref.columns)
    tot = f.sum(axis=1).replace(0, np.nan)
    f = f.div(tot, axis=0)
    mal = [c for c in f.columns if c.startswith("malignant")]
    f["malignant_frac"] = f[mal].sum(axis=1)
    f["MES_share_of_malignant"] = f["malignant_MES"] / f["malignant_frac"].replace(0, np.nan)
    return f


def zmean(expr_genes_x_samples, genes):
    g = [x for x in genes if x in expr_genes_x_samples.index]
    s = expr_genes_x_samples.loc[g]
    z = s.sub(s.mean(axis=1), axis=0).div(s.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), len(g)


# ---------------------------------------------------------------- TCGA
def tcga(mes, prog_up, ref):
    man = pd.read_csv(os.path.join(TCGA_DIR, "manifest.tsv"), sep="\t")
    files = sorted(glob.glob(os.path.join(TCGA_DIR, "*rna_seq.augmented_star_gene_counts.tsv")))
    FLOW.append({"cohort": "TCGA", "step": "STAR-count files on disk", "n": len(files)})
    fn2 = dict(zip(man["file_name"], zip(man["sample"], man["project"])))

    need = set(mes) | set(prog_up) | set(R4.COMPARTMENTS["endothelial"]) | \
        set(R4.COMPARTMENTS["pericyte"]) | set(R4.COMPARTMENTS["myeloid"]) | set(ref.index)
    rows, samples, projects = {}, [], {}
    for f in files:
        fname = os.path.basename(f)
        meta = fn2.get(fname)
        if meta is None:
            continue
        sample, proj = meta
        df = pd.read_csv(f, sep="\t", comment="#", usecols=["gene_name", "tpm_unstranded"])
        df = df[df["gene_name"].isin(need)].drop_duplicates("gene_name")
        rows[sample] = dict(zip(df["gene_name"], np.log1p(df["tpm_unstranded"])))
        samples.append(sample)
        projects[sample] = proj
    FLOW.append({"cohort": "TCGA", "step": "files mapped to a GDC sample", "n": len(samples)})
    exp = pd.DataFrame(rows).T                       # samples x genes
    exp.index.name = "sample"
    lib.log.info(f"TCGA expression: {exp.shape}")

    cases = pd.Series([s[:12] for s in exp.index], index=exp.index)
    case_level = exp.groupby(cases).mean()
    lib.log.info(f"TCGA cases with expression: {case_level.shape[0]}")
    FLOW.append({"cohort": "TCGA", "step": "unique patients (barcode prefix 12) with expression",
                 "n": int(case_level.shape[0])})

    T = case_level.T
    sc_cols = {}
    sc_cols["MES_Neftel"], n1 = zmean(T, mes)
    sc_cols["MES_prog131"], n2 = zmean(T, prog_up)
    for k in ["endothelial", "pericyte", "myeloid"]:
        sc_cols[f"{k}_content"], _ = zmean(T, R4.COMPARTMENTS[k])
    S = pd.DataFrame(sc_cols)
    S["stromal_content"] = S[["endothelial_content", "pericyte_content"]].mean(axis=1)
    f = nnls_fractions(np.log2(np.expm1(T) + 1) if False else T, ref)
    S = S.join(f[["malignant_frac", "MES_share_of_malignant", "endothelial", "pericyte"]])
    S["project"] = [projects[s] for s in case_level.index]
    S.index.name = "case"

    cl = pd.read_csv(os.path.join(lib.DATASETS, "TCGA", "TCGA_lgggbm_clinical.tsv"), sep="\t")
    cl["case"] = cl["case"].astype(str)
    S = S.join(cl.set_index("case"), how="left")
    FLOW.append({"cohort": "TCGA", "step": "cases merged with curated clinical (os_months/os_status)",
                 "n": int(S["os_months"].notna().sum())})
    # event/time construction, stated explicitly:
    #   event = 1 if os_status contains DECEASED (the GDC `demographic.vital_status` == Dead),
    #   time  = os_months x 30.4375 days.  Cross-checked below against the survival fields that
    #           survive from the original XML route (days_to_death for deaths, otherwise
    #           days_to_last_followup), which is the construction the reviewers asked to see.
    S["event"] = S["os_status"].astype(str).str.contains("DECEASED").astype(int)
    S["os_months"] = pd.to_numeric(S["os_months"], errors="coerce")
    S["os_time_days"] = S["os_months"] * 30.4375
    xmlp = os.path.join(lib.PROJ, "results", "step07_clinical", "tcga_mes_clinical.tsv")
    if os.path.exists(xmlp):
        x = pd.read_csv(xmlp, sep="\t")
        x["case"] = x["sample"].astype(str)
        x = x.groupby("case").agg(xml_time=("os_time", "max"), xml_event=("os_event", "max"))
        S = S.join(x, how="left")
        both = S.dropna(subset=["xml_time"])
        agree_t = (np.abs(both["os_time_days"] - both["xml_time"]) <= 35).mean()
        agree_e = (both["event"] == both["xml_event"]).mean()
        lib.log.info(f"TCGA event/time cross-check vs XML route: n={len(both)}, "
                     f"event agreement={agree_e:.3f}, time within 35 d={agree_t:.3f}")
        FLOW.append({"cohort": "TCGA",
                     "step": f"XML-route cross-check: event agreement {agree_e:.3f}, "
                             f"time agreement(<=35d) {agree_t:.3f} on n={len(both)}", "n": len(both)})
    S["grade_num"] = S["grade"].astype(str).str.extract(r"G(\d)").astype(float)
    S["idh_mut"] = S["idh"].astype(str).str.upper().str.startswith("MUT").astype(float)
    S["n_genes_MES_Neftel"] = n1
    S["n_genes_MES_prog131"] = n2
    S = S.reset_index()
    S.to_csv(os.path.join(OUT, "rev05_tcga_cases.tsv"), sep="\t", index=False)
    lib.log.info(f"TCGA case table: {S.shape}; events={int(S['event'].sum())}; "
                 f"grade known={int(S['grade_num'].notna().sum())}; "
                 f"projects={S['project'].value_counts().to_dict()}")
    return S


# ---------------------------------------------------------------- CGGA
def cgga(tag, mes, prog_up, ref):
    zpath = os.path.join(CGGA_DIR, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip")
    cpath = os.path.join(CGGA_DIR, f"CGGA.mRNAseq_{tag}_clinical.zip")
    tmp = os.path.join("/tmp", f"cgga_{tag}")
    os.makedirs(tmp, exist_ok=True)
    for z in (zpath, cpath):
        zipfile.ZipFile(z).extractall(tmp)
    expf = glob.glob(os.path.join(tmp, f"CGGA.mRNAseq_{tag}.RSEM-genes*.txt"))[0]
    clinf = glob.glob(os.path.join(tmp, f"CGGA.mRNAseq_{tag}_clinical*.txt"))[0]
    exp = pd.read_csv(expf, sep="\t", index_col=0)
    exp = exp[~exp.index.duplicated(keep="first")]
    cl = pd.read_csv(clinf, sep="\t")
    FLOW.append({"cohort": f"CGGA-{tag}", "step": "samples with RSEM expression", "n": int(exp.shape[1])})
    FLOW.append({"cohort": f"CGGA-{tag}", "step": "clinical records", "n": int(cl.shape[0])})
    cl = cl.rename(columns={"Censor (alive=0; dead=1)": "event", "OS": "os_days",
                            "IDH_mutation_status": "idh", "MGMTp_methylation_status": "mgmt",
                            "PRS_type": "prs"})
    cl["os_days"] = pd.to_numeric(cl["os_days"], errors="coerce")
    cl["event"] = pd.to_numeric(cl["event"], errors="coerce")
    cl["age"] = pd.to_numeric(cl["Age"], errors="coerce")
    # CGGA encodes WHO grade in Roman numerals ("WHO II"/"WHO III"/"WHO IV")
    rom = {"II": 2, "III": 3, "IV": 4}
    cl["grade_num"] = cl["Grade"].astype(str).str.replace("^WHO ", "", regex=True).map(rom)
    cl["idh_mut"] = cl["idh"].astype(str).str.lower().str.startswith("mut").astype(float)
    keep = [c for c in exp.columns if c in set(cl["CGGA_ID"])]
    FLOW.append({"cohort": f"CGGA-{tag}", "step": "samples with both expression and clinical",
                 "n": len(keep)})
    E = exp[keep]
    log = np.log2(E + 1)
    sc_cols = {}
    sc_cols["MES_Neftel"], n1 = zmean(log, mes)
    sc_cols["MES_prog131"], n2 = zmean(log, prog_up)
    for k in ["endothelial", "pericyte", "myeloid"]:
        sc_cols[f"{k}_content"], _ = zmean(log, R4.COMPARTMENTS[k])
    S = pd.DataFrame(sc_cols)
    S["stromal_content"] = S[["endothelial_content", "pericyte_content"]].mean(axis=1)
    f = nnls_fractions(log, ref)
    S = S.join(f[["malignant_frac", "MES_share_of_malignant"]])
    S.index.name = "CGGA_ID"
    S = S.reset_index().merge(cl, on="CGGA_ID", how="inner")
    S["n_genes_MES_Neftel"] = n1
    S["n_genes_MES_prog131"] = n2
    S.to_csv(os.path.join(OUT, f"rev05_cgga{tag}_cases.tsv"), sep="\t", index=False)
    lib.log.info(f"CGGA-{tag}: {S.shape}; events={int(pd.to_numeric(S['event'], errors='coerce').sum())}; "
                 f"grade known={int(S['grade_num'].notna().sum())}")
    return S


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "rev05_survival.log"))
    mes, prog_up = genesets()
    lib.log.info(f"gene sets: MES_Neftel={len(mes)}, MES_prog131={len(prog_up)}")
    ref = reference()
    t = tcga(mes, prog_up, ref)
    c325 = cgga("325", mes, prog_up, ref)
    c693 = cgga("693", mes, prog_up, ref)
    for nm, d in [("TCGA", t), ("CGGA-325", c325), ("CGGA-693", c693)]:
        ev = pd.to_numeric(d["event"], errors="coerce")
        tm = pd.to_numeric(d.get("os_time_days", d.get("os_days")), errors="coerce")
        FLOW.append({"cohort": nm, "step": "final analysis set (event + time available)",
                     "n": int((ev.notna() & tm.notna() & (tm > 0)).sum())})
    pd.DataFrame(FLOW).to_csv(os.path.join(OUT, "rev05_cohort_flow.tsv"), sep="\t", index=False)
    lib.log.info("cohort flow:\n" + pd.DataFrame(FLOW).to_string(index=False))
    lib.log.info("DONE rev05")


if __name__ == "__main__":
    main()
