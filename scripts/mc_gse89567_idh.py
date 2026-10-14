#!/usr/bin/env python3
"""IDH-genotype comparison of MES-like proportion (IDH-wildtype vs IDH-mutant).

Adds the second IDH-mutant cohort to the comparison — GSE89567 (Venteicher 2017,
IDH-mutant *astrocytoma*; ``lib.GSE89567_TXT`` was defined but never used by any
released script, so the cohort the manuscript cites as "[6,10]" was only half
present in the code) — and re-derives the mutant arm at real-patient level.

Outputs
  results/multicohort/gse89567_composition.tsv      per-patient 4-state composition
  results/multicohort/idh_patients.tsv              canonical patient-level IDH table
  results/multicohort/idh_wt_vs_mut.tsv             corrected comparison (manuscript)
  results/multicohort/idh_wt_vs_mut_superseded.tsv  the superseded numbers + why

Why this script exists (the n = 259 defect)
-------------------------------------------
The mutant arm of an earlier, unpublished run pooled "patients" as
``barcode.split("_")[0]``.  GSE89567 barcodes are *dash*-separated
(``MGH64-P16-A01``), so that parser turned 247 single cells into 247 one-cell
"patients": the arm became 8 (GSE70630) + 251 (GSE89567) = **259**, giving
mean MES 0.008 and Mann-Whitney P = 1.6e-57 — all three reproduced by this
script's legacy path (``_legacy_reference()``), which is why the manuscript
carried them.  Two of the 240 surviving single cells were MES-positive, so the
P-value was inflated by n, not by biology.  ``tumour_id()`` below fixes the
parsing; the mean is essentially unchanged (0.008).  See PITFALLS in the
repository README.

Convention notes
  * GSE70630 and GSE89567 mix ``MGH<nn>_P#`` and bare ``<nn>_P#`` barcodes for
    the same tumour (GSE70630: 97/93 == MGH97/MGH93; GSE89567: 57 == MGH57,
    mgh103 == MGH103, MGH107neg/pos = two sorted fractions of MGH107), and both
    carry sorted-fraction / multiregion suffixes.  ``tumour_id()`` collapses all
    of these to the patient.
  * Cohort labels follow the *accessions* (GEO): GSE70630 = Tirosh 2016
    oligodendroglioma (PMID 27806376) and GSE89567 = Venteicher 2017 IDH-mutant
    astrocytoma (PMID 28360267) — the reverse of the "Venteicher_GSE70630" label
    used by ``mc_cohort_state.py``.
  * The IDH-wildtype arm is read as shipped from ``composition_all_cohorts.tsv``
    (n = 69: GSE103224 + GSE182109 + Neftel) so this script changes one term of
    the comparison only.  NOTE: applying ``tumour_id()`` to the Neftel cohort
    collapses 45 sample tokens to 31 tumours; that is a known open item and is
    deliberately NOT applied here.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402
from genesets import load_neftel_modules, four_state_tables, STATE_NAMES  # noqa: E402
import mc_cohort_state as mcs  # noqa: E402

OUT = os.path.join(lib.PROJ, "results", "multicohort")
log = lib.log                                          # same logger as lib.py
WILDTYPE_COHORTS = ["GSE103224", "GSE182109", "Neftel_GSE131928"]
MUTANT_COHORTS = ["GSE70630", "GSE89567"]
STATES = [s for s in STATE_NAMES]                       # AC, NPC, OPC, MES

# sorted-fraction / multiregion / plate suffixes seen in the GBM scRNA barcodes
_FRACTION_PLATE = re.compile(r"(?:CD45|CD3)?(?:NEG|POS)P\d+$", re.I)   # MGH104negP2
_CD45 = re.compile(r"CD45(?:NEG|POS)$", re.I)                          # MGH143CD45neg
_FRACTION = re.compile(r"(?:NEG|POS)$", re.I)                          # MGH107neg
_MULTIREGION = re.compile(r"(?<=[0-9])[A-D]$")                         # MGH105A -> MGH105


def tumour_id(barcode):
    """Patient/tumour id from a cell barcode, robust to the naming variants.

    ``MGH64-P16-A01`` -> MGH64 · ``MGH53_P5_C12`` -> MGH53 · ``97_P5_H07`` -> MGH97
    ``MGH104negP2`` -> MGH104 · ``MGH143CD45neg`` -> MGH143 · ``MGH105A`` -> MGH105
    ``mgh103`` -> MGH103
    """
    tok = re.split(r"[-_]", str(barcode).strip())[0]
    tok = _FRACTION_PLATE.sub("", tok)
    tok = _CD45.sub("", tok)
    tok = _FRACTION.sub("", tok)
    tok = _MULTIREGION.sub("", tok)
    tok = tok.upper()
    if re.fullmatch(r"\d+", tok):
        tok = "MGH" + tok
    return tok


def composition_from_scores(path, label):
    """Per-patient 4-state composition from a per-cell ``*_state_scores.tsv``."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df["patient"] = [tumour_id(b) for b in df["patient"]]
    counts = df.groupby(["patient", "state"]).size().unstack(fill_value=0)
    for s in STATES:
        if s not in counts:
            counts[s] = 0
    counts = counts[STATES]
    comp = counts.div(counts.sum(axis=1), axis=0)
    comp.columns = [f"{s}_prop" for s in STATES]
    comp["n_cells"] = counts.sum(axis=1)
    comp.index.name = "patient"
    comp["cohort"] = label
    return comp


def score_gse89567():
    """Score GSE89567 (Venteicher IDH-mutant astrocytoma) with the shipped pipeline."""
    log.info(f"scoring {os.path.basename(lib.GSE89567_TXT)}")
    ad = mcs.load_log_matrix(lib.GSE89567_TXT, "GSE89567", lambda c: c, "processed")
    log.info(f"  loaded {ad.n_obs} cells x {ad.n_vars} genes")
    ad = mcs.exclude_immune(ad)
    ad = mcs.score_state(ad, mcs.state_sub)
    mal = ad[ad.obs["nonmalignant"] == 0].copy()
    log.info(f"  malignant cells kept: {mal.n_obs} of {ad.n_obs}")
    return composition_from_scores(_materialise_state_scores(mal), "GSE89567")


def _materialise_state_scores(mal):
    """Write the per-cell (patient, state, score_MES) table and return its path."""
    path = os.path.join(OUT, "gse89567_state_scores.tsv")
    mal.obs[["patient", "state", "score_MES"]].to_csv(path, sep="\t")
    return path


def compare(wt, mut, label):
    u, p = mannwhitneyu(wt, mut, alternative="two-sided")
    row = {
        "comparison": label,
        "n_IDHwt": len(wt), "n_IDHmut": len(mut),
        "MES_IDHwt_mean": round(float(wt.mean()), 3),
        "MES_IDHmut_mean": round(float(mut.mean()), 3),
        "MES_IDHwt_median": round(float(wt.median()), 3),
        "MES_IDHmut_median": round(float(mut.median()), 3),
        "mannwhitney_U": float(u), "p": float(p),
    }
    log.info(f"{label}: IDH-wt n={len(wt)} mean={wt.mean():.3f} | "
             f"IDH-mut n={len(mut)} mean={mut.mean():.3f} | MWU P={p:.3g}")
    return row


def _legacy_reference(gse89567_cells):
    """Reproduce the superseded n=259 artifact (``split('_')[0]`` parsing).

    Kept as an executable regression record of *why* the manuscript's mutant arm
    read n = 259 / 0.008 / P = 1.6e-57.
    """
    tok = gse89567_cells["patient"].astype(str).str.split("_").str[0]
    df = gse89567_cells.assign(patient=tok)
    counts = df.groupby(["patient", "state"]).size().unstack(fill_value=0)
    for s in STATES:
        if s not in counts:
            counts[s] = 0
    prop = counts[STATES].div(counts[STATES].sum(axis=1), axis=0)["MES"]
    return prop


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "mc_gse89567_idh.log"))
    os.makedirs(OUT, exist_ok=True)
    mcs.state_sub = four_state_tables(load_neftel_modules())[1]

    # ---------------- IDH-wildtype arm: as shipped, untouched ----------------
    allc = pd.read_csv(os.path.join(OUT, "composition_all_cohorts.tsv"), sep="\t")
    wt_tbl = allc[allc["cohort"].isin(WILDTYPE_COHORTS)].copy()
    log.info(f"IDH-wildtype arm (as shipped): {len(wt_tbl)} patients from "
             f"{sorted(wt_tbl['cohort'].unique())}")

    # ---------------- IDH-mutant arm: GSE70630 + GSE89567, real patients ----------------
    g70630 = composition_from_scores(
        os.path.join(OUT, "Venteicher_GSE70630_state_scores.tsv"), "GSE70630")
    log.info(f"GSE70630: {g70630.shape[0]} tumours after barcode normalisation "
             f"(8 raw sample tokens)")
    g89567 = score_gse89567()
    log.info(f"GSE89567: {g89567.shape[0]} patients after barcode normalisation")
    g89567.to_csv(os.path.join(OUT, "gse89567_composition.tsv"), sep="\t")

    mut_tbl = pd.concat([g70630, g89567]).reset_index()
    mut_tbl["genotype"] = "IDH-mut"

    # canonical patient-level table (both arms, useful for figures)
    patients = pd.concat([
        wt_tbl[["patient", "cohort"]].assign(
            genotype="IDH-wt", MES_prop=wt_tbl["MES_prop"].values),
        mut_tbl[["patient", "cohort", "genotype", "MES_prop"]],
    ], ignore_index=True)
    patients.to_csv(os.path.join(OUT, "idh_patients.tsv"), sep="\t", index=False)
    log.info(f"wrote idh_patients.tsv: {len(patients)} patients "
             f"({(patients['genotype'] == 'IDH-wt').sum()} wt / "
             f"{(patients['genotype'] == 'IDH-mut').sum()} mut)")

    # ---------------- corrected comparison ----------------
    wt = wt_tbl["MES_prop"].reset_index(drop=True)
    mut_corrected = pd.concat([g70630["MES_prop"], g89567["MES_prop"]],
                              ignore_index=True)
    corrected = compare(wt, mut_corrected, "IDH_wt_vs_IDH_mut")
    pd.DataFrame([corrected]).to_csv(
        os.path.join(OUT, "idh_wt_vs_mut.tsv"), sep="\t", index=False)

    # ---------------- superseded numbers + provenance ----------------
    cells = pd.read_csv(os.path.join(OUT, "gse89567_state_scores.tsv"), sep="\t",
                        index_col=0)
    legacy_89567 = _legacy_reference(cells)
    vent_legacy = pd.read_csv(os.path.join(OUT, "gse70630_composition.tsv"),
                              sep="\t", index_col=0)["MES_prop"]
    legacy = compare(wt, pd.concat([vent_legacy.reset_index(drop=True),
                                    legacy_89567.reset_index(drop=True)],
                                   ignore_index=True),
                     "IDH_wt_vs_IDH_mut_superseded")
    legacy["note"] = ("superseded: GSE89567 pooled with barcode.split('_')[0], "
                      "which splits MGH64-P16-A01 into one 'patient' per cell; "
                      "reproduced by the same script via _legacy_reference()")
    pd.DataFrame([legacy]).to_csv(
        os.path.join(OUT, "idh_wt_vs_mut_superseded.tsv"), sep="\t", index=False)

    log.info("done: mutant arm = GSE70630 (%d tumours) + GSE89567 (%d patients) = %d"
             % (g70630.shape[0], g89567.shape[0], len(mut_corrected)))


if __name__ == "__main__":
    main()
