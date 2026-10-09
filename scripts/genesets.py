#!/usr/bin/env python3
"""Gene set / module utilities for Paper 2.

Loads the Neftel 2019 glioblastoma meta-modules (MES1, MES2, AC, OPC, NPC1,
NPC2) plus the cell-cycle modules (G1/S, G2/M) from tables/neftel_meta_modules.tsv.
The TSV has module names as the header row; each column is a gene list (empty /
'NA' entries are trailing padding).
"""
import os
import pandas as pd

# Portable: same convention as lib.py (override with GLIOMA_PROJ if the tables/
# directory is not next to this file's parent).
PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEFTEL_TSV = os.path.join(PROJ, "tables", "neftel_meta_modules.tsv")

# four-state mapping
STATE_MODULES = {
    "MES": ["MES1", "MES2"],
    "NPC": ["NPC1", "NPC2"],
    "AC": ["AC"],
    "OPC": ["OPC"],
}
STATE_NAMES = ["MES", "NPC", "AC", "OPC"]
# submodule split (for MES1 vs MES2, NPC1 vs NPC2)
SUBMODULE_MAP = {
    "MES1": ["MES1"], "MES2": ["MES2"],
    "NPC1": ["NPC1"], "NPC2": ["NPC2"],
    "AC": ["AC"], "OPC": ["OPC"],
}


def load_neftel_modules(path=NEFTEL_TSV):
    """Return dict module_name -> list of genes (dropping empty/'NA')."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    modules = {}
    for col in df.columns:
        genes = [g.strip() for g in df[col].dropna().tolist()
                 if str(g).strip() not in ("", "NA", "nan", "N/A")]
        modules[col] = genes
    return modules


def four_state_tables(modules):
    """Return (state->genes for the 4 states, submodule->genes)."""
    state_genes = {}
    sub_genes = {}
    for state, subs in STATE_MODULES.items():
        g = []
        for s in subs:
            g += modules.get(s, [])
            sub_genes[s] = modules.get(s, [])
        state_genes[state] = sorted(set(g))
    for m, genes in modules.items():
        if m in ("G1/S", "G2/M"):
            continue
        sub_genes[m] = genes
    sub_genes = {k: v for k, v in sub_genes.items() if k in SUBMODULE_MAP}
    sub_genes["MES"] = list(set(modules.get("MES1", []) + modules.get("MES2", [])))
    sub_genes["NPC"] = list(set(modules.get("NPC1", []) + modules.get("NPC2", [])))
    sub_genes["AC"] = modules.get("AC", [])
    sub_genes["OPC"] = modules.get("OPC", [])
    return state_genes, sub_genes


if __name__ == "__main__":
    m = load_neftel_modules()
    sg, sub = four_state_tables(m)
    for k in ["MES1", "MES2", "AC", "OPC", "NPC1", "NPC2"]:
        print(f"{k:5s} {len(sub.get(k, [])):3d} genes")
    print("four-state counts:")
    for k, v in sg.items():
        print(f"  {k:4s} {len(v):3d}")
