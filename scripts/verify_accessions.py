#!/usr/bin/env python3
"""W1 — Verify GEO/EGA accessions for Paper 2 and emit results/step00_*.tsv.

Uses NCBI E-utilities (db=gds) to resolve each candidate accession to title,
organism, sample count, and PubMed link. Can only verify; it never invents data.
Writes tables/datasets.tsv (canonical list, to be reviewed by Prof. Meng).
"""
import os
import sys
import json
import time
import urllib.parse
import urllib.request
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_TSV = os.path.join(lib.PROJ, "tables", "datasets.tsv")

# (dataset_key, accession, cohort_role, note)
CANDIDATES = [
    ("GSE131928", "GSE131928", "Discovery (Neftel 2019 GBM 4-state)", "verify"),
    ("GSE103224", "GSE103224", "Discovery+internal validation (Yuan 2018)", "verify"),
    ("GSE70630", "GSE70630", "IDH-mut biological control (Venteicher 2017)", "verify"),
    ("GSE182109", "GSE182109", "candidate Abdelfattah 2022 primary vs recurrence", "TO CONFIRM"),
    ("GSE163795", "GSE163795", "candidate Ravi 2022 / GBM spatial", "TO CONFIRM"),
    ("GSE189396", "GSE189396", "backup GBM spatial (10x)", "backup"),
]


def eutils_get(endpoint, params):
    url = f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "paper2-verify/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def esearch(accession):
    try:
        d = eutils_get("esearch", {"db": "gds", "term": f"{accession}[Accession]", "retmode": "json"})
        ids = d.get("esearchresult", {}).get("idlist", [])
        return ids
    except Exception as e:
        return []


def esummary(uid):
    try:
        d = eutils_get("esummary", {"db": "gds", "id": uid, "retmode": "json"})
        res = d.get("result", {}).get(uid, {})
        return {
            "title": res.get("title", ""),
            "organism": res.get("taxon", ""),
            "n_samples": res.get("n_samples", ""),
            "pubmed": res.get("pubmedids", []),
            "gdsType": res.get("gdstype", ""),
            "summary": (res.get("summary", "") or "")[:200],
        }
    except Exception as e:
        return {"error": str(e)}


def main():
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "step00_verify.log"))
    rows = []
    for key, acc, role, note in CANDIDATES:
        ids = esearch(acc)
        info = {"dataset": key, "accession": acc, "role": role, "note": note}
        if not ids:
            info["resolved"] = "NOT_FOUND"
            info["title"] = ""
            info["organism"] = ""
            info["n_samples"] = ""
            info["pubmed"] = ""
            info["summary"] = ""
            lib.log.warning(f"{acc}: NOT FOUND on GEO")
        else:
            uid = ids[0]
            s = esummary(uid)
            info["resolved"] = "FOUND"
            info.update(s)
            info["pubmed"] = ",".join(str(p) for p in (s.get("pubmed") or []))
            info["summary"] = s.get("summary", "")
            lib.log.info(f"{acc}: FOUND -> {s.get('title')} | {s.get('organism')} "
                         f"| n={s.get('n_samples')} | PMID={info['pubmed']}")
        rows.append(info)
        time.sleep(0.5)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    lib.log.info(f"Wrote {OUT_TSV} ({len(df)} rows)")


if __name__ == "__main__":
    main()
