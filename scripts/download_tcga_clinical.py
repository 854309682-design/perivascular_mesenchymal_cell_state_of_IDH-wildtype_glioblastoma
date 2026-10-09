#!/usr/bin/env python3
"""Fetch TCGA-GBM + TCGA-LGG clinical (survival/age/grade) and IDH1/2 mutation
status via the GDC API into data/raw/tcga/clinical/ and data/raw/tcga/mutation/.
Parses a case-level clinical table and extracts per-case IDH1/2 mutation status.
Run in the background (parallel to the expression download).
"""
import os
import sys
import json
import glob
import urllib.request, urllib.parse
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

API = "https://api.gdc.cancer.gov"
CLIN = os.path.join(lib.RAW_DIR, "tcga", "clinical")
MUT = os.path.join(lib.RAW_DIR, "tcga", "mutation")
os.makedirs(CLIN, exist_ok=True); os.makedirs(MUT, exist_ok=True)
UA = {"User-Agent": "paper2/1.0"}
PROJECTS = ["TCGA-GBM", "TCGA-LGG"]
IDHG = ["IDH1", "IDH2"]


def api_files(filt, fields="file_id,file_name,file_size,data_type,data_format"):
    p = {"filters": json.dumps(filt), "size": "5000", "format": "json", "fields": fields}
    url = API + "/files?" + urllib.parse.urlencode(p)
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode()).get("data", {}).get("hits", [])


def download_file(fid, out):
    req = urllib.request.Request(API + "/data/" + fid, headers=UA)
    with urllib.request.urlopen(req, timeout=600) as r, open(out, "wb") as fh:
        while True:
            c = r.read(1 << 20)
            if not c:
                break
            fh.write(c)
    return out


def fetch(filt, outdir, tag):
    hits = api_files(filt)
    print(f"[{tag}] {len(hits)} files", flush=True)
    for h in hits:
        out = os.path.join(outdir, h["file_name"])
        if os.path.exists(out) and os.path.getsize(out) == h["file_size"]:
            continue
        try:
            download_file(h["file_id"], out)
        except Exception as e:
            print(f"  FAIL {h['file_name']}: {e}", flush=True)
    print(f"[{tag}] done", flush=True)


def main():
    import time
    for proj in PROJECTS:
        # clinical supplement
        fetch({"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Clinical"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}}]},
            CLIN, f"{proj}_clinical")
        # masked somatic mutation (for IDH status)
        fetch({"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}}]},
            MUT, f"{proj}_mutation")

    # ---- parse clinical supplement(s) into a case table ----
    clin_files = glob.glob(os.path.join(CLIN, "*.tsv")) + glob.glob(os.path.join(CLIN, "*.tsv.gz"))
    rows = []
    print("[parse] clinical files:", len(clin_files), flush=True)
    # clinical supplement format: first columns are metadata header lines; find the data header
    for cf in clin_files:
        try:
            df = pd.read_csv(cf, sep="\t", comment="#", low_memory=False)
        except Exception as e:
            print(f"  parse skip {os.path.basename(cf)}: {e}", flush=True)
            continue
        if "days_to_death" in df.columns or "vital_status" in df.columns:
            keep = [c for c in ["submitter_id", "case_submitter_id", "age_at_index", "vital_status",
                               "days_to_death", "days_to_last_follow_up", "morphology", "histological_type",
                               "tumor_grade", "primary_site", "gender"] if c in df.columns]
            sub = df[keep] if keep else df
            if "case_submitter_id" not in sub.columns:
                sub["case_submitter_id"] = sub.get("submitter_id", "")
            sub["project"] = os.path.basename(cf)
            rows.append(sub)
    if rows:
        clin = pd.concat(rows, ignore_index=True)
        clin.to_csv(os.path.join(CLIN, "tcga_clinical.tsv"), sep="\t", index=False)
        print(f"[parse] clinical cases: {len(clin)} -> tcga_clinical.tsv", flush=True)

    # ---- extract IDH1/2 mutation status from MAFs ----
    maf = glob.glob(os.path.join(MUT, "*"))
    idh_rows = []
    for mf in maf:
        try:
            d = pd.read_csv(mf, sep="\t", comment="#", low_memory=False, usecols=lambda c: c in
                            ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
        except Exception:
            continue
        if "Hugo_Symbol" not in d.columns:
            continue
        d = d[d["Hugo_Symbol"].isin(IDHG)]
        if len(d):
            idh_rows.append(d)
    if idh_rows:
        idh = pd.concat(idh_rows, ignore_index=True)
        # per case: marker if any IDH1/2 non-silent mutation
        ns = idh[idh["Variant_Classification"] != "Silent"]
        per_case = ns.groupby("Tumor_Sample_Barcode")["Hugo_Symbol"].apply(lambda x: ",".join(sorted(set(x))))
        idh_status = pd.DataFrame({"IDH_mutated_genes": per_case}).reset_index()
        idh_status["IDH_mutated"] = idh_status["IDH_mutated_genes"].notna().astype(int)
        idh_status.to_csv(os.path.join(MUT, "idh_status.tsv"), sep="\t", index=False)
        print(f"[parse] IDH status: {len(idh_status)} samples, "
              f"{int(idh_status.IDH_mutated.sum())} IDH-mutated", flush=True)
    print("CLINICAL_DONE", flush=True)


if __name__ == "__main__":
    main()
