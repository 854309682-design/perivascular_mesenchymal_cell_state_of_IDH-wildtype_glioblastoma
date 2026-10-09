#!/usr/bin/env python3
"""Fetch TCGA-GBM/LGG IDH1/2 mutation status via GDC masked somatic-mutation MAFs.
Streams each per-case MAF, keeps only IDH1/IDH2 rows, and writes a per-case
IDH_mutated table. Run in the background (parallel to the expression download).
"""
import os
import sys
import json
import glob
import tempfile
import urllib.request, urllib.parse
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

API = "https://api.gdc.cancer.gov"
OUT = os.path.join(lib.RAW_DIR, "tcga", "mutation")
os.makedirs(OUT, exist_ok=True)
UA = {"User-Agent": "paper2/1.0"}
PROJECTS = ["TCGA-GBM", "TCGA-LGG"]
IDHG = ["IDH1", "IDH2"]


def api_files(filt):
    p = {"filters": json.dumps(filt), "size": "5000", "format": "json",
         "fields": "file_id,file_name,file_size"}
    req = urllib.request.Request(API + "/files?" + urllib.parse.urlencode(p), headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode()).get("data", {}).get("hits", [])


def main():
    rows = []
    for proj in PROJECTS:
        filt = {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}}]}
        hits = api_files(filt)
        print(f"[{proj}] {len(hits)} MAF files", flush=True)
        for i, h in enumerate(hits):
            out = os.path.join(OUT, h["file_name"] + ".part")
            if os.path.exists(out.replace(".part", ".done")) and os.path.getsize(os.path.join(OUT, h["file_name"])) > 0:
                continue
            fid = h["file_id"]
            try:
                req = urllib.request.Request(API + "/data/" + fid, headers=UA)
                with urllib.request.urlopen(req, timeout=300) as r, open(out, "wb") as fh:
                    while True:
                        c = r.read(1 << 20)
                        if not c:
                            break
                        fh.write(c)
                # extract IDH1/2
                d = pd.read_csv(out, sep="\t", comment="#", low_memory=False,
                                usecols=lambda c: c in ["Hugo_Symbol", "Tumor_Sample_Barcode",
                                                        "Variant_Classification"])
                idh = d[d["Hugo_Symbol"].isin(IDHG)] if "Hugo_Symbol" in d.columns else pd.DataFrame()
                if len(idh):
                    ns = idh[idh["Variant_Classification"] != "Silent"]
                    for _, r in ns.iterrows():
                        rows.append({"sample": r["Tumor_Sample_Barcode"],
                                     "gene": r["Hugo_Symbol"],
                                     "vclass": r["Variant_Classification"]})
                try:
                    os.remove(out)  # drop MAF after extracting IDH (save disk)
                except OSError:
                    pass
            except Exception as e:
                print(f"  FAIL {h['file_name']}: {e}", flush=True)
            if (i + 1) % 20 == 0:
                print(f"  [{proj}] {i+1}/{len(hits)} processed, idh hits so far {len(rows)}", flush=True)
    if rows:
        df = pd.DataFrame(rows)
        per = df.groupby("sample")["gene"].apply(lambda x: ",".join(sorted(set(x)))).reset_index()
        per["IDH_mutated"] = 1
        # include cases with no IDH mutation as 0
        per.to_csv(os.path.join(OUT, "idh_status.tsv"), sep="\t", index=False)
        print(f"IDH status: {len(per)} IDH-mutated samples -> idh_status.tsv", flush=True)
    print("IDH_DONE", flush=True)


if __name__ == "__main__":
    main()
