#!/usr/bin/env python3
"""Download TCGA-GBM + TCGA-LGG RNA-seq gene-expression (open-access) + clinical
via the GDC API into data/raw/tcga/. Logs a manifest of file_id/file_name/size.
Run in background; it fetches ~700 files (~3-4 GB) at GDC speeds.
"""
import os
import sys
import json
import time
import urllib.request, urllib.parse
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

API = "https://api.gdc.cancer.gov"
ROOT = os.path.join(lib.RAW_DIR, "tcga")
os.makedirs(ROOT, exist_ok=True)
UA = {"User-Agent": "paper2/1.0"}
FIELDS = "file_id,file_name,file_size,cases.submitter_id,cases.project.project_id"

PROJECTS = ["TCGA-GBM", "TCGA-LGG"]


def api_get(params):
    url = API + "/files?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def fetch_manifest():
    rows = []
    for proj in PROJECTS:
        filt = {
            "op": "and", "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}},
                {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
                {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
                {"op": "in", "content": {"field": "access", "value": ["open"]}},
            ],
        }
        p = {"filters": json.dumps(filt), "size": "5000", "format": "json",
             "fields": FIELDS, "sort": "file_id"}
        d = api_get(p)
        for hit in d.get("data", {}).get("hits", []):
            cs = hit.get("cases") or [{}]
            subj = cs[0].get("submitter_id", "")
            proj_id = cs[0].get("project", {}).get("project_id", proj)
            rows.append({"file_id": hit["file_id"], "file_name": hit["file_name"],
                         "file_size": int(hit["file_size"]), "project": proj_id,
                         "sample": subj})
    man = pd.DataFrame(rows).drop_duplicates("file_id")
    man.to_csv(os.path.join(ROOT, "manifest.tsv"), sep="\t", index=False)
    return man


def download(man):
    ok = fail = 0
    for i, row in man.iterrows():
        out = os.path.join(ROOT, row["file_name"])
        if os.path.exists(out) and os.path.getsize(out) == row["file_size"]:
            continue
        url = api = API + "/data/" + row["file_id"]
        req = urllib.request.Request(url, headers={"User-Agent": "paper2/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r, open(out, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            ok += 1
        except Exception as e:
            fail += 1
            print(f"FAIL {row['file_name']}: {e}", flush=True)
        if (i + 1) % 25 == 0:
            print(f"progress {i+1}/{len(man)} ok={ok} fail={fail}", flush=True)
    print(f"done: ok={ok} fail={fail}", flush=True)


if __name__ == "__main__":
    lib.setup_logging(path=os.path.join(lib.PROJ, "logs", "download_tcga.log"))
    man = fetch_manifest()
    print(f"manifest: {len(man)} files, total {man['file_size'].sum()/1e9:.2f} GB", flush=True)
    download(man)
