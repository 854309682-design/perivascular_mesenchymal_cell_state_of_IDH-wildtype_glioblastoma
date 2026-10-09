#!/usr/bin/env python3
"""Re-download TCGA expression files missing from disk (transient SSL failures)."""
import os
import sys
import glob
import urllib.request
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

ROOT = os.path.join(lib.RAW_DIR, "tcga")
UA = {"User-Agent": "paper2/1.0"}
man = pd.read_csv(os.path.join(ROOT, "manifest.tsv"), sep="\t")
have = set(os.path.basename(f) for f in glob.glob(os.path.join(ROOT, "*.tsv")))
miss = man[~man["file_name"].isin(have)]
print("missing expression files:", len(miss))
for _, r in miss.iterrows():
    out = os.path.join(ROOT, r["file_name"])
    try:
        req = urllib.request.Request("https://api.gdc.cancer.gov/data/" + r["file_id"], headers=UA)
        with urllib.request.urlopen(req, timeout=300) as resp, open(out, "wb") as fh:
            while True:
                c = resp.read(1 << 20)
                if not c:
                    break
                fh.write(c)
        print("downloaded", r["file_name"])
    except Exception as e:
        print("FAIL", r["file_name"], e)
print("done")
