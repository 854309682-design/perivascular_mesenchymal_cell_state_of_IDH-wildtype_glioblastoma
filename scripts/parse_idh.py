#!/usr/bin/env python3
"""Parse the downloaded (gzipped) TCGA masked somatic-mutation MAF .part files,
extract IDH1/IDH2 non-silent mutations, and write a per-case IDH status table.
The earlier fetch failed because it read the .maf.gz as plain text; this re-parses
with gzip. Deletes the .part files after parsing (they were only needed for IDH).
"""
import os
import sys
import glob
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

MUT = os.path.join(lib.RAW_DIR, "tcga", "mutation")
OUT = os.path.join(MUT, "idh_status.tsv")
IDHG = ["IDH1", "IDH2"]

parts = sorted(glob.glob(os.path.join(MUT, "*.part")))
print("part files:", len(parts))
rows = []
fails = 0
for i, p in enumerate(parts):
    try:
        d = pd.read_csv(p, sep="\t", comment="#", low_memory=False, compression="gzip",
                        usecols=lambda c: c in ["Hugo_Symbol", "Tumor_Sample_Barcode",
                                                "Variant_Classification"])
        idh = d[d["Hugo_Symbol"].isin(IDHG)] if "Hugo_Symbol" in d.columns else pd.DataFrame()
        if len(idh):
            ns = idh[idh["Variant_Classification"] != "Silent"]
            for _, r in ns.iterrows():
                rows.append({"sample": r["Tumor_Sample_Barcode"], "gene": r["Hugo_Symbol"],
                             "vclass": r["Variant_Classification"]})
        os.remove(p)
    except Exception as e:
        fails += 1
        print(f"  fail {os.path.basename(p)}: {e}")
    if (i + 1) % 200 == 0:
        print(f"  {i+1}/{len(parts)} idh hits so far {len(rows)}")

print("total IDH-mutated calls:", len(rows), "| parse fails:", fails)
if rows:
    df = pd.DataFrame(rows).drop_duplicates()
    per = df.groupby("sample")["gene"].apply(lambda x: ",".join(sorted(set(x)))).reset_index()
    per["IDH_mutated"] = 1
    per.to_csv(OUT, sep="\t", index=False)
    print(f"IDH status table: {len(per)} samples -> {OUT}")
else:
    print("no IDH mutations found in parsed MAFs")
