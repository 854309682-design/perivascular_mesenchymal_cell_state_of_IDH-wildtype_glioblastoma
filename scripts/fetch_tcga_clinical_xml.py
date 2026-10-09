#!/usr/bin/env python3
"""Fetch per-case GDC BCR clinical XML for TCGA-GBM/LGG and extract follow-up time
(days_to_last_followup / days_to_death / vital_status) to enable a censored OS.
Parses XML namespaces generically (local-name match). Run in the background.
"""
import os
import sys
import json
import re
import glob
import urllib.request, urllib.parse
import xml.etree.ElementTree as ET
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

API = "https://api.gdc.cancer.gov"
CLIN = os.path.join(lib.RAW_DIR, "tcga", "clinical_xml")
os.makedirs(CLIN, exist_ok=True)
UA = {"User-Agent": "paper2/1.0"}
PROJECTS = ["TCGA-GBM", "TCGA-LGG"]


def api_files(filt):
    p = {"filters": json.dumps(filt), "size": "5000", "format": "json", "fields": "file_id,file_name,file_size"}
    req = urllib.request.Request(API + "/files?" + urllib.parse.urlencode(p), headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode()).get("data", {}).get("hits", [])


def parse_xml(path):
    """Return dict of case-level survival fields (max followup, vital, days_to_death)."""
    tree = ET.parse(path)
    root = tree.getroot()
    d = {"days_to_last_followup": [], "days_to_death": [], "vital_status": None,
         "submitter_id": None}
    for el in root.iter():
        tag = el.tag.rsplit("}", 1)[-1].split(":")[-1]  # local name
        txt = (el.text or "").strip()
        if not txt:
            continue
        if tag in ("days_to_last_followup", "days_to_last_follow_up"):
            try:
                d["days_to_last_followup"].append(float(txt))
            except ValueError:
                pass
        elif tag == "days_to_death":
            try:
                d["days_to_death"].append(float(txt))
            except ValueError:
                pass
        elif tag == "vital_status" and d["vital_status"] is None:
            d["vital_status"] = txt
        elif tag == "bcr_patient_barcode" and d["submitter_id"] is None:
            d["submitter_id"] = txt
    return {"submitter_id": d["submitter_id"],
            "vital_status": d["vital_status"],
            "days_to_death": max(d["days_to_death"]) if d["days_to_death"] else None,
            "days_to_last_followup": max(d["days_to_last_followup"]) if d["days_to_last_followup"] else None}


def main():
    rows = []
    for proj in PROJECTS:
        filt = {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Clinical"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}}]}
        hits = api_files(filt)
        print(f"[{proj}] {len(hits)} clinical XML files", flush=True)
        for i, h in enumerate(hits):
            out = os.path.join(CLIN, h["file_name"])
            if os.path.exists(out) and os.path.getsize(out) == h["file_size"]:
                try:
                    rows.append(parse_xml(out)); continue
                except Exception:
                    pass
            try:
                req = urllib.request.Request(API + "/data/" + h["file_id"], headers=UA)
                with urllib.request.urlopen(req, timeout=300) as r, open(out, "wb") as fh:
                    while True:
                        c = r.read(1 << 20)
                        if not c:
                            break
                        fh.write(c)
                rows.append(parse_xml(out))
            except Exception as e:
                print(f"  FAIL {h['file_name']}: {e}", flush=True)
            if (i + 1) % 50 == 0:
                print(f"  [{proj}] {i+1}/{len(hits)} parsed rows={len(rows)}", flush=True)
    df = pd.DataFrame(rows).drop_duplicates("submitter_id")
    df = df[df["submitter_id"].notna()]
    df.to_csv(os.path.join(CLIN, "clinical_xml_survival.tsv"), sep="\t", index=False)
    print(f"clinical XML survival: {len(df)} cases | vital_status counts: "
          f"{df['vital_status'].value_counts(dropna=False).to_dict()}", flush=True)
    print("XML_DONE", flush=True)


if __name__ == "__main__":
    main()
