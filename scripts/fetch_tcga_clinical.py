#!/usr/bin/env python3
"""Efficient TCGA-GBM + TCGA-LGG clinical fetch via the GDC case API (paginated).
Returns a case-level table: submitter_id, age_at_index, gender, vital_status,
days_to_death, days_to_last_follow_up, tumor_grade, histological_type, primary_site.
Run in the background (fast, a few HTTP calls).
"""
import os
import sys
import json
import urllib.request, urllib.parse
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import lib

API = "https://api.gdc.cancer.gov"
OUT = os.path.join(lib.RAW_DIR, "tcga", "clinical")
os.makedirs(OUT, exist_ok=True)
UA = {"User-Agent": "paper2/1.0"}
PROJECTS = ["TCGA-GBM", "TCGA-LGG"]
FIELDS = ("submitter_id,demographic.age_at_index,demographic.gender,"
          "demographic.vital_status,demographic.days_to_death,"
          "demographic.days_to_last_follow_up,diagnoses.tumor_grade,"
          "diagnoses.histological_type,diagnoses.primary_site,"
          "diagnoses.morphology")


def cases(proj):
    data = []
    for page in range(0, 100):
        filt = {"op": "in", "content": {"field": "cases.project.project_id", "value": [proj]}}
        p = {"filters": json.dumps(filt), "size": "5000", "from": str(page * 5000),
             "format": "json", "fields": FIELDS}
        url = API + "/cases?" + urllib.parse.urlencode(p)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read().decode())
        hits = d.get("data", {}).get("hits", [])
        data.extend(hits)
        if len(hits) < 5000:
            break
    return data


def flatten(hits, proj):
    rows = []
    for h in hits:
        dem = h.get("demographic") or {}
        dx = (h.get("diagnoses") or [{}])[0]
        rows.append({
            "submitter_id": h.get("submitter_id", ""),
            "project": proj,
            "age_at_index": dem.get("age_at_index", ""),
            "gender": dem.get("gender", ""),
            "vital_status": dem.get("vital_status", ""),
            "days_to_death": dem.get("days_to_death", ""),
            "days_to_last_follow_up": dem.get("days_to_last_follow_up", ""),
            "tumor_grade": dx.get("tumor_grade", ""),
            "histological_type": dx.get("histological_type", ""),
            "primary_site": dx.get("primary_site", ""),
            "morphology": dx.get("morphology", ""),
        })
    return rows


def main():
    allr = []
    for proj in PROJECTS:
        h = cases(proj)
        allr += flatten(h, proj)
        print(f"[{proj}] {len(h)} cases", flush=True)
    df = pd.DataFrame(allr).drop_duplicates("submitter_id")
    df.to_csv(os.path.join(OUT, "tcga_clinical.tsv"), sep="\t", index=False)
    print(f"clinical: {len(df)} cases -> tcga_clinical.tsv; "
          f"vital_status values: {df['vital_status'].value_counts().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
