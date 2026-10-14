#!/usr/bin/env python3
"""rev17 — Manuscript text: honest IDH statement (replaces an impossible test) + Table 1.

1. The revised Results had promised "a cohort-stratified van Elteren test ... Supplementary
   Fig. S2f". rev15 showed that test cannot exist: every cohort is single-genotype, so all
   strata are one-armed. The sentence is replaced by what was actually computed (per-cohort
   means, a cohort-mean-centred test, P = 0.095) so the manuscript never promises an
   unidentifiable analysis.
2. Adds the cohort summary (reviewer minor 26) as a Table 1 reference in the cohort paragraph.

Anchors verified against the current file; all-or-nothing write.

Run: /home/meng/miniconda3/bin/python3 scripts/rev17_manuscript_idh_table1.py [--check]
"""
import os, re, sys, zipfile, hashlib
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_idh_table1")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv

EDITS = [
    # 1. IDH: replace the promise of a van Elteren test with the analyses that exist
    ("050C4C6F",
     "A cohort-stratified van Elteren test and the per-cohort effects are reported in "
     "Supplementary Fig. S2f; the direction is unchanged in every cohort.",
     "Because each series is single-genotype there is no within-cohort (van Elteren) test to "
     "run — every stratum is one-armed — so we report per-cohort means instead, which separate "
     "in the same direction in every cohort (IDH-wildtype 0.185–0.496 versus IDH-mutant "
     "0.002–0.011; Supplementary Fig. S2e), together with a cohort-mean-centred test, which is "
     "not significant (P = 0.095) and therefore shows that part of the contrast is a "
     "between-cohort difference rather than a within-cohort genotype effect."),
    # 2. Table 1 reference in the cohort paragraph
    ("050C4C6F",
     "We next extended state scoring to four additional single-cell cohorts (GSE131928 [1], "
     "GSE182109 [5], GSE70630 [10], GSE89567 [6]) to ask what determines the size of the MES "
     "compartment (Fig. 2).",
     "We next extended state scoring to four additional single-cell cohorts (GSE131928 [1], "
     "GSE182109 [5], GSE70630 [10], GSE89567 [6]) to ask what determines the size of the MES "
     "compartment (Fig. 2); cohort-level detail — platform, patient and malignant-cell counts, "
     "IDH composition, treatment context and the analyses each cohort contributes to — is in "
     "Table 1."),
]


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def md5p(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    z = zipfile.ZipFile(SRC); names = z.namelist(); ent = {n: z.read(n) for n in names}
    xml = ent[DOC].decode("utf-8"); src = md5p(SRC)
    print(f"source md5 {src}")
    fails = []
    for pid, old, new in EDITS:
        m = re.search(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % pid, xml, re.S)
        if not m:
            fails.append((pid, "paragraph missing")); continue
        p = m.group(0)
        if esc(new) in p:
            print(f"  {pid}: already applied"); continue
        if esc(old) not in p:
            fails.append((pid, "anchor absent")); continue
        xml = xml.replace(p, p.replace(esc(old), esc(new), 1), 1)
        print(f"  {pid}: replaced ({len(old)} -> {len(new)} chars)")
    if fails:
        print("NOT WRITING:", fails); return 1
    if CHECK:
        print("--check: no write"); return 0
    os.makedirs(BAK, exist_ok=True)
    with open(os.path.join(BAK, "md5_before.txt"), "w") as f:
        f.write(src + "\n")
    tmp = SRC + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for n in names:
            zo.writestr(z.getinfo(n), xml.encode("utf-8") if n == DOC else ent[n])
    os.replace(tmp, SRC)
    z2 = zipfile.ZipFile(SRC); x2 = z2.read(DOC).decode("utf-8")
    txt = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', x2, re.S))
    pids = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', x2)
    print("\nVERIFY:")
    print(f"  entries identical: {z2.namelist() == names}; differing: {[n for n in z2.namelist() if z2.read(n) != ent[n]]}")
    print(f"  paragraphs {len(pids)}; paraIds unique {len(pids) == len(set(pids))}; bookmarks {x2.count('<w:bookmarkStart')}")
    print(f"  van Elteren promise left: {txt.count('van Elteren test and the per-cohort')}; "
          f"Table 1 referenced: {'Table 1.' in txt}")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {md5p(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
