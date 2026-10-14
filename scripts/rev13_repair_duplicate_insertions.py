#!/usr/bin/env python3
"""rev13 — Repair: collapse insertions that were applied more than once.

Cause: the idempotency probe of `rev11_manuscript_revision.py` compared the *whole* new string
against the paragraph text, normalised for typography, and was run again after the typographic
pass. When an edit's replacement RETAINS its own anchor (`new` contains `old`), a probe
mismatch makes the edit look unapplied, and the anchor still being present lets the append
happen again — so the appended tail was inserted three times in two paragraphs:

    pid=050C4C6F  the IDH-confound sentences of Results
    pid=07505CF4  the "Further limitations specific to the revised analyses" block

This script removes the extra copies (keeping the first), verifying paragraph count, paraId
uniqueness, bookmark count, entry integrity and that each tail now occurs exactly once.

Run: /home/meng/miniconda3/bin/python3 scripts/rev13_repair_duplicate_insertions.py [--check]
"""
import os
import re
import sys
import zipfile
import hashlib
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_dup_repair")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv

# (paragraph paraId, the duplicated tail)
TARGETS = [
    ("050C4C6F", "A cohort-stratified van Elteren test and the per-cohort effects are reported in "
                 "Supplementary Fig. S2f; the direction is unchanged in every cohort. "),
    ("07505CF4", "Further limitations specific to the revised analyses: the de-novo arm was re-run "
                 "with a smaller state-stratified subsample (2,000 cells per cohort) than our first "
                 "analysis because the original code was not recoverable, so its ranking supersedes "
                 "rather than reproduces the submitted one; the IvyGAP deconvolution is a "
                 "signature-matrix non-negative least-squares fit, not a full probabilistic "
                 "deconvolution; the bulk MES score correlates with WHO grade and with stromal "
                 "content, and its prognostic value in TCGA survives neither adjustment while in "
                 "CGGA it survives grade adjustment with a smaller effect size; the TCGA arm includes "
                 "low-grade glioma and is therefore not a GBM-only analysis; the submitted IDH "
                 "comparison and the submitted survival numbers were built from clinical fields "
                 "whose event/time coding disagrees with the curated tables for about a fifth of "
                 "cases, which is why every number in this version is recomputed and stated here; "
                 "and per-cell statistics remain descriptive throughout, with patient-level "
                 "estimates reported alongside them. "),
]


def md5p(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fragments():
    """(pid, duplicated fragment) pairs derived from the edit table itself."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rev11_edits import EDITS
    out = []
    for e in EDITS:
        if e[0] != "sub":
            continue
        _, pid, old, new = e
        if old in new:
            if new.startswith(old):                 # appended tail
                out.append((pid, new[len(old):]))
            elif new.endswith(old):                 # prepended head
                out.append((pid, new[: len(new) - len(old)]))
    return out


def main():
    src_md5 = md5p(SRC)
    z = zipfile.ZipFile(SRC)
    names = z.namelist()
    ent = {n: z.read(n) for n in names}
    xml = ent[DOC].decode("utf-8")
    print(f"source md5 {src_md5}")
    total_fixed = 0
    TARGETS = fragments()
    for pid, tail in TARGETS:
        m = re.search(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % pid, xml, re.S)
        if not m:
            print(f"  {pid}: paragraph not found"); continue
        p = m.group(0)
        esc = tail.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        c = p.count(esc)
        print(f"  {pid}: tail occurs {c}x")
        if c <= 1:
            continue
        p2 = p.replace(esc + esc, esc)
        while esc + esc in p2:
            p2 = p2.replace(esc + esc, esc)
        xml = xml.replace(p, p2, 1)
        total_fixed += 1
    print(f"paragraphs repaired: {total_fixed}")
    if CHECK:
        print("--check: no write"); return 0
    if total_fixed == 0:
        return 0
    os.makedirs(BAK, exist_ok=True)
    with open(os.path.join(BAK, "md5_before.txt"), "w") as f:
        f.write(src_md5 + "\n")
    tmp = SRC + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for n in names:
            zo.writestr(z.getinfo(n), xml.encode("utf-8") if n == DOC else ent[n])
    os.replace(tmp, SRC)
    z2 = zipfile.ZipFile(SRC)
    x2 = z2.read(DOC).decode("utf-8")
    pids = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', x2)
    print("\nVERIFY:")
    print(f"  entries identical: {z2.namelist() == names}; differing: "
          f"{[n for n in z2.namelist() if z2.read(n) != ent[n]]}")
    print(f"  paragraphs {len(pids)}; paraIds unique: {len(pids) == len(set(pids))}; "
          f"bookmarks {x2.count('<w:bookmarkStart')}")
    TARGETS = fragments()
    for pid, tail in TARGETS:
        esc = tail.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        m = re.search(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % pid, x2, re.S)
        print(f"  {pid}: tail now occurs {m.group(0).count(esc)}x")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {md5p(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
