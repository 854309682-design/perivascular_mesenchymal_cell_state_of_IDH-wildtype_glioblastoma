#!/usr/bin/env python3
"""rev12 — Typographic normalisation of the revised manuscript (own zip-surgery pass).

1. Gives the four appended reference paragraphs fresh, unique `w14:paraId` values
   (they were cloned from the last reference and would otherwise duplicate its id).
2. Normalises the notation introduced by the revision to the manuscript's existing
   conventions: proportional to `× 10⁻⁷` superscripts, `C/EBPβ`, `NF-κB`, and en dashes
   in numeric ranges (`1.02-1.30` -> `1.02–1.30`).

Verified separately: the manuscript contains **no** literal LaTeX-style double hyphens
(`--` count = 0 on the raw XML), so review item 5.1 needs no change; the review report's
`C/EBP--AP-1--TWIST1` / `Mann--Whitney` examples appear to be an artefact of the tool the
reviewer used to read the file.

Run: /home/meng/miniconda3/bin/python3 scripts/rev12_manuscript_normalise.py [--check]
"""
import os
import re
import sys
import zipfile
import hashlib
import random
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_normalise")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv
SUP = {"-3": "⁻³", "-5": "⁻⁵", "-7": "⁻⁷", "-14": "⁻¹⁴", "-30": "⁻³⁰"}
REPL = [
    ("C/EBPbeta", "C/EBPβ"),
    ("NF-kB", "NF-κB"),
    ("1.2 x 10-7", "1.2 × 10⁻⁷"),
]


def md5p(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fresh_pid(used):
    while True:
        v = "%08X" % random.randrange(0x10000000, 0xFFFFFFFF)
        if v not in used:
            used.add(v)
            return v


def main():
    src_md5 = md5p(SRC)
    z = zipfile.ZipFile(SRC)
    names = z.namelist()
    ent = {n: z.read(n) for n in names}
    xml = ent[DOC].decode("utf-8")
    print(f"source md5 {src_md5}")

    # ---- 1. unique paraIds for the appended references ----
    pids = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', xml)
    used = set(pids)
    dups = {p for p in pids if pids.count(p) > 1}
    n_fixed = 0
    if dups:
        # keep the first occurrence of each duplicated id; renumber the rest
        seen = set()
        out = []
        idx = 0
        for m in re.finditer(r'w14:paraId="([0-9A-Fa-f]+)"', xml):
            pid = m.group(1)
            if pid in dups and pid in seen:
                new = fresh_pid(used)
                out.append((m.start(1), m.end(1), new))
                n_fixed += 1
            else:
                seen.add(pid)
        for s, e, new in reversed(out):
            xml = xml[:s] + new + xml[e:]
    print(f"duplicate paraIds renumbered: {n_fixed} (was {sorted(dups)})")

    # ---- 2. notation normalisation ----
    counts = {}
    for old, new in REPL:
        c = xml.count(old)
        counts[f"{old} -> {new}"] = c
        xml = xml.replace(old, new)
    # numeric ranges -> en dash (only where both sides are numbers, never after '10-7' style exponents)
    def en_dash(m):
        return f"{m.group(1)}\u2013{m.group(2)}"
    before = len(re.findall(r"(?<!10)\b(\d+\.\d+)-(\d+\.\d+)\b", xml))
    xml = re.sub(r"(?<!10)\b(\d+\.\d+)-(\d+\.\d+)\b", en_dash, xml)
    counts["decimal range a-b -> en dash"] = before
    b2 = len(re.findall(r"\b(3\.11)-(3\.12)\b", xml))
    xml = re.sub(r"\b(3\.11)-(3\.12)\b", en_dash, xml)
    counts["3.11-3.12 -> en dash"] = b2

    for k, v in counts.items():
        print(f"  {k}: {v}")
    if CHECK:
        print("--check: no write")
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
    names2 = z2.namelist()
    x2 = z2.read(DOC).decode("utf-8")
    pids2 = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', x2)
    txt = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', x2, re.S))
    print("\nVERIFY:")
    print(f"  entries identical: {names2 == names}; differing: "
          f"{[n for n in names2 if z2.read(n) != ent[n]]}")
    print(f"  paragraph count {len(pids2)}; paraIds unique: {len(pids2) == len(set(pids2))}")
    print(f"  bookmarks {x2.count('<w:bookmarkStart')} (was 22)")
    print(f"  residual 'x 10-' patterns: {len(re.findall(r'[x×] 10-', txt))}; "
          f"'--' count: {txt.count('--')}; 'C/EBPbeta' left: {txt.count('C/EBPbeta')}")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {md5p(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
