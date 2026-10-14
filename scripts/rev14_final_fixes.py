#!/usr/bin/env python3
"""rev14 — Final targeted fixes after the duplication repair.

1. Collapse the triplicated Intro positioning sentence (C/EBPbeta-STAT3 ... [29]) in pid=6109BE8D.
2. Apply the two edits that were reported as applied but did not land, taken verbatim from
   scripts/rev11_edits.py: the Methods NicheNet gene-set name (73-gene -> 128-gene) and the
   Figure 5 legend, which must describe the rebuilt panels.

Run: /home/meng/miniconda3/bin/python3 scripts/rev14_final_fixes.py [--check]
"""
import os, re, sys, zipfile, hashlib
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_final_fixes")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv
sys.path.insert(0, os.path.join(PROJ, "scripts"))
from rev11_edits import EDITS

CARRO = ("The transcriptional basis of this state was first defined by C/EBP\u03b2-STAT3 as master "
         "regulators of mesenchymal transformation [26], extended by NF-\u03baB-driven differentiation "
         "and radioresistance [27], linked in bulk tumours to microenvironmental rather than genetic "
         "evolution [28], and shown to be inducible and reversible after irradiation in preclinical "
         "models [29]. ")


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
    log = []

    # 1. collapse the triplicated Carro sentence
    m = re.search(r'<w:p\b[^>]*w14:paraId="6109BE8D".*?</w:p>', xml, re.S)
    if m:
        p = m.group(0); c = p.count(esc(CARRO))
        p2 = p
        while esc(CARRO) + esc(CARRO) in p2:
            p2 = p2.replace(esc(CARRO) + esc(CARRO), esc(CARRO))
        if p2 != p:
            xml = xml.replace(p, p2, 1)
        log.append(f"6109BE8D Carro sentence {c}x -> {p2.count(esc(CARRO))}x")

    # 2. the two unapplied edits from the table
    for e in EDITS:
        if e[0] != "sub":
            continue
        _, pid, old, new = e
        if pid not in ("4DF903D0", "50242B74"):
            continue
        m = re.search(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % pid, xml, re.S)
        if not m:
            log.append(f"{pid}: paragraph missing"); continue
        p = m.group(0)
        if esc(old) not in p:
            log.append(f"{pid}: anchor absent (old count {p.count(esc(old))}) — nothing to do")
            continue
        p2 = p.replace(esc(old), esc(new), 1)
        xml = xml.replace(p, p2, 1)
        log.append(f"{pid}: replaced ({len(old)} -> {len(new)} chars)")

    for l in log:
        print("  " + l)
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
    print(f"  73-gene left: {txt.count('73-gene')}; 128-gene programme: {txt.count('128-gene MES-transition programme')}")
    print(f"  new legend present: {'IvyGAP anatomic regions: MES (MES1+MES2)' in txt}")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {md5p(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
