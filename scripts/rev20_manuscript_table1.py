#!/usr/bin/env python3
"""rev20 — Manuscript: insert Table 1 (cohort summary) as a real Word table.

Reviewer item *minor 26* asked for a table that summarises the cohorts. rev17 already added
the in-text reference ("cohort-level detail ... is in Table 1"), so the table itself has to
exist in the file; it did not (`<w:tbl>` count was 0).

The table is built from `results/revision/rev15_table1_cohorts.tsv` (written by rev15) and is
appended on its own page after the reference list, which is the layout most journals ask for.
Seven columns keep it readable in portrait: cohort, platform, IDH genotype, patients,
malignant cells, mean MES-like proportion, treatment context. The long "analyses each cohort
contributes to" strings become a footnote paragraph under the table instead of a wrapped
column.

Same rules as every other batch in this project: paragraph-level zip surgery, only
`word/document.xml` may differ, a dated backup first, all-or-nothing, idempotent (a sentinel
check on the whole document), then structural verification. The manuscript is never rebuilt
from markdown and never re-saved with python-docx.

Run:
  /home/meng/miniconda3/bin/python3 scripts/rev20_manuscript_table1.py [--check]
"""
import os
import re
import sys
import shutil
import zipfile
import hashlib
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
TSV = os.path.join(PROJ, "results", "revision", "rev15_table1_cohorts.tsv")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_table1")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv

# tail anchor: the last reference paragraph of the list (Minata 2019 = ref 29)
TAIL_PID = "201B6ED6"
SENTINEL = "Table 1 | Cohort summary"

CAPTION_BOLD = "Table 1 | Cohort summary."
CAPTION_REST = (" Platform, IDH genotype, patient and malignant-cell counts, mean MES-like "
                "proportion and treatment context of the five single-cell cohorts; the "
                "analyses each cohort contributes to are listed beneath the table.")

COLS = [("Cohort", 1150), ("Platform", 1600), ("IDH genotype", 1250), ("Patients", 800),
        ("Malignant cells", 1150), ("Mean MES-like", 1000), ("Treatment context", 2050)]
TOTAL = sum(w for _, w in COLS)


def esc(t):
    return (str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_USED_PIDS = set()


def pid_for(text):
    """A deterministic, unique 8-hex `w14:paraId` (every paragraph in this file has one)."""
    h = hashlib.md5(("rev20|" + str(text)).encode("utf-8")).hexdigest()[:8].upper()
    cand, n = h, 0
    while cand in _USED_PIDS:
        n += 1
        cand = (h[:7] + format((int(h[7], 16) + n) % 16, "X"))
        if n > 15:
            cand = hashlib.md5((h + str(n)).encode()).hexdigest()[:8].upper()
    _USED_PIDS.add(cand)
    return cand


def seed_used_pids(xml):
    _USED_PIDS.update(re.findall(r'w14:paraId="([0-9A-Fa-f]{8})"', xml))


def md5p(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def cell(text, w, bold=False, size=18):
    rpr = "<w:rPr>" + ("<w:b/><w:bCs/>" if bold else "") + \
          '<w:color w:val="000000" w:themeColor="text1"/>' \
          f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/></w:rPr>'
    return (f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/></w:tcPr>'
            f'<w:p w14:paraId="{pid_for(text)}"><w:pPr><w:spacing w:before="20" '
            f'w:after="20" w:line="240" w:lineRule="auto"/></w:pPr>'
            f'<w:r>{rpr}<w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p></w:tc>')


def table(rows):
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for _, w in COLS)
    b = ('<w:top w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
         '<w:left w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
         '<w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
         '<w:right w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
         '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
         '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="000000"/>')
    out = [f'<w:tbl><w:tblPr><w:tblW w:w="{TOTAL}" w:type="dxa"/><w:tblBorders>{b}'
           '</w:tblBorders><w:tblLayout w:type="fixed"/>'
           '<w:tblCellMar><w:left w:w="70" w:type="dxa"/><w:right w:w="70" w:type="dxa"/>'
           '</w:tblCellMar></w:tblPr>', f'<w:tblGrid>{grid}</w:tblGrid>']
    # header
    out.append('<w:tr><w:trPr><w:tblHeader/></w:trPr>' +
               "".join(cell(h, w, bold=True) for h, w in COLS) + '</w:tr>')
    for r in rows:
        out.append('<w:tr>' + "".join(cell(v, w) for (v, (_, w)) in zip(r, COLS)) + '</w:tr>')
    out.append('</w:tbl>')
    return "".join(out)


def para(text, style="3", bold_prefix=None, page_break=False, size=None, space_after=120):
    """A body paragraph; `bold_prefix` is emitted as a separate bold run.

    `w:pPr` children must follow the schema order (pStyle, pageBreakBefore, spacing, rPr).
    """
    ppr = "<w:pPr>"
    if style:
        ppr += f'<w:pStyle w:val="{style}"/>'
    if page_break:
        ppr += "<w:pageBreakBefore/>"
    ppr += f'<w:spacing w:after="{space_after}"/>' \
           '<w:rPr><w:color w:val="000000" w:themeColor="text1"/></w:rPr></w:pPr>'
    sz = f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>' if size else ""
    runs = ""
    if bold_prefix:
        runs += ('<w:r><w:rPr><w:b/><w:bCs/>'
                 '<w:color w:val="000000" w:themeColor="text1"/>' + sz + '</w:rPr>'
                 f'<w:t xml:space="preserve">{esc(bold_prefix)}</w:t></w:r>')
    runs += ('<w:r><w:rPr><w:color w:val="000000" w:themeColor="text1"/>' + sz + '</w:rPr>'
             f'<w:t xml:space="preserve">{esc(text)}</w:t></w:r>')
    return f'<w:p w14:paraId="{pid_for(text[:48] + str(len(runs)))}">{ppr}{runs}</w:p>'


def read_rows():
    import csv
    with open(TSV, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    def num(v, nd=0):
        try:
            return f"{float(v):,.{nd}f}"
        except (TypeError, ValueError):
            return str(v).strip() or "–"

    tbl, notes = [], []
    for r in rows:
        tbl.append([r["cohort"], r["platform"], r["genotype"], num(r["n_patients"]),
                    num(r["n_malignant_cells"]), num(r["mean_MES_like"], 3),
                    r["treatment_context"]])
        notes.append(f'{r["cohort"]}: {r["analyses"]}')
    return tbl, notes


def main():
    src = md5p(SRC)
    z = zipfile.ZipFile(SRC)
    names = z.namelist()
    ent = {n: z.read(n) for n in names}
    xml = ent[DOC].decode("utf-8")
    text_all = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', xml, re.S))
    print(f"source {SRC}\n  md5 {src}; <w:tbl> before = {xml.count('<w:tbl>')}")

    if SENTINEL in text_all:
        print("\nalready applied (sentinel present) — nothing to do.")
        return 0

    m = re.search(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % TAIL_PID, xml, re.S)
    if not m:
        print(f"\nNOT WRITING: tail anchor paragraph {TAIL_PID} not found.")
        return 1
    tail = m.group(0)
    if "Minata M. et al." not in tail:
        print(f"\nNOT WRITING: paragraph {TAIL_PID} is not the expected reference.")
        return 1

    seed_used_pids(xml)
    tbl, notes = read_rows()
    if len(tbl) != 5:
        print(f"\nNOT WRITING: expected 5 cohort rows, got {len(tbl)}.")
        return 1

    block = (para(CAPTION_REST, bold_prefix=CAPTION_BOLD, page_break=True, space_after=160) +
             table(tbl) +
             para("Analyses contributed. " + " ".join(n + "." for n in notes), space_after=0))
    new_xml = xml.replace(tail, tail + block, 1)
    if new_xml == xml:
        print("\nNOT WRITING: replacement produced no change.")
        return 1
    try:
        import xml.etree.ElementTree as ET
        ET.fromstring(new_xml)
    except Exception as e:                                   # malformed XML must never ship
        print(f"\nNOT WRITING: generated document.xml is not well-formed ({e}).")
        return 1

    print(f"\ninserting after {TAIL_PID}: caption + {len(tbl)}x{len(COLS)} table + footnote "
          f"({len(block)} XML chars)")
    if CHECK:
        print("--check: no write performed.")
        return 0

    os.makedirs(BAK, exist_ok=True)
    shutil.copy2(SRC, os.path.join(BAK, os.path.basename(SRC)))
    with open(os.path.join(BAK, "md5_before.txt"), "w") as f:
        f.write(src + "\n")

    tmp = SRC + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for n in names:
            zo.writestr(z.getinfo(n), new_xml.encode("utf-8") if n == DOC else ent[n])
    os.replace(tmp, SRC)

    z2 = zipfile.ZipFile(SRC)
    names2 = z2.namelist()
    x2 = z2.read(DOC).decode("utf-8")
    txt2 = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', x2, re.S))
    pids = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', x2)
    print("\nVERIFY:")
    print(f"  entries identical: {names2 == names}; differing: "
          f"{[n for n in names2 if z2.read(n) != ent[n]]}")
    print(f"  <w:tbl> {xml.count('<w:tbl>')} -> {x2.count('<w:tbl>')}; "
          f"<w:p> {len(re.findall(r'<w:p\b', xml))} -> {len(re.findall(r'<w:p\b', x2))}; "
          f"paraIds {len(re.findall(r'w14:paraId=', xml))} -> {len(pids)}; "
          f"paraIds unique {len(pids) == len(set(pids))}; "
          f"paragraphs without paraId {len(re.findall(r'<w:p\b', x2)) - len(pids)}")
    print(f"  bookmarks {xml.count('<w:bookmarkStart')} -> {x2.count('<w:bookmarkStart')}")
    print(f"  sentinel present: {SENTINEL in txt2}; footnote present: "
          f"{'Analyses contributed.' in txt2}")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {md5p(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
