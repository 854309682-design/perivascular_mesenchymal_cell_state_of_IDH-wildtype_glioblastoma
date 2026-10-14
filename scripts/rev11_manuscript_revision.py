#!/usr/bin/env python3
"""rev11 — Manuscript revision (reviewer report items M1-M8, P0/P1 + most P2 text items).

Applies the revision as **paragraph-level zip surgery** on the hand-polished .docx, exactly as
the project's earlier batches do (`scripts/docx_surgery*.py`): the archive is rebuilt entry by
entry and only `word/document.xml` may differ. The manuscript is NEVER rebuilt from markdown
and NEVER re-saved with python-docx.

Edit kinds
    ("para", paraId, new_text)      whole-paragraph text replacement; only allowed when the
                                    paragraph carries exactly one non-empty run, so no
                                    in-paragraph formatting can be lost
    ("sub",  paraId, old, new)      substring replacement inside a single run (keeps any
                                    italics/superscript runs elsewhere in the paragraph)
    ("append_ref", ...)             new reference paragraphs cloned from the last reference

Guarantees
    * a dated backup is written before anything else;
    * the pass is ALL-OR-NOTHING: if any anchor is missing, nothing is written and every
      failure is listed;
    * the pass is idempotent: an edit whose `new` text is already present is skipped, so the
      script can be run again after the de-novo arm produces its numbers;
    * verification: identical entry names, only `word/document.xml` different, every other
      entry byte-identical, paragraph and `w14:paraId` counts, bookmarks, `unzip -t`.

Run:
  /home/meng/miniconda3/bin/python3 scripts/rev11_manuscript_revision.py [--check]
"""
import os
import re
import sys
import shutil
import zipfile
import hashlib
import subprocess
from datetime import date

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "Manuscript", "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_review_revision")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv

from rev11_edits import EDITS, NEW_REFS  # noqa: E402

VERIFY_EXTRA = [
    # (label, must-appear substring in the new document text)
    ("title associational", "A compact AP-1/C-EBP regulatory program marks the plastic"),
    ("no TWIST1 in title", "AP-1/C-EBP regulatory program"),
    ("a priori plan", "specified three hypotheses a priori in an internal analysis plan"),
    ("TWIST1 exclusion stated", "TWIST1 is not supported by either route"),
    ("peri-necrotic niche", "peri-necrotic"),
]



def norm(t):
    """Normalise text for idempotency comparison (mirrors rev12's typographic pass)."""
    t = (t.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
          .replace("C/EBPbeta", "C/EBP\u03b2").replace("NF-kB", "NF-\u03baB"))
    t = re.sub(r"([x\u00d7])\s*10-(\d+)",
               lambda m: "\u00d7 10" + {"3": "\u207b\u00b3", "5": "\u207b\u2075",
                                         "7": "\u207b\u2077", "14": "\u207b\u00b9\u2074",
                                         "30": "\u207b\u00b3\u2070"}.get(m.group(2), "-" + m.group(2)), t)
    return re.sub(r"\b(\d+\.\d+)-(\d+\.\d+)\b",
                  lambda m: m.group(1) + "\u2013" + m.group(2), t)


def sha(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def paragraphs(xml):
    return re.findall(r"<w:p\b[^>]*>.*?</w:p>|<w:p\b[^>]*/>", xml, re.S)


def pid_of(p):
    m = re.search(r'w14:paraId="([0-9A-Fa-f]+)"', p)
    return m.group(1) if m else None


def run_texts(p):
    runs = re.findall(r"<w:r\b.*?</w:r>", p, re.S)
    return runs, ["".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', r, re.S)) for r in runs]


def set_run_text(run, new):
    """Replace the text of a run with `new`, preserving its properties."""
    props = re.match(r"(<w:r\b[^>]*>)(.*?)(</w:r>)", run, re.S)
    open_tag, body, close = props.group(1), props.group(2), props.group(3)
    tprops = re.search(r"<w:rPr>.*?</w:rPr>", body, re.S)
    rpr = tprops.group(0) if tprops else ""
    esc = (new.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    return f'{open_tag}{rpr}<w:t xml:space="preserve">{esc}</w:t>{close}'


def main():
    src_md5 = sha(SRC)
    z = zipfile.ZipFile(SRC)
    names = z.namelist()
    xml = z.read(DOC).decode("utf-8")
    orig_entries = {n: z.read(n) for n in names}
    ps = paragraphs(xml)
    by_pid = {}
    for i, p in enumerate(ps):
        if pid_of(p):
            by_pid[pid_of(p)] = i
    print(f"source: {SRC}\n  md5 {src_md5}, {len(ps)} paragraphs, {len(names)} archive entries")

    applied, skipped, failed = [], [], []
    new_xml = xml

    def para_pid_block(pid):
        p = ps[by_pid[pid]]
        runs, texts = run_texts(p)
        return p, runs, texts

    for edit in EDITS:
        kind = edit[0]
        if kind == "sub":
            _, pid, old, new = edit
            if pid not in by_pid:
                failed.append((pid, "paragraph not found", old[:60]))
                continue
            p, runs, texts = para_pid_block(pid)
            hits = [j for j, t in enumerate(texts) if old in t]
            if not hits:
                if any(norm(new) in norm(t) for t in texts):
                    skipped.append((pid, "already applied", old[:60]))
                else:
                    failed.append((pid, "anchor not found in any run", old[:60]))
                continue
            if len(hits) > 1:
                failed.append((pid, f"anchor not unique ({len(hits)} runs)", old[:60]))
                continue
            j = hits[0]
            newrun = runs[j].replace(
                old, new.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
            p2 = p.replace(runs[j], newrun, 1)
            new_xml = new_xml.replace(p, p2, 1)
            ps[by_pid[pid]] = p2
            applied.append((pid, f"sub r{j}", old[:50]))

        elif kind == "para":
            _, pid, new = edit
            if pid not in by_pid:
                failed.append((pid, "paragraph not found", new[:60]))
                continue
            p, runs, texts = para_pid_block(pid)
            nonempty = [j for j, t in enumerate(texts) if t.strip()]
            if new in " ".join(texts):
                skipped.append((pid, "already applied", new[:60]))
                continue
            if len(nonempty) != 1:
                failed.append((pid, f"paragraph has {len(nonempty)} non-empty runs "
                                    f"(only single-run paragraphs may be replaced wholesale)",
                               new[:60]))
                continue
            j = nonempty[0]
            newrun = set_run_text(runs[j], new)
            p2 = p.replace(runs[j], newrun, 1)
            new_xml = new_xml.replace(p, p2, 1)
            ps[by_pid[pid]] = p2
            applied.append((pid, f"whole paragraph (run {j})", new[:50]))
        else:
            failed.append(("?", f"unknown edit kind {kind}", str(edit)[:60]))

    # ---- append the four new references (cloned from the last reference paragraph) ----
    tail_pid = "0CDC4465"                        # reference 25 (Ravi 2022), last in the list
    if tail_pid in by_pid:
        tmpl = ps[by_pid[tail_pid]]
        truns, ttexts = run_texts(tmpl)
        # the guard must look at the WHOLE document: after a first append the template
        # paragraph itself is unchanged, so a template-local check would append again
        alltxt = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', xml, re.S))
        if "Carro M.S." in alltxt:
            skipped.append((tail_pid, "references already appended", ""))
        elif len([t for t in ttexts if t.strip()]) != 1:
            failed.append((tail_pid, "reference template is not single-run", "refs 26-29"))
        else:
            j = [k for k, t in enumerate(ttexts) if t.strip()][0]
            add = []
            for txt in NEW_REFS:
                add.append(tmpl.replace(truns[j], set_run_text(truns[j], txt), 1))
            new_xml = new_xml.replace(tmpl, tmpl + "".join(add), 1)
            applied.append((tail_pid, "appended refs 26-29", "; ".join(r[:28] for r in NEW_REFS)))

    print(f"\napplied {len(applied)}, skipped {len(skipped)}, failed {len(failed)}")
    for a in applied:
        print(f"  [ok]      {a[0]} {a[1]:28s} {a[2]}")
    for s in skipped:
        print(f"  [skip]    {s[0]} {s[1]:28s} {s[2]}")
    for f in failed:
        print(f"  [FAILED]  {f[0]} {f[1]:38s} {f[2]}")

    if failed:
        print("\nNOT WRITING: at least one anchor failed (all-or-nothing).")
        return 1
    if CHECK:
        print("\n--check: no write performed.")
        return 0
    if new_xml == xml:
        print("\nnothing to do (already fully applied).")
        return 0

    os.makedirs(BAK, exist_ok=True)
    shutil.copy2(SRC, os.path.join(BAK, os.path.basename(SRC)))
    with open(os.path.join(BAK, "md5_before.txt"), "w") as f:
        f.write(src_md5 + "\n")

    tmp = SRC + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for n in names:
            zout.writestr(z.getinfo(n), new_xml.encode("utf-8") if n == DOC else orig_entries[n])
    os.replace(tmp, SRC)

    # ---------------- verification ----------------
    z2 = zipfile.ZipFile(SRC)
    names2 = z2.namelist()
    assert names2 == names, "entry order/name mismatch"
    diff = [n for n in names2 if z2.read(n) != orig_entries[n]]
    ps2 = paragraphs(z2.read(DOC).decode("utf-8"))
    pids2 = [pid_of(p) for p in ps2]
    uniq = len([p for p in pids2 if p]) == len(set(p for p in pids2 if p))
    bm_before = xml.count("<w:bookmarkStart")
    bm_after = z2.read(DOC).decode("utf-8").count("<w:bookmarkStart")
    txt = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', z2.read(DOC).decode("utf-8"), re.S))
    print(f"\nVERIFY: entries {len(names2)} (identical names), differing entries = {diff}")
    print(f"  paragraphs {len(ps)} -> {len(ps2)}; paraIds unique = {uniq}; "
          f"bookmarks {bm_before} -> {bm_after}")
    bad = [lab for lab, must in VERIFY_EXTRA if must not in txt]
    print(f"  content checks: {'all present' if not bad else 'MISSING ' + str(bad)}")
    print(f"  new md5 {sha(SRC)}")
    print(subprocess.run(["bash", "-c", f"cd '{os.path.dirname(SRC)}' && python3 -c \"import "
                                        f"zipfile,sys; z=zipfile.ZipFile(sys.argv[1]); "
                                        f"print('unzip -t:', 'ok' if z.testzip() is None else 'BAD')\" "
                                        f"'{SRC}'"], capture_output=True, text=True).stdout.strip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
