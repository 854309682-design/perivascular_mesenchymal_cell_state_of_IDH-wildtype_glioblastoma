#!/usr/bin/env python3
"""rev21 — Final manuscript text pass: de-novo arm numbers, Figure 4 numbers, S1/S2 captions.

Applies the last round of reviewer-driven text edits as **paragraph-level zip surgery** on the
hand-polished .docx, exactly like rev11/rev20: the archive is rebuilt entry by entry and only
`word/document.xml` may differ. The manuscript is never rebuilt from markdown and never re-saved
with python-docx.

Edit kinds
    ("para", paraId, new_text)          whole paragraph text replacement; only allowed when the
                                        paragraph carries exactly one non-empty run
    ("run",  paraId, anchor, new_text)  replace the text of the *one* run containing `anchor`
                                        (used for the caption paragraphs, which have a bold
                                        label run plus a body run)
    ("sub",  paraId, old, new)          substring replacement inside a single run

Guarantees
    * dated backup written before anything else;
    * ALL-OR-NOTHING: if any anchor is missing nothing is written and every failure is listed;
    * idempotent: an edit whose new text is already present (whole-string, normalised
      comparison) is skipped, so the pass can be re-run after the de-novo arm or Figure 4 moves;
    * verification: identical entry names, only `word/document.xml` different, every other entry
      byte-identical, `<w:p`/paraId counts and uniqueness, bookmarks, `unzip -t`.

Run:
  /home/meng/miniconda3/bin/python3 scripts/rev21_manuscript_final.py [--check]
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
BAK = os.path.join(PROJ, "Manuscript", "backups", f"{date.today()}_pre_final_text_pass")
DOC = "word/document.xml"
CHECK = "--check" in sys.argv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rev21_edits import EDITS, VERIFY_EXTRA      # noqa: E402


def norm(t):
    """Normalise typographic variants so an "already applied" check cannot be fooled."""
    t = (t.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
          .replace("C/EBPbeta", "C/EBP\u03b2").replace("NF-kB", "NF-\u03baB"))
    t = re.sub(r"([x\u00d7])\s*10-(\d+)",
               lambda m: "\u00d7 10" + {"3": "\u207b\u00b3", "5": "\u207b\u2075",
                                         "7": "\u207b\u2077", "14": "\u207b\u00b9\u2074",
                                         "30": "\u207b\u00b3\u2070"}.get(m.group(2), "-" + m.group(2)), t)
    return re.sub(r"\s+", " ", t).strip()


def sha(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def para_re(pid):
    return re.compile(r'<w:p\b[^>]*w14:paraId="%s".*?</w:p>' % pid, re.S)


def runs_of(p):
    return re.findall(r"<w:r\b.*?</w:r>", p, re.S)


def text_of(r):
    return "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', r, re.S))


def set_run_text(run, new):
    """Replace a run's text, preserving its run properties (XML-escapes `new`)."""
    m = re.match(r"(<w:r\b[^>]*>)(.*?)(</w:r>)", run, re.S)
    open_tag, body, close = m.group(1), m.group(2), m.group(3)
    tprops = re.search(r"<w:rPr>.*?</w:rPr>", body, re.S)
    rpr = tprops.group(0) if tprops else ""
    return f'{open_tag}{rpr}<w:t xml:space="preserve">{esc(new)}</w:t>{close}'


def unesc(t):
    """Inverse of esc() for the entities that occur in this document (only &lt; in practice)."""
    return (t.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&"))


def main():
    src = sha(SRC)
    z = zipfile.ZipFile(SRC)
    names = z.namelist()
    ent = {n: z.read(n) for n in names}
    xml = ent[DOC].decode("utf-8")
    n_p = len(re.findall(r"<w:p\b", xml))
    print(f"source {SRC}\n  md5 {src}, {len(names)} entries, {n_p} paragraphs")

    new_xml = xml
    applied, skipped, failed = [], [], []

    for edit in EDITS:
        kind, pid = edit[0], edit[1]
        m = para_re(pid).search(new_xml)
        if not m:
            failed.append((pid, kind, "paragraph not found"))
            continue
        p = m.group(0)
        runs = runs_of(p)

        if kind == "para_sub":
            # apply several (old -> new) substring replacements inside ONE paragraph's single
            # non-empty run. Anchors may contain plain "<" (they are matched against the
            # unescaped text); everything not named in `pairs` is preserved verbatim. Used for
            # the long paragraphs (e.g. the abstract) where retyping risks silent drift.
            pairs = edit[2]
            nonempty = [i for i, r in enumerate(runs) if text_of(r).strip()]
            if len(nonempty) != 1:
                failed.append((pid, kind, f"{len(nonempty)} non-empty runs (needs exactly 1)"))
                continue
            j = nonempty[0]
            txt = unesc(text_of(runs[j]))
            new_text, bad = txt, []
            for old, rep in pairs:
                # check the replacement FIRST: for an append-style pair the anchor survives the
                # replacement, so an old-first test would append the clause on every re-run
                if norm(rep) in norm(new_text):
                    pass                             # already applied
                elif old in new_text:
                    new_text = new_text.replace(old, rep, 1)
                else:
                    bad.append(old[:44])
            if bad:
                failed.append((pid, kind, f"{len(bad)} anchor(s) missing: {bad}"))
                continue
            if new_text == txt:
                skipped.append((pid, kind, "already applied"))
                continue
            p2 = p.replace(runs[j], set_run_text(runs[j], new_text), 1)
        elif kind in ("para", "run"):
            new_text = edit[2] if kind == "para" else edit[3]
            if kind == "para":
                nonempty = [i for i, r in enumerate(runs) if text_of(r).strip()]
                if len(nonempty) != 1:
                    failed.append((pid, kind, f"{len(nonempty)} non-empty runs (needs exactly 1)"))
                    continue
                j = nonempty[0]
            else:
                anchor = edit[2]
                hits = [i for i, r in enumerate(runs) if anchor in text_of(r)]
                if not hits:
                    if any(norm(new_text) in norm(text_of(r)) for r in runs):
                        skipped.append((pid, kind, "already applied"))
                    else:
                        failed.append((pid, kind, f"run anchor not found: {anchor[:40]!r}"))
                    continue
                if len(hits) > 1:
                    failed.append((pid, kind, f"run anchor not unique ({len(hits)} runs)"))
                    continue
                j = hits[0]
            if norm(text_of(runs[j])) == norm(new_text):
                skipped.append((pid, kind, "already applied"))
                continue
            p2 = p.replace(runs[j], set_run_text(runs[j], new_text), 1)
        elif kind == "remove":
            # delete an exact substring from one run (used to repair duplicated text). Idempotent:
            # if the text is absent it is treated as already removed and skipped.
            _, _, text = edit
            hits = [i for i, r in enumerate(runs) if text in text_of(r)]
            if not hits:
                skipped.append((pid, kind, "text already absent"))
                continue
            if len(hits) > 1:
                failed.append((pid, kind, f"text to remove appears in {len(hits)} runs"))
                continue
            j = hits[0]
            if text_of(runs[j]).count(text) > 1:
                failed.append((pid, kind, "text to remove occurs more than once in the run"))
                continue
            p2 = p.replace(runs[j], runs[j].replace(text, "", 1), 1)
        elif kind == "sub":
            _, _, old, new_text = edit
            # check the NEW text FIRST: an append-style replacement (old is a substring of new,
            # e.g. "sentence." -> "sentence. Extra clause.") keeps its anchor after the first
            # application, so an old-first check would apply it again on every re-run
            if any(norm(new_text) in norm(text_of(r)) for r in runs):
                skipped.append((pid, kind, "already applied"))
                continue
            hits = [i for i, r in enumerate(runs) if old in text_of(r)]
            if not hits:
                failed.append((pid, kind, f"anchor not found: {old[:48]!r}"))
                continue
            if len(hits) > 1:
                failed.append((pid, kind, f"anchor not unique ({len(hits)} runs)"))
                continue
            j = hits[0]
            p2 = p.replace(runs[j], runs[j].replace(old, esc(new_text), 1), 1)
        else:
            failed.append((pid, kind, f"unknown edit kind {kind}"))
            continue

        if p2 == p:
            failed.append((pid, kind, "replacement produced no change"))
            continue
        new_xml = new_xml.replace(p, p2, 1)
        applied.append((pid, kind, str(edit[2])[:56]))

    print(f"\napplied {len(applied)}, skipped {len(skipped)}, failed {len(failed)}")
    for a in applied:
        print(f"  [ok]     {a[0]} {a[1]:5s} {a[2]!r}")
    for s in skipped:
        print(f"  [skip]   {s[0]} {s[1]:5s} {s[2]}")
    for f in failed:
        print(f"  [FAILED] {f[0]} {f[1]:5s} {f[2]}")

    if failed:
        print("\nNOT WRITING: at least one anchor failed (all-or-nothing).")
        return 1
    if CHECK:
        print("\n--check: no write performed.")
        return 0
    if new_xml == xml:
        print("\nnothing to do (already fully applied).")
        return 0
    try:
        import xml.etree.ElementTree as ET
        ET.fromstring(new_xml)
    except Exception as e:
        print(f"\nNOT WRITING: generated document.xml is not well-formed ({e}).")
        return 1

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
    x2 = z2.read(DOC).decode("utf-8")
    pids = re.findall(r'w14:paraId="([0-9A-Fa-f]+)"', x2)
    txt = "".join(re.findall(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', x2, re.S))
    print("\nVERIFY:")
    print(f"  entries identical: {z2.namelist() == names}; differing: "
          f"{[n for n in z2.namelist() if z2.read(n) != ent[n]]}")
    n_p2 = len(re.findall(r"<w:p\b", x2))
    print(f"  <w:p> {n_p} -> {n_p2}; paraIds {len(pids)}; unique {len(pids) == len(set(pids))}; "
          f"bookmarks {xml.count('<w:bookmarkStart')} -> {x2.count('<w:bookmarkStart')}")
    bad = [lab for lab, must in VERIFY_EXTRA if must not in txt]
    print(f"  content checks: {'all present' if not bad else 'MISSING ' + str(bad)}")
    print(f"  unzip -t: {'ok' if z2.testzip() is None else 'BAD'}")
    print(f"  new md5 {sha(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
