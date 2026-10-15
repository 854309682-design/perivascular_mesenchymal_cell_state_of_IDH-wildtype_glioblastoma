#!/usr/bin/env python3
"""Renumber the manuscript's references into order of first citation (Vancouver).

The reference list is a Word auto-numbered list (no literal digits in the text), so the
displayed numbers follow the list order. This script therefore:

  1. walks the body paragraphs in document order and records, for every in-text marker
     ``[n]`` / ``[n,m]`` / ``[a-b]``, the reference numbers in the order they are *read*
     (paragraph, then position in the paragraph, then position inside the marker);
  2. derives the old -> new permutation from the first appearance of each reference;
  3. rewrites every in-text marker with the new numbers (sorted ascending, runs of >= 3
     collapsed to ``a-b`` to match the manuscript's existing style);
  4. permutes the 29 reference entries so that list position i holds the entry that is
     now number i (each entry is a single run, so only its text moves - formatting,
     pPr and numPr stay with the paragraph).

Nothing else in ``word/document.xml`` is touched: no paragraph is moved, no number of a
reference entry is typed by hand, and no reference text is edited.

Usage:
    python rev29_reference_order.py            # dry run: print the permutation + checks
    python rev29_reference_order.py --apply    # apply, then verify
"""
from __future__ import annotations

import datetime as _dt
import os
import re
import shutil
import sys
import zipfile
from xml.sax.saxutils import escape, unescape

DOCX = "Manuscript/Manuscript_perivascular_mesenchymal_cell_in_glioma.docx"
PARA_RE = re.compile(r'<w:p\b[^>]*?w14:paraId="([0-9A-Fa-f]+)".*?</w:p>', re.S)
RUN_RE = re.compile(r"<w:r\b.*?</w:r>", re.S)
TEXT_RE = re.compile(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", re.S)
MARKER_RE = re.compile(r"\[([0-9][0-9,\u2013\u2014\- ]*)\]")


def para_text(raw: str) -> str:
    return unescape("".join(TEXT_RE.findall(raw)))


def is_numbered_ref(raw: str) -> bool:
    return bool(re.search(r"<w:numPr>", raw))


def set_run_text(raw: str, new_text: str) -> str:
    """Replace the text of the paragraph's single run, keeping every other byte."""
    runs = list(RUN_RE.finditer(raw))
    if len(runs) != 1:
        raise ValueError(f"expected a single run, found {len(runs)}")
    r = runs[0]
    body = r.group(0)
    tm = TEXT_RE.search(body)
    if tm is None:
        raise ValueError("run has no <w:t>")
    fixed = body[: tm.start(1)] + escape(new_text) + body[tm.end(1):]
    return raw[: r.start()] + fixed + raw[r.end():]


def renumber_paragraph(raw: str, new_of: dict[int, int]) -> tuple[str, int]:
    """Renumber every in-text marker of one paragraph, run by run.

    Markers live inside single runs in this manuscript. If one ever spans two runs it is
    left untouched here, which the caller's verification then catches (a marker still
    carrying an old number breaks the cited-entry identity check), so the write is
    refused rather than half-applied.
    """
    out = raw
    n = 0
    for r in list(RUN_RE.finditer(raw)):
        body = r.group(0)
        tm = TEXT_RE.search(body)
        if tm is None:
            continue
        t = unescape(tm.group(1))

        def sub(m: re.Match) -> str:
            nonlocal n
            n += 1
            return render_marker([new_of[x] for x in expand_marker(m.group(1))])

        t2 = MARKER_RE.sub(sub, t)
        if t2 == t:
            n -= len(MARKER_RE.findall(t))      # processed but unchanged
        if t2 != t:
            fixed = body[: tm.start(1)] + escape(t2) + body[tm.end(1):]
            out = out.replace(body, fixed, 1)
    return out, n


def expand_marker(inner: str) -> list[int]:
    nums: list[int] = []
    for part in inner.split(","):
        part = part.strip()
        if re.fullmatch(r"\d+", part):
            nums.append(int(part))
        elif re.fullmatch(r"\d+\s*[\u2013\u2014-]\s*\d+", part):
            a, b = (int(x) for x in re.findall(r"\d+", part))
            nums.extend(range(a, b + 1))
        else:
            raise ValueError(f"unparsable marker part {part!r}")
    return nums


def render_marker(new_nums: list[int]) -> str:
    """Ascending, runs of >= 3 collapsed to 'a-b' (the manuscript's existing style)."""
    ns = sorted(set(new_nums))
    out: list[str] = []
    i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        run = ns[i: j + 1]
        if len(run) >= 3:
            out.append(f"{run[0]}-{run[-1]}")
        else:
            out.extend(str(n) for n in run)
        i = j + 1
    return "[" + ",".join(out) + "]"


def load():
    zf = zipfile.ZipFile(DOCX)
    xml = zf.read("word/document.xml").decode("utf-8")
    names = zf.namelist()
    zf.close()
    return xml, names


def parse(xml: str):
    return [(m.group(1), m.start(), m.end(), m.group(0)) for m in PARA_RE.finditer(xml)]


def main() -> int:
    apply = "--apply" in sys.argv
    xml, names = load()
    paras = parse(xml)

    head_idx = next(
        i for i, p in enumerate(paras)
        if para_text(p[3]).strip().lower() == "references"
    )
    # reference list = the run of auto-numbered paragraphs straight after the heading
    ref_idx: list[int] = []
    for i in range(head_idx + 1, len(paras)):
        if is_numbered_ref(paras[i][3]):
            ref_idx.append(i)
        elif ref_idx:
            break
    if len(ref_idx) != 29:
        print(f"!! expected 29 numbered reference paragraphs, found {len(ref_idx)}")
        return 1

    ref_pid = {n: paras[i][0] for n, i in enumerate(ref_idx, start=1)}
    ref_text = {n: para_text(paras[i][3]) for n, i in enumerate(ref_idx, start=1)}

    # ---- citation order: (paragraph, marker position, position inside the marker) ----
    body_idx = [i for i in range(len(paras)) if i not in ref_idx]
    events: list[tuple[int, int, int, int]] = []   # (para, pos, k-th number, ref)
    for i in body_idx:
        t = para_text(paras[i][3])
        for m in MARKER_RE.finditer(t):
            for k, n in enumerate(expand_marker(m.group(1))):
                events.append((i, m.start(), k, n))
    cited = sorted({n for *_x, n in events})
    if cited != list(range(1, 30)):
        print(f"!! references cited are {cited}\n   expected 1..29 with no gaps")
        return 1

    order: list[int] = []
    for _i, _p, _k, n in sorted(events):
        if n not in order:
            order.append(n)
    if len(order) != 29:
        print(f"!! first-appearance order has {len(order)} entries, expected 29")
        return 1
    new_of = {old: new for new, old in enumerate(order, start=1)}

    print("First-citation order and the resulting renumbering\n")
    print(f"  {'new':>4} {'old':>4}   entry")
    for new, old in enumerate(order, start=1):
        print(f"  {new:>4} {old:>4}   {ref_text[old][:84]}")

    # ---- rewrite markers ----
    new_xml = xml
    changed_markers = 0
    for pid, s, e, raw in paras:
        if is_numbered_ref(raw):
            continue
        if not MARKER_RE.search(para_text(raw)):
            continue
        new_raw, n = renumber_paragraph(raw, new_of)
        changed_markers += n
        if new_raw != raw:
            new_xml = new_xml.replace(raw, new_raw, 1)

    # ---- permute the reference entries ----
    for new, old in enumerate(order, start=1):
        raw = next(p[3] for p in paras if p[0] == ref_pid[new])
        new_xml = new_xml.replace(raw, set_run_text(raw, ref_text[old]), 1)

    # ---- verification on the rewritten XML ----
    xp = parse(new_xml)
    xref_idx: list[int] = []
    for i in range(head_idx + 1, len(xp)):
        if is_numbered_ref(xp[i][3]):
            xref_idx.append(i)
        elif xref_idx:
            break
    new_ref_text = {n: para_text(xp[i][3]) for n, i in enumerate(xref_idx, start=1)}
    ok = True
    if sorted(new_ref_text.values()) != sorted(ref_text.values()):
        print("!! reference entry texts are not a permutation of the originals"); ok = False
    if len(xref_idx) != 29 or not all(is_numbered_ref(xp[i][3]) for i in xref_idx):
        print("!! reference list length/numbering changed"); ok = False

    # identity of every marker: same set of cited *entries* before and after
    def marker_ids(pxml) -> list[frozenset[str]]:
        out = []
        for i in range(len(pxml)):
            if is_numbered_ref(pxml[i][3]):
                continue
            for m in MARKER_RE.finditer(para_text(pxml[i][3])):
                out.append(frozenset(ref_text[x] if pxml is paras else new_ref_text[x]
                                     for x in expand_marker(m.group(1))))
        return out
    before_ids = marker_ids(paras)
    after_ids = marker_ids(xp)
    if before_ids != after_ids:
        print("!! a marker no longer cites the same articles"); ok = False
    if len(after_ids) != len(before_ids):
        print(f"!! marker count changed {len(before_ids)} -> {len(after_ids)}"); ok = False

    # the requirement itself: first appearance must now be 1,2,3,...
    ev2 = []
    for i in range(len(xp)):
        if is_numbered_ref(xp[i][3]):
            continue
        for m in MARKER_RE.finditer(para_text(xp[i][3])):
            for k, n in enumerate(expand_marker(m.group(1))):
                ev2.append((i, m.start(), k, n))
    seq: list[int] = []
    for _i, _p, _k, n in sorted(ev2):
        if n not in seq:
            seq.append(n)
    if seq != list(range(1, 30)):
        print(f"!! first-appearance order is not 1..29: {seq}"); ok = False

    print(f"\n  markers rewritten      : {changed_markers}")
    print(f"  reference entries moved: {sum(1 for new, old in enumerate(order, start=1) if new != old)}")
    print(f"  first-appearance order : {seq}")
    print(f"  checks                 : {'ALL PASS' if ok else 'FAILED'}")

    if not apply:
        print("\nDRY RUN - nothing written (pass --apply)")
        return 0 if ok else 1
    if not ok:
        print("\nNOT WRITING: verification failed")
        return 1

    # ---- write (zip surgery: only word/document.xml differs) ----
    today = _dt.date.today().isoformat()
    bak = os.path.join(os.path.dirname(DOCX), "backups", f"{today}_pre_reference_order")
    os.makedirs(bak, exist_ok=True)
    shutil.copy2(DOCX, os.path.join(bak, os.path.basename(DOCX)))
    tmp = DOCX + ".tmp"
    with zipfile.ZipFile(DOCX) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                data = new_xml.encode("utf-8")
            zout.writestr(item, data)
    os.replace(tmp, DOCX)
    print(f"\napplied. backup: {bak}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
