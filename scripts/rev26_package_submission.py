#!/usr/bin/env python3
"""rev26 — package the submission: sync figures, render, merge, rebuild the archives.

One command, in the order the project's handoff documents it, so the shipped package cannot
drift from `figures/`:

  1. copy the seven rebuilt figures from `figures/` to `Manuscript/` (both formats) — the merge
     below reads `Manuscript/`, and two of them were still the September versions;
  2. render the manuscript .docx to PDF with LibreOffice;
  3. `pdfunite` the manuscript + Figures 1-5 + Figures S1-S2, then rewrite through `qpdf`
     (never wrapped in `set -e`: qpdf exits 3 "succeeded with warnings" and `pdfunite` output
     carries stale object numbers);
  4. validate the merged PDF with `qpdf --check` **in /tmp** (`qpdf --check` cannot open files
     on the `/run/media` mount) and report its page count / md5;
  5. copy the merged PDF into the 备案 folder and run the archive builder (`rev25`).

Run:
  GLIOMA_PROJ=/run/media/meng/Seqdata/glioma/Paper_2 \
  /home/meng/miniconda3/bin/python3 scripts/rev26_package_submission.py
"""
import glob
import hashlib
import os
import shutil
import subprocess
import sys

PROJ = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(PROJ, "Manuscript")
FIG = os.path.join(PROJ, "figures")
BEIAN = glob.glob(os.path.join(PROJ, "*备案*"))[0]
DOCX = os.path.join(MAN, "Manuscript_perivascular_mesenchymal_cell_in_glioma.docx")
FIGS = ["Figure1_main", "Figure2_main", "Figure3_main", "Figure4_main", "Figure5_main",
        "FigureS1_supp", "FigureS2_supp"]


def md5p(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    return r.returncode, (r.stdout + r.stderr).strip()


def main():
    # 1. sync the figures the merge reads
    print("1. syncing figures -> Manuscript/")
    for stem in FIGS:
        for ext in ("pdf", "jpeg"):
            src = os.path.join(FIG, f"{stem}.{ext}")
            dst = os.path.join(MAN, f"{stem}.{ext}")
            if not os.path.exists(src):
                # Figures that this revision did not regenerate (Figure 1) live only in
                # Manuscript/; keep the shipped copy rather than failing the package
                if os.path.exists(dst):
                    print(f"   {stem}.{ext}  unchanged (not regenerated; {os.path.getsize(dst)/1e3:.0f} kB)")
                    continue
                print(f"   MISSING {src} and no shipped copy at {dst}")
                return 1
            shutil.copy2(src, dst)
            print(f"   {stem}.{ext}  {os.path.getsize(src)/1e3:.0f} kB  md5 {md5p(src)}")

    # 2. render
    print("2. rendering the manuscript")
    rc, out = run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", MAN, DOCX])
    pdf = DOCX[:-5] + ".pdf"
    if not os.path.exists(pdf):
        print(f"   render FAILED: {out}")
        return 1
    print(f"   {os.path.basename(pdf)}  {os.path.getsize(pdf)/1e3:.0f} kB  md5 {md5p(pdf)}")

    # 3. merge
    print("3. merging")
    parts = [pdf] + [os.path.join(MAN, f"{f}.pdf") for f in FIGS]
    raw = "/tmp/paper2_raw.pdf"
    merged = "/tmp/paper2_all.pdf"
    rc, out = run(["pdfunite"] + parts + [raw])
    if rc != 0:
        print(f"   pdfunite FAILED: {out}")
        return 1
    rc, out = run(["qpdf", raw, merged])           # exit 3 = "succeeded with warnings"
    print(f"   qpdf rewrite: rc={rc} {out.splitlines()[0] if out else ''}")

    # 4. validate in /tmp
    rc, out = run(["qpdf", "--check", merged])
    if "No syntax or stream encoding errors found" not in out:
        print(f"   qpdf --check reported a problem:\n{out}")
        return 1
    rc, out = run(["pdfinfo", merged])
    pages = [l for l in out.splitlines() if l.startswith("Pages")]
    print(f"   merged: {pages}, {os.path.getsize(merged)/1e6:.2f} MB, md5 {md5p(merged)}")

    # 5. ship + archive
    dst = os.path.join(BEIAN, "Paper_2_all.pdf")
    shutil.copy2(merged, dst)
    print(f"5. {os.path.relpath(dst, PROJ)}  md5 {md5p(dst)}")
    rc, out = run([sys.executable, os.path.join(PROJ, "scripts", "rev25_build_submission_archive.py")],
                  env={**os.environ, "GLIOMA_PROJ": PROJ})
    print(out)
    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
