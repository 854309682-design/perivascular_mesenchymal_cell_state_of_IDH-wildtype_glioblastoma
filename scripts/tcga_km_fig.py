#!/usr/bin/env python3
"""Standalone KM figure from the saved tcga_survival.tsv (IDH-stratified)."""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

OUT = os.path.join(lib.PROJ, "results", "step07_clinical")
d = pd.read_csv(os.path.join(OUT, "tcga_survival.tsv"), sep="\t")


def km(times, events):
    o = np.argsort(times); t = np.asarray(times)[o]; e = np.asarray(events)[o]
    n = len(t); s = 1.0; xs = [0.0]; ys = [1.0]
    for i in range(n):
        if e[i] == 1:
            s *= (1 - 1.0 / (n - i)); xs.append(float(t[i])); ys.append(s)
    return xs, ys


fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
for col, idhv, label in [(0, 1, "IDH-mutant"), (1, 0, "IDH-wildtype")]:
    sub = d[d["IDH_mutated"] == idhv]
    ax[col].set_title(f"{label} (n={len(sub)})")
    for grp, c in [("high", "#C0392B"), ("low", "#2980B9")]:
        g = sub[sub["MES_group"] == grp]
        if len(g) < 3:
            continue
        xs, ys = km(g["os_time"].values, g["os_event"].values)
        ax[col].step(xs, ys, where="post", c=c, label=f"MES {grp} (n={len(g)})")
    ax[col].set_xlabel("days"); ax[col].set_ylabel("overall survival"); ax[col].legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "tcga_km.png"), dpi=150)
print("KM figure saved", os.path.getsize(os.path.join(OUT, "tcga_km.png")), "bytes")
