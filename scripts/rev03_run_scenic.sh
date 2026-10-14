#!/usr/bin/env bash
# rev03 — de-novo regulatory-network inference (pySCENIC 0.12.1) for the four cohorts.
#
# Reviewer item M2/DA2/DA3: the submitted prior-independent ("de-novo") arm behind
# Fig. 3b and Results §3.3 was produced by code that no longer exists. This driver
# rebuilds it from the looms of rev01, with the standard three-step pySCENIC chain and
# every parameter recorded here:
#
#   1. GRNBoost2 co-expression modules      (arboreto, seed 42)
#   2. cisTarget motif pruning to regulons  (hg38 refseq-r80 mc9nr, 10kb + 500bp)
#   3. AUCell per-cell regulon activity
#
# Environment: the Python 3.11 environment that carries pySCENIC 0.12.1 + ctxcore +
# arboreto (Paper_1/env/.venv, interpreter restored 2026-10-14 with `uv python install 3.11`).
#
# Usage:  bash scripts/rev03_run_scenic.sh [cohort ...]      (default: all four)
set -u

PY=${PY:-/run/media/meng/Seqdata/glioma/Paper_1/env/.venv/bin/python}
PROJ=${GLIOMA_PROJ:-/run/media/meng/Seqdata/glioma/Paper_2}
DATA=${GLIOMA_DATA:-/run/media/meng/Seqdata/glioma/Datasets}
RES="$DATA/scenic_resources"
SCENIC="$PROJ/results/scenic"
LOGS="$PROJ/logs"
SHIM="$PROJ/scripts/rev03_scenic_cli.py"
NWORK=12

TFS="$RES/allTFs_hg38.txt"
DB10="$RES/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
DB500="$RES/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIFS="$RES/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

COHORTS=${*:-"GSE103224 GSE182109 Neftel_GSE131928"}
mkdir -p "$SCENIC" "$LOGS"

for C in $COHORTS; do
  LOOM="$SCENIC/rev01_${C}.loom"
  ADJ="$SCENIC/rev03_${C}_adj.tsv"
  REG="$SCENIC/rev03_${C}_regulons.csv"
  AUC="$SCENIC/rev03_${C}_aucell.csv"
  [ -f "$LOOM" ] || { echo "MISSING $LOOM"; continue; }

  if [ ! -s "$ADJ" ]; then
    echo "=== [$C] GRNBoost2 ==="
    "$PY" "$SHIM" grn "$LOOM" "$TFS" --method grnboost2 --output "$ADJ" \
        --num_workers "$NWORK" --seed 42 >> "$LOGS/rev03_${C}.log" 2>&1 \
        || { echo "[$C] GRN failed"; continue; }
  fi

  if [ ! -s "$REG" ]; then
    echo "=== [$C] cisTarget ctx ==="
    "$PY" "$SHIM" ctx "$ADJ" "$DB10" "$DB500" --annotations_fname "$MOTIFS" \
        --expression_mtx_fname "$LOOM" --output "$REG" \
        --num_workers "$NWORK" --nes_threshold 3.0 --rank_threshold 5000 \
        >> "$LOGS/rev03_${C}.log" 2>&1 \
        || { echo "[$C] ctx failed"; continue; }
  fi

  if [ ! -s "$AUC" ]; then
    echo "=== [$C] AUCell ==="
    "$PY" "$SHIM" aucell "$LOOM" "$REG" --output "$AUC" \
        --num_workers "$NWORK" --seed 42 >> "$LOGS/rev03_${C}.log" 2>&1 \
        || { echo "[$C] aucell failed"; continue; }
  fi
  echo "[$C] done: $(wc -l < "$ADJ") adjacencies, $(wc -l < "$REG") motif rows, $(wc -l < "$AUC") aucell rows"
done
echo "ALL DONE"
