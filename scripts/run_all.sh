#!/usr/bin/env bash
# Paper 2 pipeline — run the GSE103224 single-cohort analysis in order.
# Usage: bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PY:-python3}
echo "[run_all] Step 0: build the GSE103224 raw object (skipped if it already exists)"
echo "[run_all]         environment/accessions: see env/ and tables/datasets.tsv"
$PY scripts/00_build_gse103224_raw.py
echo "[run_all] Step 1: QC"
$PY scripts/01_qc.py
echo "[run_all] Step 2: integration + inferCNV"
$PY scripts/02_integration_infercnv.py
echo "[run_all] Step 3: four-state scoring"
$PY scripts/03_state_scoring_shift.py
echo "[run_all] Step 4: MES transition trajectory"
$PY scripts/04_trajectory.py
echo "[run_all] Step 5: regulon / TF activity"
$PY scripts/05_regulon.py
echo "[run_all] Done. Review results/step01_qc..step05_regulon."
