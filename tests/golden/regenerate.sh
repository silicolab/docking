#!/usr/bin/env bash
#
# Regenerate the golden ground-truth values used to validate the `docking` crate
# against official AutoDock Vina.
#
# These golden files are committed to the repo. Regenerate ONLY when intentionally
# updating the ground truth (e.g. bumping the reference Vina version), and review the
# diff carefully — the whole point of the harness is that Rust output matches THESE
# numbers, so changing them silently defeats the test.
#
#   Reference engine : AutoDock Vina v1.2.7 (official precompiled Windows binary)
#   Determinism      : verified bit-for-bit reproducible with a fixed --seed and --cpu 1
#   Inputs           : prepared receptor/ligand PDBQT taken verbatim from Vina's own
#                      example/basic_docking/solution (system "1iep").
#
# Usage (from repo root):  bash tests/golden/regenerate.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VINA="${VINA:-vina}"

if [[ ! -x "$VINA" && ! -f "$VINA" ]]; then
  echo "error: official vina binary not found at: $VINA" >&2
  echo "set \$VINA to the official AutoDock Vina v1.2.7 binary." >&2
  exit 1
fi
echo "Using: $("$VINA" --version | head -1)"

SEED=42
EXH=8
NUM_MODES=9

regen_system() {
  local G="$HERE/$1"
  local COMMON=(--receptor "$G/receptor.pdbqt" --ligand "$G/ligand.pdbqt" --config "$G/config.txt" --cpu 1)
  echo ">>> $1: score_only"
  "$VINA" "${COMMON[@]}" --score_only > "$G/score_only.out.txt" 2>&1
  echo ">>> $1: local_only (BFGS from input pose)"
  "$VINA" "${COMMON[@]}" --local_only --out "$G/local_only.pdbqt" > "$G/local_only.out.txt" 2>&1
  echo ">>> $1: full dock (seed=$SEED exhaustiveness=$EXH)"
  "$VINA" "${COMMON[@]}" --seed "$SEED" --exhaustiveness "$EXH" --num_modes "$NUM_MODES" \
      --out "$G/dock_seed${SEED}.pdbqt" > "$G/dock_seed${SEED}.out.txt" 2>&1
}

regen_system 1iep

echo "done."
