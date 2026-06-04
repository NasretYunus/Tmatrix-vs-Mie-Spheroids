#!/usr/bin/env bash
# Reproduce the headline results of the paper.
#
# Assumes pytmatrix and Meep are already built (see scripts/build_pytmatrix.sh
# and scripts/build_meep.sh) and the Python deps from requirements.txt are
# installed in the active environment.
#
# Usage:  bash reproduce.sh [quick|full|figs]
#   quick : validation + band-centre grid only (a few CPU-minutes)
#   full  : the complete (xi,kappa) grid at all three wavelengths (CPU-hours)
#   figs  : regenerate all figures from the shipped data/ (seconds; no T-matrix)
#
# Results land in build/.  The shipped data/ and figures/ are the reference.
set -euo pipefail
MODE="${1:-quick}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/build"
mkdir -p "$BUILD"
cd "$BUILD"

run() { echo ">>> $*"; "$@"; }

case "$MODE" in
  figs)
    # No solver needed: rebuild figures from the reference data.
    cp "$ROOT/data/full_grid.json" .
    run python "$ROOT/src/make_fig3_maps.py"
    echo "Figures written to $BUILD (compare against $ROOT/figures/)."
    ;;

  quick)
    # 1. Validate the T-matrix engine against Lorenz-Mie in the sphere limit.
    run python "$ROOT/src/validate_tmat.py"
    # 2. Compute the depolarization-factor shape correction.
    run python "$ROOT/src/shape_correction.py"
    # 3. Band-centre grid only: 550 nm, all 12 xi-rows (chunked).
    cp "$ROOT/data/full_grid.json" full_grid.json 2>/dev/null || true
    for xi0 in 0 3 6 9; do
      run python "$ROOT/src/full_grid.py" 550 "$xi0" 3
    done
    run python "$ROOT/src/make_fig3_maps.py"
    echo "Quick reproduction done. Compare $BUILD against $ROOT/data and $ROOT/figures."
    ;;

  full)
    # The complete grid: three wavelengths x twelve xi-rows, chunked so each
    # invocation finishes quickly and persists to full_grid.json.
    : > full_grid.json
    for lam in 400 550 800; do
      for xi0 in 0 2 4 6 8 10; do
        run python "$ROOT/src/full_grid.py" "$lam" "$xi0" 2
      done
    done
    run python "$ROOT/src/make_fig3_maps.py"
    echo "Full grid done (note: the xi=5/400nm row is beyond reliable EBCM and is"
    echo "flagged ebcm_limited in full_grid.json, as discussed in the paper)."
    ;;

  *)
    echo "Unknown mode '$MODE'. Use: quick | full | figs"; exit 1;;
esac
