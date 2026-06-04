#!/usr/bin/env bash
# Build pytmatrix from source against a Python >= 3.12 environment.
#
# Why source: pytmatrix wraps Mishchenko's Fortran T-matrix code (ampld.lp.f) via
# f2py. The PyPI build assumes numpy.distutils, which was removed in Python 3.12,
# so the wheel install fails on modern Python. We build the Fortran extension with
# the meson f2py backend instead, and patch a deprecated SciPy call.
#
# Prerequisites (apt): gfortran, and a working numpy/scipy in the active venv.
#   sudo apt-get install -y gfortran
#   pip install numpy scipy meson ninja
#
# Run inside your activated virtualenv:  bash scripts/build_pytmatrix.sh
set -euo pipefail

WORK=$(mktemp -d)
echo "Building pytmatrix in $WORK"
cd "$WORK"

# 1. Fetch source
pip download --no-deps --no-binary :all: pytmatrix -d .
tar xzf pytmatrix-*.tar.gz
cd pytmatrix-*/

# 2. Inline ampld.par.f into the Fortran source so f2py sees the parameters,
#    then build the extension module with the meson backend (Python 3.12-safe).
PYF_DIR=pytmatrix/fortran_tm
if [ -f "$PYF_DIR/ampld.par.f" ]; then
  # ampld.lp.f includes ampld.par.f via an INCLUDE; f2py+meson needs it inlined.
  python - <<'PY'
import re, pathlib
d = pathlib.Path("pytmatrix/fortran_tm")
src = (d/"ampld.lp.f").read_text()
par = (d/"ampld.par.f").read_text()
src = src.replace("      INCLUDE 'ampld.par.f'\n", par)
(d/"ampld.lp.f").write_text(src)
print("inlined ampld.par.f")
PY
fi

# 3. Patch deprecated scipy.integrate.trapz -> trapezoid (removed in SciPy >=1.14)
grep -rl "trapz" pytmatrix/*.py | while read -r f; do
  sed -i 's/\btrapz\b/trapezoid/g' "$f" || true
done

# 4. Compile the Fortran extension and install
f2py -c -m pytmatrix.fortran_tm.pytmatrix \
     "$PYF_DIR/ampld.lp.f" --backend meson 2>/dev/null || \
  echo "NOTE: if the f2py line fails, build via 'pip install . --no-build-isolation' after the patches above"
pip install . --no-build-isolation

# 5. Smoke test
python - <<'PY'
from pytmatrix.tmatrix import Scatterer
s = Scatterer(radius=0.5, wavelength=1.0, m=complex(1.5,0), axis_ratio=1.0001)
print("pytmatrix OK; Z[0,0] at default geometry =", s.get_Z()[0,0])
PY
echo "pytmatrix build complete."
