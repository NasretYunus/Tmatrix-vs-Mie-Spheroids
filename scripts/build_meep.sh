#!/usr/bin/env bash
# Build Meep (FDTD) 1.29 from source with Python bindings, against Python >= 3.12.
#
# Why source: Meep is distributed through conda-forge, not PyPI. If conda is
# unavailable (or its channels are unreachable), the only route is a source build.
# All build dependencies are in the Ubuntu apt repositories; the Python bindings
# are generated with SWIG (--enable-maintainer-mode).
#
# Prerequisites (apt):
#   sudo apt-get install -y libctl-dev libhdf5-dev guile-3.0-dev libgsl-dev \
#        libharminv-dev libfftw3-dev liblapack-dev swig automake libtool \
#        pkg-config python3-dev
#   pip install numpy h5py mpi4py
#
# Run inside your activated virtualenv:  bash scripts/build_meep.sh
set -euo pipefail

MEEP_VER=1.29.0
WORK=$(mktemp -d)
echo "Building Meep $MEEP_VER in $WORK"
cd "$WORK"

# 1. Fetch source
wget -q "https://github.com/NanoComp/meep/archive/refs/tags/v${MEEP_VER}.tar.gz" -O meep.tar.gz
tar xzf meep.tar.gz
cd "meep-${MEEP_VER}"

# 2. Configure: serial (no MPI), no MPB, Python bindings on. The CPPFLAGS expose
#    Python.h and the numpy headers; maintainer-mode lets autotools run SWIG to
#    generate the wrapper (meep_wrap.cxx), which is not shipped in the tarball.
export PYTHON=$(which python)
export CPPFLAGS="-I$(python -c 'import sysconfig; print(sysconfig.get_path("include"))') -I$(python -c 'import numpy; print(numpy.get_include())')"
sh autogen.sh --without-mpi --without-mpb --with-libctl=/usr/share/libctl \
   --enable-shared --enable-maintainer-mode PYTHON="$PYTHON"

# 3. Compile (SWIG wrapper is large; -j speeds it up) and install
make -j"$(nproc)"
make install
ldconfig 2>/dev/null || true

# 4. Smoke test. Meep installs to the system site-packages and libmeep to
#    /usr/local/lib; expose both at run time.
export PYTHONPATH=/usr/local/lib/python3.12/site-packages:${PYTHONPATH:-}
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH:-}
python -c "import meep as mp; print('Meep', mp.__version__, 'OK')"
echo "Meep build complete."
echo "Remember to export PYTHONPATH and LD_LIBRARY_PATH (see above) before running FDTD scripts."
