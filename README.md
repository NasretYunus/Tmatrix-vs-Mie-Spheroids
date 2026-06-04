#  Computational pipeline for the paper *When shape and alignment conspire: polarization
in spheroidal zeolite films

The study compares a rigorous single-spheroid **T-matrix** computation against a fast
**equivalent-sphere Mie** model (carrying a quasi-static depolarization-factor shape
correction) for zeolite inclusions (n_p = 1.47) in a PMMA host (n_h = 1.50), across
aspect ratio ξ ∈ [0.2, 5.0] and orientation concentration κ ∈ [0, 50] (von Mises–Fisher),
at three visible wavelengths. An independent **FDTD** (Meep) cross-check anchors the
rigorous result at both a sphere and the extreme prolate needle.

**Headline result.** For this near-index-matched system the two methods agree to within
~1% across **87% of the (ξ,κ) plane**. They diverge only at the aligned-needle corner
(high ξ and high κ together), where peak degree of linear polarization collapses from
99.5% to ~58% at band centre and the symmetry-breaking metric rises three orders of
magnitude. The practical consequence is a permissive design rule: the cheap model is
safe almost everywhere; rigor is mandatory only for aligned high-aspect-ratio fillers.

## Repository layout

```
src/                 analysis and simulation scripts (relative paths; run from build/)
  validate_tmat.py     T-matrix vs Lorenz-Mie sphere-limit validation
  mie_baseline.py      Lorenz-Mie baseline (Figs. 1-2 numbers)
  full_grid.py         resumable (xi,kappa,lambda) T-matrix grid generator
  shape_correction.py  quasi-static depolarization-factor shape model
  make_fig3_maps.py    builds the (xi,kappa) contour maps (Fig. 3)
  make_all_figs.py     regenerates all four figures in PNG (1000 dpi), PDF, SVG, EPS
  ff_run.py            one cached Meep far-field run (sphere or spheroid)
  combine_dolp.py      assembles DOLP from cached far-field runs; FDTD vs reference
data/
  full_grid.json       reference grid: 525 computed cells across 3 wavelengths
  shape_corr.json      shape-corrected peak-DOLP vs aspect ratio
figures/             the four figures as published, each in PNG (1000 dpi), PDF, SVG, EPS
manuscript/          manuscript in Markdown, LaTeX, compiled PDF, and Word
scripts/
  build_pytmatrix.sh   build pytmatrix from source (Python >= 3.12 safe)
  build_meep.sh        build Meep 1.29 (FDTD) from source with Python bindings
reproduce.sh         one-command reproduction (quick | full | figs)
requirements.txt     pip dependencies
LICENSE              MIT (plus notes on third-party solver licenses)
```

## Environment setup

Tested on Ubuntu 24.04 with Python 3.12.

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

Two dependencies are **not** on PyPI and are built from source:

- **pytmatrix** wraps Mishchenko's `ampld.lp.f` Fortran T-matrix code. The PyPI build
  relies on `numpy.distutils`, removed in Python 3.12, so it is built with the meson
  f2py backend instead:
  ```bash
  sudo apt-get install -y gfortran
  pip install meson ninja
  bash scripts/build_pytmatrix.sh
  ```

- **Meep** (FDTD) ships through conda-forge, not PyPI. If conda is unavailable, build
  from source; every dependency is in apt:
  ```bash
  sudo apt-get install -y libctl-dev libhdf5-dev guile-3.0-dev libgsl-dev \
       libharminv-dev libfftw3-dev liblapack-dev swig automake libtool \
       pkg-config python3-dev
  bash scripts/build_meep.sh
  ```
  Meep installs `libmeep` to `/usr/local/lib` and the module to the system
  site-packages; export both before running FDTD scripts:
  ```bash
  export PYTHONPATH=/usr/local/lib/python3.12/site-packages:$PYTHONPATH
  export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
  ```

The T-matrix and Mie parts (Figs. 1–3, Table 1, the whole grid) need **only** pytmatrix;
Meep is required solely for the FDTD cross-check (Fig. 4).

## Reproducing the results

```bash
bash reproduce.sh figs     # rebuild figures from shipped data (seconds, no solver)
bash reproduce.sh quick    # validation + band-centre grid (a few CPU-minutes)
bash reproduce.sh full     # complete grid, all wavelengths (CPU-hours)
```

Outputs land in `build/`. The shipped `data/` and `figures/` are the reference to
compare against.

### FDTD cross-check (Fig. 4)

The FDTD comparison uses a two-incident-polarization scattered-field calculation with
empty-run subtraction. Each panel needs four cached runs (perpendicular/parallel ×
empty/particle); `ff_run.py` does one run and caches it, `combine_dolp.py` assembles
the DOLP curve. For the sphere validation at reduced size parameter:

```bash
cd build
for pol in y x; do for part in 0 1; do
  python ../src/ff_run.py sphere 1.0 $pol $part 20
done; done
python ../src/combine_dolp.py sphere 1.0 20    # prints FDTD peak vs Mie
```

For the extreme prolate needle, broadside, at the full band-centre size parameter
x = 8.57 (heavier; lower resolution keeps it tractable):

```bash
for pol in y x; do for part in 0 1; do
  python ../src/ff_run.py prolate 5.0 $pol $part 10 8.57 90
done; done
```

## Validation status

Three independent checks pass:

1. **Sphere limit** — at axis ratio → 1 the T-matrix reproduces the independent
   Lorenz–Mie peak DOLP (99.5%) and scattering efficiency to all reported digits, with
   the symmetry-breaking metric B → 0.
2. **Isotropy recovery** — a spheroid at random orientation gives B ≈ 1e-4, as required.
3. **FDTD** — peak DOLP 98.5% (FDTD) vs 96.2% (Mie) for the sphere; 99.2% (FDTD) vs
   99.9% (T-matrix) for the broadside needle at x = 8.57.

## Known scope limits (not bugs)

- The extreme-prolate row at 400 nm (ξ = 5, the largest size parameter) lies beyond the
  reliable convergence range of the null-field/EBCM method; those 15 cells are flagged
  `ebcm_limited` in `full_grid.json` rather than reported.
- The cliff's ~58% collapse is an **orientation-averaging** effect, not a single-particle
  one. A full orientation-averaged FDTD reconstruction would require integrating dozens of
  fixed-orientation runs over the von Mises–Fisher distribution; the per-orientation FDTD
  agreement (above) is what the present cross-check establishes.

## License

MIT (see `LICENSE`). Third-party solvers retain their own licenses: pytmatrix (MIT)
over Mishchenko's public-domain Fortran, miepython (MIT), Meep (GPL v2+).
