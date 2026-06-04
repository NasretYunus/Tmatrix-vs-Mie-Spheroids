"""Validate the compiled Mishchenko T-matrix engine against the independent
miepython Lorenz-Mie baseline at the band centre (lambda0 = 550 nm).

At axis_ratio -> 1 the T-matrix MUST reproduce the sphere: same peak DOLP,
same scattering efficiency, and the symmetry-breaking metric B -> 0.
"""
import numpy as np
from pytmatrix.tmatrix import Scatterer
from pytmatrix import orientation, scatter
import miepython

# --- physical system -------------------------------------------------------
n_h   = 1.50               # PMMA host
n_p   = 1.47               # zeolite inclusion
m_rel = n_p / n_h          # 0.980 relative index (particle rarer than host)
d_eqv = 1.0                # volume-equivalent diameter, micron
a_eqv = d_eqv / 2.0        # equal-volume-sphere radius
lam0  = 0.550              # vacuum wavelength, micron
lam_m = lam0 / n_h         # wavelength in the host medium

# --- independent miepython baseline ----------------------------------------
x_host = 2*np.pi*a_eqv/lam_m
theta = np.linspace(0, np.pi, 361)
mu = np.cos(theta)
qext, qsca, qback, g = miepython.efficiencies_mx(complex(m_rel, 0.0), x_host)
S1, S2 = miepython.S1_S2(complex(m_rel, 0.0), x_host, mu)
S11_mie = 0.5*(np.abs(S1)**2 + np.abs(S2)**2)
S12_mie = 0.5*(np.abs(S2)**2 - np.abs(S1)**2)
dolp_mie = -S12_mie / S11_mie
imax = np.argmax(dolp_mie)
print("=== miepython baseline (independent) ===")
print(f"  x_host          = {x_host:.4f}")
print(f"  Q_sca           = {qsca:.4f}")
print(f"  g               = {g:.4f}")
print(f"  peak DOLP       = {100*dolp_mie[imax]:.3f}%  at theta = {np.degrees(theta[imax]):.1f} deg")

# --- T-matrix in the sphere limit ------------------------------------------
def scattering_matrix(scat, n_ang=361):
    """Return theta[deg], S11, S12, S22 in the conventional scattering-plane
    convention: incidence along +z (thet0->0), scattered at thet=Theta in the
    phi=0 meridian plane (which then coincides with the scattering plane)."""
    th = np.linspace(0.0, 180.0, n_ang)
    S11 = np.zeros(n_ang); S12 = np.zeros(n_ang); S22 = np.zeros(n_ang)
    scat.thet0 = 1e-4; scat.phi0 = 0.0; scat.phi = 0.0
    for i, t in enumerate(th):
        scat.thet = max(min(t, 179.9999), 1e-4)
        Z = scat.get_Z()
        S11[i] = Z[0, 0]; S12[i] = Z[0, 1]; S22[i] = Z[1, 1]
    return th, S11, S12, S22

def metric_B(th_deg, S11, S22):
    th = np.radians(th_deg)
    integ = np.abs(S11 - S22) / S11 * np.sin(th)
    return np.trapezoid(integ, th) / np.pi

sph = Scatterer(radius=a_eqv, wavelength=lam_m, m=complex(m_rel, 0.0),
                axis_ratio=1.000001, ndgs=2, ddelt=1e-3)
sph.orient = orientation.orient_single
th, S11, S12, S22 = scattering_matrix(sph)
dolp_t = -S12 / S11
jmax = np.argmax(dolp_t)
B = metric_B(th, S11, S22)

# cross sections from the T-matrix
Csca = scatter.sca_xsect(sph)
Cext = scatter.ext_xsect(sph)
geom_area = np.pi*a_eqv**2
print("\n=== T-matrix, axis_ratio = 1.000001 (sphere limit) ===")
print(f"  Q_sca           = {Csca/geom_area:.4f}")
print(f"  Q_ext           = {Cext/geom_area:.4f}")
print(f"  peak DOLP       = {100*dolp_t[jmax]:.3f}%  at theta = {th[jmax]:.1f} deg")
print(f"  metric B        = {B:.3e}   (must be ~0 for a sphere)")

print("\n=== agreement ===")
print(f"  peak DOLP diff  = {abs(100*dolp_t[jmax]-100*dolp_mie[imax]):.3f} pts")
print(f"  Q_sca diff      = {abs(Csca/geom_area - qsca):.5f}")
