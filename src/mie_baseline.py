"""
Volume-equivalent Mie baseline for spheroidal zeolite-doped PMMA film (xi = 1 row).
Zeolite n_p = 1.47, PMMA host n_h = 1.50  ->  relative index m = 0.98 (particle rarer than host).
Equivalent-sphere diameter d = 1.0 um. Visible band 400-800 nm.
Size parameter uses the HOST medium: x = pi * d * n_h / lambda0.
"""
import numpy as np
import miepython as mp

n_h = 1.50            # PMMA host
n_p = 1.47            # zeolite inclusion
m   = n_p / n_h       # relative refractive index seen by Mie theory
d   = 1.0e-6          # equivalent-sphere diameter (m)

def size_param(lam0_nm):
    lam0 = lam0_nm * 1e-9
    return np.pi * d * n_h / lam0   # size parameter in the host medium

theta = np.linspace(0, np.pi, 1801)          # 0.1 deg resolution
mu = np.cos(theta)

def analyse(lam0_nm):
    x = size_param(lam0_nm)
    qext, qsca, qback, g = mp.efficiencies_mx(m, x)
    P = mp.phase_matrix(m, x, mu, norm='bohren')   # 4x4xN Mueller matrix
    S11 = P[0, 0, :]
    S12 = P[0, 1, :]
    S22 = P[1, 1, :]
    S33 = P[2, 2, :]
    S34 = P[2, 3, :]
    S44 = P[3, 3, :]
    dolp = -S12 / S11                              # DOLP for unpolarised incidence
    pk = np.nanmax(dolp)
    pk_ang = np.degrees(theta[np.nanargmax(dolp)])
    # sphere symmetry checks (should be ~0): these define B = 0 baseline
    sym_22 = np.nanmax(np.abs(S22 - S11) / S11)
    sym_44 = np.nanmax(np.abs(S44 - S33) / np.nanmax(np.abs(S33)))
    Csca = qsca * np.pi * (d/2)**2                 # cross-section, m^2
    Cext = qext * np.pi * (d/2)**2
    return dict(lam=lam0_nm, x=x, qext=qext, qsca=qsca, g=g,
                Cext=Cext, Csca=Csca, peak_dolp=pk, peak_ang=pk_ang,
                dolp90=float(-S12[np.argmin(np.abs(theta-np.pi/2))]/S11[np.argmin(np.abs(theta-np.pi/2))]),
                sym22=sym_22, sym44=sym_44)

print(f"Relative index m = n_p/n_h = {m:.4f}  (particle optically RARER than host)")
print("="*92)
hdr = f"{'lam0(nm)':>8} {'x':>7} {'Qext':>8} {'Qsca':>8} {'g':>7} {'peakDOLP':>9} {'@deg':>6} {'DOLP90':>8} {'symB':>9}"
print(hdr); print("-"*92)
rows = []
for lam in [400, 450, 500, 550, 600, 650, 700, 750, 800]:
    r = analyse(lam)
    rows.append(r)
    print(f"{r['lam']:8d} {r['x']:7.3f} {r['qext']:8.4f} {r['qsca']:8.4f} {r['g']:7.4f} "
          f"{r['peak_dolp']:9.4f} {r['peak_ang']:6.1f} {r['dolp90']:8.4f} {r['sym22']:9.2e}")
print("="*92)
c = [r for r in rows if r['lam']==550][0]
print(f"\nBAND-CENTRE 550 nm (primary baseline, xi = 1):")
print(f"  size parameter x            = {c['x']:.4f}")
print(f"  Q_ext                       = {c['qext']:.4f}")
print(f"  Q_sca                       = {c['qsca']:.4f}   (near-equal to Q_ext: negligible absorption)")
print(f"  asymmetry g                 = {c['g']:.4f}   (strongly forward-peaked)")
print(f"  scattering cross-section    = {c['Csca']*1e12:.4f} um^2")
print(f"  peak DOLP                   = {c['peak_dolp']*100:.2f} %  at theta = {c['peak_ang']:.1f} deg")
print(f"  symmetry-breaking metric B  = {c['sym22']:.2e}  (== 0 to round-off: sphere symmetry confirmed)")
print(f"  S44=S33 check               = {c['sym44']:.2e}")

# spectral spread of peak DOLP for section 4.6
pk = np.array([r['peak_dolp'] for r in rows])
print(f"\nSpectral spread of peak DOLP across 400-800 nm: "
      f"{pk.min()*100:.2f}% to {pk.max()*100:.2f}%  (tilt {100*(pk.max()-pk.min()):.2f} pts)")

np.save('mie_rows.npy', rows, allow_pickle=True)
