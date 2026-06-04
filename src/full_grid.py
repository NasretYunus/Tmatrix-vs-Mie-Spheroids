"""Resumable full (xi, kappa, lambda) grid for the zeolite-PMMA system.

Grid:
  xi    : 12 log-spaced levels, 0.2 -> 5.0
  kappa : 15 levels, 0 -> 50
  lambda: 400, 550, 800 nm  (band edges + centre)

Engine: Mishchenko single-particle null-field T-matrix (pytmatrix), validated
to machine precision against miepython in the sphere limit.

Runs one chunk (one wavelength, a slice of xi-rows) per invocation so each call
finishes inside the time limit; results persist to full_grid.json after every
cell and the script skips cells already present, so it is fully resumable.

Usage:  python full_grid.py <lam_nm> <xi_start_idx> <xi_count>
"""
import numpy as np, json, time, sys, os
from pytmatrix.tmatrix import Scatterer
from pytmatrix import orientation
from scipy.integrate import quad

n_h, n_p = 1.50, 1.47
m_rel = complex(n_p/n_h, 0.0)
a_eqv = 0.5
TH = np.linspace(0.0, 180.0, 181)
THR = np.radians(TH)
OUT = 'full_grid.json'

XI    = [round(x,4) for x in np.geomspace(0.2, 5.0, 12)]
KAPPA = [round(k,3) for k in np.linspace(0.0, 50.0, 15)]
LAMS  = [400, 550, 800]

def vmf_pdf(kappa):
    if kappa < 1e-9:
        return orientation.uniform_pdf()
    def base(bd):
        b=np.radians(bd); return np.exp(kappa*np.cos(b))*np.sin(b)
    Z=quad(base,0.0,180.0)[0]
    return lambda bd: base(bd)/Z

def angular(s):
    S11=np.zeros(181); S12=np.zeros(181); S22=np.zeros(181)
    s.thet0=1e-4; s.phi0=0.0; s.phi=0.0
    for i,t in enumerate(TH):
        s.thet=max(min(t,179.9999),1e-4); Z=s.get_Z()
        S11[i]=Z[0,0]; S12[i]=Z[0,1]; S22[i]=Z[1,1]
    return S11,S12,S22

def metric_B(S11,S22):
    return float(np.trapezoid(np.abs(S11-S22)/S11*np.sin(THR),THR)/np.pi)

def peak(S11,S12):
    d=-S12/S11; j=int(np.argmax(d)); return float(100*d[j]), float(TH[j])

def cell_tmat(xi, kp, lam_nm):
    lam_m = (lam_nm/1000.0)/n_h
    s=Scatterer(radius=a_eqv, wavelength=lam_m, m=m_rel,
                axis_ratio=1.0/xi, ndgs=2, ddelt=1e-3)
    s.or_pdf = vmf_pdf(kp)
    s.orient = orientation.orient_averaged_fixed
    s.n_alpha=20; s.n_beta=20
    S11,S12,S22 = angular(s)
    P,thm = peak(S11,S12); B = metric_B(S11,S22)
    return P, thm, B

def load():
    if os.path.exists(OUT):
        return json.load(open(OUT))
    return {"_meta":{"XI":XI,"KAPPA":KAPPA,"LAMS":LAMS,"m_rel":0.98,
                     "a_eqv_um":a_eqv,"n_alpha":20,"n_beta":20,"n_angles":181},
            "cells":{}}

if __name__=='__main__':
    lam = int(sys.argv[1]); xi0 = int(sys.argv[2]); xin = int(sys.argv[3])
    d = load()
    done = 0; ran = 0; t_start = time.time()
    for xi in XI[xi0:xi0+xin]:
        for kp in KAPPA:
            key = f"l{lam}_xi{xi}_k{kp}"
            if key in d["cells"]:
                done += 1; continue
            t0=time.time()
            P,thm,B = cell_tmat(xi, kp, lam)
            d["cells"][key] = {"lam":lam,"xi":xi,"kappa":kp,
                               "P":round(P,3),"theta_max":thm,"B":B}
            json.dump(d, open(OUT,'w'))
            ran += 1
        print(f"  xi={xi} done ({len(KAPPA)} kappa) "
              f"[{time.time()-t_start:.0f}s elapsed]")
    print(f"chunk lam={lam} xi[{xi0}:{xi0+xin}]: ran {ran}, skipped {done}, "
          f"total cells now {len(d['cells'])}/{len(XI)*len(KAPPA)*len(LAMS)}")
