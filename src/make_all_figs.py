"""Regenerate all four manuscript figures in publication formats:
PNG at 1000 dpi, plus vector PDF, SVG, and EPS.

Reads the same data the paper uses: live Lorenz-Mie for Figs 1-2, full_grid.json
for Fig 3, and the cached FDTD far-field .npy files for Fig 4. Run from a directory
that contains full_grid.json and the ff_*.npy caches (e.g. the repo build/ dir), or
adjust DATA below.
"""
import numpy as np, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import miepython as mie

OUTDIR = os.environ.get('FIGDIR', '.')
FMTS = ['png', 'pdf', 'svg', 'eps']
DPI_PNG = 1000

def save_all(fig, stem):
    for f in FMTS:
        kw = {'dpi': DPI_PNG} if f == 'png' else {}
        fig.savefig(os.path.join(OUTDIR, f"{stem}.{f}"), bbox_inches='tight', **kw)
    plt.close(fig)
    print(f"wrote {stem}.{{{','.join(FMTS)}}}")

# ----------------------------------------------------------------- Fig 1 & 2
n_h, n_p, d = 1.50, 1.47, 1.0e-6
m = n_p/n_h
theta = np.linspace(0, np.pi, 1801); mu = np.cos(theta); deg = np.degrees(theta)
x_of = lambda l: np.pi*d*n_h/(l*1e-9)

# Fig 1: Mueller elements + DOLP at 550 nm
x = x_of(550); P = mie.phase_matrix(m, x, mu, norm='bohren')
S11, S12, S33, S34 = P[0,0], P[0,1], P[2,2], P[2,3]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].semilogy(deg, S11, label='S11'); ax[0].semilogy(deg, np.abs(S12), '--', label='|S12|')
ax[0].semilogy(deg, np.abs(S33), ':', label='|S33|'); ax[0].semilogy(deg, np.abs(S34), '-.', label='|S34|')
ax[0].set_xlabel('scattering angle (deg)'); ax[0].set_ylabel('Mueller element (a.u.)')
ax[0].set_title('Mueller matrix, sphere baseline, 550 nm'); ax[0].legend(fontsize=8); ax[0].set_xlim(0,180)
ax[1].plot(deg, -S12/S11, 'k'); ax[1].set_xlabel('scattering angle (deg)'); ax[1].set_ylabel('DOLP = -S12/S11')
ax[1].set_title('Degree of linear polarization (xi=1)'); ax[1].set_xlim(0,180); ax[1].set_ylim(-1,1.05); ax[1].grid(alpha=.3)
fig.tight_layout(); save_all(fig, 'fig1_mueller_550')

# Fig 2: peak DOLP + Qsca vs wavelength
lams = np.arange(400, 801, 5); pk = []; qs = []
for l in lams:
    xx = x_of(l); _, qsca, _, _ = mie.efficiencies_mx(m, xx); qs.append(qsca)
    Pm = mie.phase_matrix(m, xx, mu, norm='bohren'); pk.append(np.nanmax(-Pm[0,1]/Pm[0,0]))
fig, ax1 = plt.subplots(figsize=(7, 4.2))
ax1.plot(lams, np.array(pk)*100, 'b-'); ax1.set_xlabel('vacuum wavelength (nm)')
ax1.set_ylabel('peak DOLP (%)', color='b'); ax1.tick_params(axis='y', labelcolor='b'); ax1.set_ylim(90,100.5)
ax2 = ax1.twinx(); ax2.plot(lams, qs, 'r--'); ax2.set_ylabel('Q_sca', color='r'); ax2.tick_params(axis='y', labelcolor='r')
ax1.set_title('Spectral dependence, sphere baseline (m=0.98)'); fig.tight_layout()
save_all(fig, 'fig2_spectral')

# ----------------------------------------------------------------- Fig 3
d3 = json.load(open('full_grid.json')); c = d3['cells']
XI = d3['_meta']['XI']; KAPPA = d3['_meta']['KAPPA']; P_MIE = 99.505
def grid_at(lam, field):
    M = np.full((len(KAPPA), len(XI)), np.nan)
    for ik, kp in enumerate(KAPPA):
        for ix, xi in enumerate(XI):
            cell = c.get(f"l{lam}_xi{xi}_k{kp}")
            if cell and cell.get(field) is not None: M[ik, ix] = cell[field]
    return M
logxi = np.log10(XI); KK = np.array(KAPPA); XX, YY = np.meshgrid(logxi, KK)
P550 = grid_at(550,'P'); B550 = grid_at(550,'B'); D550 = np.abs(P550-P_MIE)/P550*100
fig = plt.figure(figsize=(14, 4.4))
axA = fig.add_subplot(1,3,1)
cf = axA.contourf(XX, YY, D550, levels=[0,1,2,5,10,20,40,60,70], cmap='YlOrRd', extend='max')
c5 = axA.contour(XX, YY, D550, levels=[5], colors='k', linewidths=2); axA.clabel(c5, fmt='5%%', fontsize=9)
axA.set_xticks(np.log10([0.2,0.5,1,2,5])); axA.set_xticklabels(['0.2','0.5','1','2','5'])
axA.set_xlabel('aspect ratio $\\xi$'); axA.set_ylabel('orientation concentration $\\kappa$')
axA.set_title('Divergence $D(\\xi,\\kappa)$ at 550 nm (%)'); fig.colorbar(cf, ax=axA, fraction=0.046, pad=0.04)
axB = fig.add_subplot(1,3,2)
Blog = np.log10(np.clip(B550, 1e-5, None))
cf2 = axB.contourf(XX, YY, Blog, levels=np.linspace(-5,0,11), cmap='viridis', extend='both')
axB.set_xticks(np.log10([0.2,0.5,1,2,5])); axB.set_xticklabels(['0.2','0.5','1','2','5'])
axB.set_xlabel('aspect ratio $\\xi$'); axB.set_ylabel('orientation concentration $\\kappa$')
axB.set_title('Symmetry breaking $\\log_{10} B(\\xi,\\kappa)$ at 550 nm'); fig.colorbar(cf2, ax=axB, fraction=0.046, pad=0.04)
axC = fig.add_subplot(1,3,3); ik50 = len(KAPPA)-1
for lam, col in [(400,'tab:blue'),(550,'tab:green'),(800,'tab:red')]:
    Pr = grid_at(lam,'P')[ik50]; mask = ~np.isnan(Pr)
    axC.plot(np.array(logxi)[mask], Pr[mask], 'o-', color=col, lw=2, ms=5, label=f'{lam} nm')
axC.axhline(P_MIE, ls='--', c='gray', lw=1, label='Mie ref')
axC.set_xticks(np.log10([0.2,0.5,1,2,5])); axC.set_xticklabels(['0.2','0.5','1','2','5'])
axC.set_xlabel('aspect ratio $\\xi$'); axC.set_ylabel('peak DOLP (%) at $\\kappa=50$')
axC.set_title('Spectral dependence (aligned)'); axC.legend(fontsize=8); axC.grid(alpha=0.3); axC.set_ylim(50,102)
fig.tight_layout(); save_all(fig, 'fig3_grid_maps')

# ----------------------------------------------------------------- Fig 4
try:
    from pytmatrix.tmatrix import Scatterer
    from pytmatrix import orientation
    th = np.linspace(30,150,13)
    def dolp_from(tag):
        Pe=np.load(f"ff_{tag}_y_0.npy") if os.path.exists(f"ff_{tag}_y_0.npy") else None
        return None
    # sphere x=pi
    def load4(prefix):
        return (np.load(f"{prefix}_y_0.npy"), np.load(f"{prefix}_y_1.npy"),
                np.load(f"{prefix}_x_0.npy"), np.load(f"{prefix}_x_1.npy"))
    def dolp(Pe,Pt,Qe,Qt):
        Sp=Pt-Pe; Sq=Qt-Qe; out=[]
        for i,t in enumerate(np.radians(th)):
            Ip=abs(Sp[i][1])**2; ex,ey,ez=Sq[i]; Ipar=abs(ex*np.cos(t)-ez*np.sin(t))**2
            out.append((Ip-Ipar)/(Ip+Ipar) if (Ip+Ipar)>0 else 0)
        return np.array(out)
    ds=dolp(*load4('ff_sphere_1.0')) if os.path.exists('ff_sphere_1.0_y_0.npy') else None
    # sphere may be stored with resolution suffix; try common ones
    if ds is None:
        for r in ['20','14','18','16']:
            if os.path.exists(f'ff_sphere_1.0_y_0_{r}.npy'):
                ds=dolp(np.load(f'ff_sphere_1.0_y_0_{r}.npy'),np.load(f'ff_sphere_1.0_y_1_{r}.npy'),
                        np.load(f'ff_sphere_1.0_x_0_{r}.npy'),np.load(f'ff_sphere_1.0_x_1_{r}.npy')); break
    mu2=np.cos(np.radians(th)); S1,S2=mie.S1_S2(complex(0.98,0),np.pi,mu2)
    S11m=0.5*(abs(S1)**2+abs(S2)**2); S12m=0.5*(abs(S2)**2-abs(S1)**2); dms=-S12m/S11m
    tag="10_x8.57_b90"
    dn=dolp(np.load(f"ff_prolate_5.0_y_0_{tag}.npy"),np.load(f"ff_prolate_5.0_y_1_{tag}.npy"),
            np.load(f"ff_prolate_5.0_x_0_{tag}.npy"),np.load(f"ff_prolate_5.0_x_1_{tag}.npy"))
    a=8.57/(2*np.pi); s=Scatterer(radius=a,wavelength=1.0,m=complex(0.98,0),axis_ratio=0.2,ndgs=3,ddelt=1e-3)
    s.orient=orientation.orient_single; s.thet0=1e-4;s.phi0=0;s.phi=0;s.alpha=0;s.beta=90.0
    dnt=[]
    for t in th: s.thet=t; Z=s.get_Z(); dnt.append(-Z[0,1]/Z[0,0])
    dnt=np.array(dnt)
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    ax[0].plot(th,100*dms,'k-',lw=2,label='Mie (exact)')
    if ds is not None: ax[0].plot(th,100*ds,'o',color='tab:red',ms=7,label='FDTD (Meep)')
    ax[0].set_title('Validation: sphere, x=3.14'); ax[0].set_xlabel('scattering angle (deg)')
    ax[0].set_ylabel('DOLP (%)'); ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    ax[1].plot(th,100*dnt,'k-',lw=2,label='T-matrix')
    ax[1].plot(th,100*dn,'s',color='tab:blue',ms=7,label='FDTD (Meep)')
    ax[1].set_title('Cross-check: prolate needle $\\xi$=5 broadside, x=8.57')
    ax[1].set_xlabel('scattering angle (deg)'); ax[1].set_ylabel('DOLP (%)'); ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
    fig.tight_layout(); save_all(fig,'fig4_fdtd')
except Exception as e:
    print(f"Fig 4 skipped (needs pytmatrix + cached FDTD .npy files): {e}")

print("done.")
