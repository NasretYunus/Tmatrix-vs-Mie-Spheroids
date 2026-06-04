"""Figure 3 (replacement): the full computed (xi,kappa) maps.
Panel A: divergence D(xi,kappa) filled contour at 550 nm.
Panel B: symmetry-breaking metric B(xi,kappa) filled contour (log) at 550 nm.
Panel C: spectral comparison -- peak DOLP vs xi for the aligned population
         (kappa=50) at 400/550/800 nm, showing how the cliff deepens to the blue.
All from full_grid.json, 12x15 grid at each wavelength.
"""
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

d = json.load(open('full_grid.json'))
c = d['cells']
XI = d['_meta']['XI']
KAPPA = d['_meta']['KAPPA']
P_MIE = 99.505  # validated kappa-independent equivalent-sphere reference at 550

def grid_at(lam, field):
    M = np.full((len(KAPPA), len(XI)), np.nan)
    for ik, kp in enumerate(KAPPA):
        for ix, xi in enumerate(XI):
            cell = c.get(f"l{lam}_xi{xi}_k{kp}")
            if cell and cell.get(field) is not None:
                M[ik, ix] = cell[field]
    return M

logxi = np.log10(XI)
KK = np.array(KAPPA)
XX, YY = np.meshgrid(logxi, KK)

P550 = grid_at(550, 'P')
B550 = grid_at(550, 'B')
D550 = np.abs(P550 - P_MIE) / P550 * 100.0

fig = plt.figure(figsize=(14, 4.4))

# Panel A: divergence D
axA = fig.add_subplot(1, 3, 1)
lev = [0, 1, 2, 5, 10, 20, 40, 60, 70]
cf = axA.contourf(XX, YY, D550, levels=lev, cmap='YlOrRd', extend='max')
c5 = axA.contour(XX, YY, D550, levels=[5], colors='k', linewidths=2)
axA.clabel(c5, fmt='5%%', fontsize=9)
axA.set_xticks(np.log10([0.2,0.5,1,2,5]))
axA.set_xticklabels(['0.2','0.5','1','2','5'])
axA.set_xlabel('aspect ratio $\\xi$'); axA.set_ylabel('orientation concentration $\\kappa$')
axA.set_title('Divergence $D(\\xi,\\kappa)$ at 550 nm (%)')
fig.colorbar(cf, ax=axA, fraction=0.046, pad=0.04)

# Panel B: symmetry metric B (log)
axB = fig.add_subplot(1, 3, 2)
Blog = np.log10(np.clip(B550, 1e-5, None))
cf2 = axB.contourf(XX, YY, Blog, levels=np.linspace(-5, 0, 11), cmap='viridis', extend='both')
axB.set_xticks(np.log10([0.2,0.5,1,2,5]))
axB.set_xticklabels(['0.2','0.5','1','2','5'])
axB.set_xlabel('aspect ratio $\\xi$'); axB.set_ylabel('orientation concentration $\\kappa$')
axB.set_title('Symmetry breaking $\\log_{10} B(\\xi,\\kappa)$ at 550 nm')
cb2 = fig.colorbar(cf2, ax=axB, fraction=0.046, pad=0.04)

# Panel C: spectral comparison of the aligned (kappa=50) cut
axC = fig.add_subplot(1, 3, 3)
ik50 = len(KAPPA) - 1  # kappa = 50
colors = {400:'tab:blue', 550:'tab:green', 800:'tab:red'}
for lam in [400, 550, 800]:
    P = grid_at(lam, 'P')[ik50]
    mask = ~np.isnan(P)
    axC.plot(np.array(logxi)[mask], P[mask], 'o-', color=colors[lam],
             lw=2, ms=5, label=f'{lam} nm')
axC.axhline(P_MIE, ls='--', c='gray', lw=1, label='Mie ref')
axC.set_xticks(np.log10([0.2,0.5,1,2,5]))
axC.set_xticklabels(['0.2','0.5','1','2','5'])
axC.set_xlabel('aspect ratio $\\xi$'); axC.set_ylabel('peak DOLP (%) at $\\kappa=50$')
axC.set_title('Spectral dependence (aligned)')
axC.legend(fontsize=8); axC.grid(alpha=0.3); axC.set_ylim(50, 102)

plt.tight_layout()
plt.savefig('fig3_grid_maps.png', dpi=150, bbox_inches='tight')
print("saved fig3_grid_maps.png")

# numeric summary for the manuscript
print("\n=== 550 nm summary ===")
print(f"max divergence: {np.nanmax(D550):.1f}% at xi={XI[np.unravel_index(np.nanargmax(D550),D550.shape)[1]]}, "
      f"kappa={KAPPA[np.unravel_index(np.nanargmax(D550),D550.shape)[0]]}")
# fraction of plane under 5%
under5 = np.sum(D550 < 5) / np.sum(~np.isnan(D550)) * 100
print(f"fraction of (xi,kappa) plane with D<5%: {under5:.0f}%")
# how cliff deepens with wavelength at xi=5, kappa=50
for lam in [800,550,400]:
    v = grid_at(lam,'P')[ik50,-1]
    print(f"  xi=5 kappa=50 at {lam}nm: peak DOLP = {v if not np.isnan(v) else 'EBCM-limited'}")
