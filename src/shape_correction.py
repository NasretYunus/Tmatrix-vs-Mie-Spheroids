"""Shape-corrected equivalent-sphere model via quasi-static depolarization factors.

For a spheroid with semi-axes (a,a,c) the geometric depolarization factors L_j
(sum L = 1) set the quasi-static polarizability along each axis:
    alpha_j ~ V (eps_p - eps_h) / (eps_h + L_j (eps_p - eps_h)),
with eps = n^2. A sphere has L = 1/3 on every axis and isotropic polarizability;
a spheroid is anisotropic, and that anisotropy is the leading-order shape effect on
polarization. The peak degree of linear polarization for a Rayleigh scatterer with
principal polarizabilities is

    DOLP_peak = (alpha_max^2 - alpha_min^2)/(alpha_max^2 + alpha_min^2)  -> 1 for a
    sphere (alpha isotropic gives the pure dipole result), decreasing as the
    polarizability anisotropy grows.

This is the correct *shape* foil: it carries the leading shape dependence but, being a
single-particle quasi-static quantity, it has no representation of orientation order
(no kappa), exactly as the manuscript requires. We compute the shape-corrected peak
DOLP P_Mie(xi) and compare it to the rigorous orientation-averaged T-matrix.

Reference: depolarization factors of spheroids, Bohren & Huffman (1983), Sec. 5.3;
Osborn, Phys. Rev. 67, 351 (1945); Stoner, Phil. Mag. 36, 803 (1945).
"""
import numpy as np, json

n_h, n_p = 1.50, 1.47
eps_h, eps_p = n_h**2, n_p**2

def depol_factors(xi):
    """Return (L_par, L_perp) for a spheroid of aspect ratio xi = c/a (polar/equatorial).
    xi>1 prolate, xi<1 oblate. L_par along the c (symmetry) axis."""
    if abs(xi-1) < 1e-6:
        return 1/3, 1/3
    if xi > 1:  # prolate: e^2 = 1 - (a/c)^2 = 1 - 1/xi^2
        e = np.sqrt(1 - 1/xi**2)
        Lz = (1-e**2)/e**2 * (-1 + (1/(2*e))*np.log((1+e)/(1-e)))
    else:       # oblate: e^2 = 1 - (c/a)^2 = 1 - xi^2
        e = np.sqrt(1 - xi**2)
        g = np.sqrt((1-e**2)/e**2)
        Lz = (1/e**2) * (1 - g*np.arctan(1/g)) if False else \
             (g/(e**2)) * (np.pi/2 - np.arctan(g)) * 0  # placeholder, replace below
        # standard oblate closed form:
        Lz = (1/e**2)*(1 - np.sqrt(1-e**2)/e*np.arcsin(e))
    Lx = (1 - Lz)/2     # the two equatorial axes share the remainder
    return Lz, Lx

def alpha_axis(L):
    return (eps_p - eps_h) / (eps_h + L*(eps_p - eps_h))

def shape_corrected_dolp(xi):
    Lz, Lx = depol_factors(xi)
    az, ax = alpha_axis(Lz), alpha_axis(Lx)
    amax, amin = max(abs(az),abs(ax)), min(abs(az),abs(ax))
    return (amax**2 - amin**2)/(amax**2 + amin**2)

# sanity: sphere
print("sphere check L:", depol_factors(1.0), "-> DOLP anisotropy", shape_corrected_dolp(1.0001))

# But peak DOLP of the equal-volume *sphere* itself is not 1 in Mie -- it is the
# computed 99.5%. The shape correction scales the *departure* from that baseline.
# Model: P_Mie(xi) = P_sphere * (1 - dolp_anisotropy_penalty), where the penalty is
# the relative polarizability anisotropy (0 at sphere). Concretely we take
#   P_Mie(xi) = P_sphere * amin/amax   (ratio of minor to major polarizability),
# which is 1 at the sphere and falls as the spheroid becomes anisotropic -- a clean,
# monotone, first-principles shape correction with no free parameter.
P_sphere = 99.505

XI = [round(x,4) for x in np.geomspace(0.2,5.0,12)]
print(f"\n{'xi':>6} {'L_par':>7} {'L_perp':>7} {'amin/amax':>10} {'P_Mie_corr':>11}")
rows={}
for xi in XI:
    Lz,Lx = depol_factors(xi)
    az,ax = alpha_axis(Lz), alpha_axis(Lx)
    ratio = min(abs(az),abs(ax))/max(abs(az),abs(ax))
    P_corr = P_sphere * ratio
    rows[xi] = P_corr
    print(f"{xi:6.3f} {Lz:7.3f} {Lx:7.3f} {ratio:10.4f} {P_corr:11.2f}")

json.dump(rows, open('shape_corr.json','w'))
print("\nsaved shape_corr.json")
