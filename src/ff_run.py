"""One FDTD far-field run, cached to disk. Args: shape axratio pol withpart res
   pol in {x,y}; withpart in {0,1}. Saves far-field array to a .npy keyed by args.
"""
import meep as mp, numpy as np, sys, os

n_h,n_p=1.50,1.47; lam_h=1.0; fcen=1.0; df=0.15
shape=sys.argv[1]; axratio=float(sys.argv[2]); pol=sys.argv[3]
withpart=int(sys.argv[4]); resolution=int(sys.argv[5])
xparam=float(sys.argv[6]) if len(sys.argv)>6 else np.pi   # size parameter
beta=float(sys.argv[7]) if len(sys.argv)>7 else 0.0       # 0=axis along k, 90=broadside
a_eq=xparam/(2*np.pi)
key=f"ff_{shape}_{axratio}_{pol}_{withpart}_{resolution}_x{xparam:.2f}_b{beta:.0f}.npy"
if os.path.exists(key):
    print(f"CACHED {key}"); sys.exit(0)

dpml=lam_h; pad=0.8*lam_h
if shape=='sphere':
    rmax=a_eq; geom=[mp.Sphere(radius=a_eq,material=mp.Medium(index=n_p))] if withpart else []
else:
    a_perp=a_eq*axratio**(-1/3); c=a_eq*axratio**(2/3); rmax=max(a_perp,c)
    if beta==0.0:      # long (c) axis along z = incidence direction
        size=mp.Vector3(2*a_perp,2*a_perp,2*c)
    else:              # beta=90: long axis along y (broadside to incidence)
        size=mp.Vector3(2*a_perp,2*c,2*a_perp)
    geom=[mp.Ellipsoid(size=size,material=mp.Medium(index=n_p))] if withpart else []
s=2*rmax+2*pad+2*dpml
cell=mp.Vector3(s,s,s); pml=[mp.PML(dpml)]; host=mp.Medium(index=n_h)
r_n2f=rmax+0.4*pad; zsrc=-s/2+dpml+0.2*pad
comp=mp.Ex if pol=='x' else mp.Ey
src=[mp.Source(mp.GaussianSource(fcen,fwidth=df),component=comp,
     center=mp.Vector3(0,0,zsrc),size=mp.Vector3(s,s,0))]
sim=mp.Simulation(cell_size=cell,boundary_layers=pml,sources=src,resolution=resolution,
    default_material=host,geometry=geom)
n2f=sim.add_near2far(fcen,0,1,
  mp.Near2FarRegion(mp.Vector3(0,0,+r_n2f),size=mp.Vector3(2*r_n2f,2*r_n2f,0),weight=+1),
  mp.Near2FarRegion(mp.Vector3(0,0,-r_n2f),size=mp.Vector3(2*r_n2f,2*r_n2f,0),weight=-1),
  mp.Near2FarRegion(mp.Vector3(+r_n2f,0,0),size=mp.Vector3(0,2*r_n2f,2*r_n2f),weight=+1),
  mp.Near2FarRegion(mp.Vector3(-r_n2f,0,0),size=mp.Vector3(0,2*r_n2f,2*r_n2f),weight=-1),
  mp.Near2FarRegion(mp.Vector3(0,+r_n2f,0),size=mp.Vector3(2*r_n2f,0,2*r_n2f),weight=+1),
  mp.Near2FarRegion(mp.Vector3(0,-r_n2f,0),size=mp.Vector3(2*r_n2f,0,2*r_n2f),weight=-1))
sim.run(until_after_sources=mp.stop_when_fields_decayed(20,comp,mp.Vector3(0,0,r_n2f),1e-4))
Rff=1000.0; thetas=np.linspace(30,150,13)
F=[]
for t in thetas:
    tr=np.radians(t); p=mp.Vector3(Rff*np.sin(tr),0,Rff*np.cos(tr))
    ff=sim.get_farfield(n2f,p); F.append([ff[0],ff[1],ff[2]])
np.save(key,np.array(F))
print(f"SAVED {key}")
