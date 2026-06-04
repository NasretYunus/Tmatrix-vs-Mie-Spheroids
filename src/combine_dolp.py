import numpy as np, sys
shape=sys.argv[1]; ax=sys.argv[2]; res=sys.argv[3]
thetas=np.linspace(30,150,13)
Pe=np.load(f"ff_{shape}_{ax}_y_0_{res}.npy")
Pt=np.load(f"ff_{shape}_{ax}_y_1_{res}.npy")
Qe=np.load(f"ff_{shape}_{ax}_x_0_{res}.npy")
Qt=np.load(f"ff_{shape}_{ax}_x_1_{res}.npy")
Sperp=Pt-Pe; Spar=Qt-Qe
dolp=[]
for i,t in enumerate(np.radians(thetas)):
    I_perp=abs(Sperp[i][1])**2          # E_y scattered, perp incidence
    ex,ey,ez=Spar[i]
    Epar=ex*np.cos(t)-ez*np.sin(t)      # in-plane theta component, par incidence
    I_par=abs(Epar)**2
    dolp.append((I_perp-I_par)/(I_perp+I_par) if (I_perp+I_par)>0 else 0)
dolp=np.array(dolp)
j=np.argmax(np.abs(dolp))
print(f"FDTD {shape} ax={ax} res={res}: peak DOLP={100*dolp[j]:.1f}% at {thetas[j]:.0f}deg")
print("curve:", "  ".join(f"{t:.0f}:{100*d:.0f}%" for t,d in zip(thetas,dolp)))
np.save(f"dolp_{shape}_{ax}_{res}.npy",np.column_stack([thetas,dolp]))

if shape=='sphere':
    import miepython
    x=np.pi; m=complex(0.98,0)
    mu=np.cos(np.radians(thetas)); S1,S2=miepython.S1_S2(m,x,mu)
    S11=0.5*(abs(S1)**2+abs(S2)**2); S12=0.5*(abs(S2)**2-abs(S1)**2)
    md=-S12/S11
    print("Mie :  ", "  ".join(f"{t:.0f}:{100*d:.0f}%" for t,d in zip(thetas,md)))
    print(f"Mie peak in window: {100*np.max(np.abs(md)):.1f}%")
