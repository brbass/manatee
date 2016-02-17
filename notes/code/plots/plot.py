import numpy as np
from matplotlib import pyplot as plt

class mc:
    cell = 0
    group = 1
    
    cur = 2
    sur = 3
    vol = 4
    col = 5
    
    cur_var = 6
    sur_var = 7
    vol_var = 8
    col_var = 9
    
    cur_err = 10
    sur_err = 11
    vol_err = 12
    col_err = 13

    cur_rel = 14
    sur_rel = 15
    vol_rel = 16
    col_rel = 17

    cur_fom = 18
    sur_fom = 19
    vol_fom = 20
    col_fom = 21

    g0 = range(300)
    g1 = range(300, 600)
    
class dt:
    cell = 0
    group = 1

    sp1ad = 2
    sp1fw = 3

    sn = 4

    g0 = range(300)
    g1 = range(300, 600)
    
case = "thick"
dett = np.loadtxt(case + "/deterministic.txt", skiprows=2)
im0ad0t = np.loadtxt(case + "/im0ad0.txt", skiprows=2)
im1ad0t = np.loadtxt(case + "/im1ad0.txt", skiprows=2)
im0ad1t = np.loadtxt(case + "/im0ad1.txt", skiprows=2)
im1ad1t = np.loadtxt(case + "/im1ad1.txt", skiprows=2)

case = "scattering"
dets = np.loadtxt(case + "/deterministic.txt", skiprows=2)
im0ad0s = np.loadtxt(case + "/im0ad0.txt", skiprows=2)
im1ad0s = np.loadtxt(case + "/im1ad0.txt", skiprows=2)
im0ad1s = np.loadtxt(case + "/im0ad1.txt", skiprows=2)
im1ad1s = np.loadtxt(case + "/im1ad1.txt", skiprows=2)

xt = np.linspace(0, 10, 300)
xs = np.linspace(0, 1, 300)

w, h = plt.figaspect(1/2)

# estimators: implicit and adjoint
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im1ad1t[mc.g0,mc.sur], label="Surface")
plt.plot(xt, im1ad1t[mc.g0,mc.vol], label="Volume")
plt.plot(xt, im1ad1t[mc.g0,mc.col], label="Collision")
plt.plot(xt, dett[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend()

plt.subplot(122)
plt.plot(xt, im1ad1t[mc.g1,mc.sur], label="Surface")
plt.plot(xt, im1ad1t[mc.g1,mc.vol], label="Volume")
plt.plot(xt, im1ad1t[mc.g1,mc.col], label="Collision")
plt.plot(xt, dett[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend()
plt.tight_layout()
plt.savefig("thick_estimators_ad.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im1ad1s[mc.g0,mc.sur], label="Surface")
plt.plot(xs, im1ad1s[mc.g0,mc.vol], label="Volume")
plt.plot(xs, im1ad1s[mc.g0,mc.col], label="Collision")
plt.plot(xs, dets[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g0,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im1ad1s[mc.g1,mc.sur], label="Surface")
plt.plot(xs, im1ad1s[mc.g1,mc.vol], label="Volume")
plt.plot(xs, im1ad1s[mc.g1,mc.col], label="Collision")
plt.plot(xs, dets[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g1,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_estimators_ad.pdf")

# estimators: implicit with no adjoint
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im1ad0t[mc.g0,mc.sur], label="Surface")
plt.plot(xt, im1ad0t[mc.g0,mc.vol], label="Volume")
plt.plot(xt, im1ad0t[mc.g0,mc.col], label="Collision")
plt.plot(xt, dett[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend()

plt.subplot(122)
plt.plot(xt, im1ad0t[mc.g1,mc.sur], label="Surface")
plt.plot(xt, im1ad0t[mc.g1,mc.vol], label="Volume")
plt.plot(xt, im1ad0t[mc.g1,mc.col], label="Collision")
plt.plot(xt, dett[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend()
plt.tight_layout()
plt.savefig("thick_estimators_norm.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im1ad0s[mc.g0,mc.sur], label="Surface")
plt.plot(xs, im1ad0s[mc.g0,mc.vol], label="Volume")
plt.plot(xs, im1ad0s[mc.g0,mc.col], label="Collision")
plt.plot(xs, dets[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g0,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im1ad0s[mc.g1,mc.sur], label="Surface")
plt.plot(xs, im1ad0s[mc.g1,mc.vol], label="Volume")
plt.plot(xs, im1ad0s[mc.g1,mc.col], label="Collision")
plt.plot(xs, dets[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g1,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_estimators_norm.pdf")

# implicit and adjoint options: scalar flux
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im0ad0t[mc.g0,mc.vol], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g0,mc.vol], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g0,mc.vol], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g0,mc.vol], label="Imp=1, Ad=1")
plt.plot(xt, dett[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend()

plt.subplot(122)
plt.plot(xt, im0ad0t[mc.g1,mc.vol], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g1,mc.vol], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g1,mc.vol], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g1,mc.vol], label="Imp=1, Ad=1")
plt.plot(xt, dett[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend()
plt.tight_layout()
plt.savefig("thick_imad.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im0ad0s[mc.g0,mc.vol], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g0,mc.vol], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g0,mc.vol], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g0,mc.vol], label="Imp=1, Ad=1")
plt.plot(xs, dets[dt.g0,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g0,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im0ad0s[mc.g1,mc.vol], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g1,mc.vol], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g1,mc.vol], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g1,mc.vol], label="Imp=1, Ad=1")
plt.plot(xs, dets[dt.g1,dt.sp1fw], label=r"SP$_1$")
plt.plot(xs, dets[dt.g1,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_imad.pdf")

# implicit and adjoint options: figure of merit
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im0ad0t[mc.g0,mc.vol_fom], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g0,mc.vol_fom], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g0,mc.vol_fom], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g0,mc.vol_fom], label="Imp=1, Ad=1")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$FOM(\phi(x))$')
plt.title('Group 1')
plt.legend(loc='lower left')

plt.subplot(122)
plt.plot(xt, im0ad0t[mc.g1,mc.vol_fom], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g1,mc.vol_fom], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g1,mc.vol_fom], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g1,mc.vol_fom], label="Imp=1, Ad=1")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$FOM(\phi(x))$')
plt.title('Group 2')
plt.legend(loc='lower left')
plt.tight_layout()
plt.savefig("thick_fom.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im0ad0s[mc.g0,mc.vol_fom], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g0,mc.vol_fom], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g0,mc.vol_fom], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g0,mc.vol_fom], label="Imp=1, Ad=1")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$FOM(\phi(x))$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im0ad0s[mc.g1,mc.vol_fom], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g1,mc.vol_fom], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g1,mc.vol_fom], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g1,mc.vol_fom], label="Imp=1, Ad=1")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$FOM(\phi(x))$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_fom.pdf")

# error in volume flux
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im0ad0t[mc.g0,mc.vol_err], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g0,mc.vol_err], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g0,mc.vol_err], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g0,mc.vol_err], label="Imp=1, Ad=1")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$Err(\phi(x))$')
plt.title('Group 1')
plt.legend()

plt.subplot(122)
plt.plot(xt, im0ad0t[mc.g1,mc.vol_err], label="Imp=0, Ad=0")
plt.plot(xt, im1ad0t[mc.g1,mc.vol_err], label="Imp=1, Ad=0")
plt.plot(xt, im0ad1t[mc.g1,mc.vol_err], label="Imp=0, Ad=1")
plt.plot(xt, im1ad1t[mc.g1,mc.vol_err], label="Imp=1, Ad=1")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$Err(\phi(x))$')
plt.title('Group 2')
plt.legend()
plt.tight_layout()
plt.savefig("thick_err.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im0ad0s[mc.g0,mc.vol_err], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g0,mc.vol_err], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g0,mc.vol_err], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g0,mc.vol_err], label="Imp=1, Ad=1")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$Err(\phi(x))$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im0ad0s[mc.g1,mc.vol_err], label="Imp=0, Ad=0")
plt.plot(xs, im1ad0s[mc.g1,mc.vol_err], label="Imp=1, Ad=0")
plt.plot(xs, im0ad1s[mc.g1,mc.vol_err], label="Imp=0, Ad=1")
plt.plot(xs, im1ad1s[mc.g1,mc.vol_err], label="Imp=1, Ad=1")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$Err(\phi(x))$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_err.pdf")

# adjoint functions
plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xt, im1ad1t[mc.g0,mc.vol], label="MC")
plt.plot(xt, dett[dt.g0,dt.sp1fw], label=r"SP$_1$ Forward")
plt.plot(xt, dett[dt.g0,dt.sp1ad], label=r"SP$_1$ Adjoint")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend(loc='right')

plt.subplot(122)
plt.plot(xt, im1ad1t[mc.g1,mc.vol], label="MC")
plt.plot(xt, dett[dt.g1,dt.sp1fw], label=r"SP$_1$ Forward")
plt.plot(xt, dett[dt.g1,dt.sp1ad], label=r"SP$_1$ Adjoint")
plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend(loc='right')
plt.tight_layout()
plt.savefig("thick_adj.pdf")

plt.figure(figsize=(w,h))
plt.subplot(121)
plt.plot(xs, im1ad1s[mc.g0,mc.vol], label="MC")
plt.plot(xs, dets[dt.g0,dt.sp1fw], label=r"SP$_1$ Forward")
plt.plot(xs, dets[dt.g0,dt.sp1ad], label=r"SP$_1$ Adjoint")
plt.plot(xs, dets[dt.g0,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 1')
plt.legend(loc='lower center')

plt.subplot(122)
plt.plot(xs, im1ad1s[mc.g1,mc.vol], label="MC")
plt.plot(xs, dets[dt.g1,dt.sp1fw], label=r"SP$_1$ Forward")
plt.plot(xs, dets[dt.g1,dt.sp1ad], label=r"SP$_1$ Adjoint")
plt.plot(xs, dets[dt.g1,dt.sn], label=r"S$_N$")
#plt.yscale('log')
plt.xlabel('x')
plt.ylabel(r'$\phi(x)$')
plt.title('Group 2')
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("scattering_adj.pdf")

#plt.show()
