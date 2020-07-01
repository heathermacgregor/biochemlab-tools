#!/usr/bin/python
########################################################################################################################
# Simple Newton-Rhapson minimization using the harmonic potential.
# Not working at the moment.
########################################################################################################################
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    sys.exit("Usage: Harmonic_NR.py param_file initial_r")

par_file = open(sys.argv[1])
r_ini = float(sys.argv[2])

parameters = par_file.readlines()

for parameter in parameters:
    fields = parameter.split()
    mol = (fields[0])  # Name of the molecule or bond
    k = float(fields[1])  # Harmonic force constant
    r_eq = float(fields[2])  # Potential minimum distance


def harmonic_potential(r):
    return 0.5 * k * (r - r_eq) ** 2


def force(r): # First derivative
    return k * (r - r_eq)


def hessian(r): # Second derivative
    return k, 1/k
    #hes = np.array([[k,0,0],[0,0,0],[0,0,0]])
    #inv = np.linalg.inv(hes)
    #return hes, inv

def hessian_bfgs(r, H_old, iH_old):
    p_new = -1*force(r)*iH_old
    alpha = .01
    s_new = alpha * p_new
    r_new = r + s_new
    y_new = force(r_new) - force(r)
    H_new = y_new * y_new / (y_new * s_new)
    iH_new = 1/H_new
    return H_new, iH_new



x = [0]
y = [r_ini]
def newton_rhapson(r):
    E = harmonic_potential(r)
    F = force(r)
    H, iH = hessian(r)
    print("STEP  DISTANCE   ENERGY   FORCE ")
    print("  0   %6.4f      %8.3f    %8.3f " % (r, E, F))
    for i in range(1, 100):
        r -= F*iH
        E = harmonic_potential(r)
        F = force(r)
        H, iH = hessian(r)
        x.append(i)
        y.append(r)
        print("%3d   %6.4f      %8.3f    %8.3f " % (i, r, E, F))
        if (abs(F) < 0.01):
            break
    return

newton_rhapson(r_ini)
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set(xlabel='step', ylabel='r (A)', title=str('NR Minimization with Harmonic Potential from '+str(r_ini)+' A.'))

fig.savefig(str("Harmonic_NR_" + str(r_ini) + ".png"))
plt.show()

def bfgs(r):
    E = harmonic_potential(r)
    F = force(r)
    H, iH = hessian_bfgs(r, 1, 1)
    print("STEP  DISTANCE   ENERGY   FORCE ")
    print("  0   %6.4f      %8.3f    %8.3f " % (r, E, F))
    for i in range(1, 100):
        r -= F * iH
        E = harmonic_potential(r)
        F = force(r)
        x.append(i)
        y.append(r)
        print("%3d   %6.4f      %8.3f    %8.3f " % (i, r, E, F))
        if (abs(F) < 0.01):
            break
        H, iH = hessian_bfgs(r, H, iH)

    return

bfgs(r_ini)