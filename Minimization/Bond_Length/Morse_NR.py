#!/usr/bin/python
########################################################################################################################
# Simple Newton-Rhapson minimization using the Morse potential.
########################################################################################################################
import sys
import math
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
    sys.exit("Usage: Morse_NR.py param_file initial_r")

par_file = open(sys.argv[1])
r_ini = float(sys.argv[2])

parameters = par_file.readlines()

for parameter in parameters:
    fields = parameter.split()
    mol = (fields[0])  # Name of the molecule or bond
    d_0 = float(fields[1])  # Dissociation energy
    a = float(fields[2])    # Morse width parameter
    r_eq = float(fields[3])  # Potential minimum distance


def morse_potential(r):
    return d_0 * (1 - math.exp(-a*(r - r_eq)))**2


def force(r): # First derivative
    return 2 * a * d_0 * math.exp(-a*(r - r_eq))*(1 - math.exp(-a*(r - r_eq)))


def hessian(r): # Second derivative
    h = d_0*(2 * a**2 * math.exp(-2 * a * (r - r_eq)) - 2 * a**2 * math.exp(-a * (r - r_eq)) * (1 - math.exp(-a * (r - r_eq))))
    ih = 1/h
    return h, ih


def newton_rhapson(r):
    x = [0]
    y = [r]
    E = morse_potential(r)
    F = force(r)
    H, iH = hessian(r)
    print("STEP  DISTANCE   ENERGY       FORCE ")
    print("%3d  %8.4f    %5.3f   %8.3f " % (0, r, E, F))
    for i in range(1, 10000000):
        r -= F * iH
        E = morse_potential(r)
        F = force(r)
        H, iH = hessian(r)
        print("%3d  %8.4f     %5.3f   %8.3f " % (i, r, E, F))
        x.append(i)
        y.append(r)
        if abs(F) < 0.01 or r < 0:
            break
    return x, y

x, y = newton_rhapson(r_ini)



fig, ax = plt.subplots()
ax.plot(x, y)
ax.set(xlabel='step', ylabel='r (A)', title=str('NR Minimization with Morse Potential from '+str(r_ini)+' A.'))

fig.savefig(str("Morse_NR_" + str(r_ini) + ".png"))
plt.show()
