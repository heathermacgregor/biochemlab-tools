#!/usr/bin/python
########################################################################################################################
# 2D steepest descent search using the harmonic potential.
########################################################################################################################
import sys

if len(sys.argv) != 3:
   sys.exit("Usage: Harmonic_SD2.py param_file initial_r")

par_file = open(sys.argv[1])
r_ini   = float(sys.argv[2])

parameters = par_file.readlines()

for parameter in parameters:
    fields = parameter.split()
    mol  = (fields[0])                        # Name of the molecule or bond
    k    = float(fields[1])                  # Harmonic force constant
    r_eq = float(fields[2])                  # Potential minimum distance

def Hamiltonian(r):
    pot = 0.5 * k * (r - r_eq)**2
    return pot

def Force(r):
    grad = k * (r - r_eq)                            # Force   (first derivative)
    return grad

def SDS(r):
    step = 0.0001
    E = Hamiltonian(r)
    F = Force(r)
    print("STEP  DISTANCE   ENERGY   FORCE ")
    print("  0   %6.4f      %8.3f    %8.3f " % (r, E, F))
    for i in range(1, 50):                                       # Iterate up to 50 times
        r -= step*F                                             # Steepest Descent update
        E = Hamiltonian(r)
        F = Force(r)
        print("%3d   %6.4f      %8.3f    %8.3f " % (i, r, E, F))
        if (abs(F) < 0.01):                                     # Are we at the stationary point already?
            break
    return

SDS(r_ini)