#/usr/bin/python
########################################################################################################################
# 1D steepest descent search using the harmonic potential.
########################################################################################################################
import sys

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

# Steepest Descent search
r = r_ini
step = 0.0001
E = 0.5 * k * (r - r_eq)**2                   # Energy
F = k * (r - r_eq)                            # Force   (first derivative)
print("Step Distance Energy Force Hessian ")
print("0", r, E, F)                            # Initial values
for i in range(1,50):                         # Iterate up to 50 times
   r = r - step*F                             # Steepest Descent update
   E = 0.5 * k * (r - r_eq)**2
   F = k * (r - r_eq)
   print(i, r, E, F)
   if abs(F) < 0.01:                        # Are we at the stationary point already?
      break