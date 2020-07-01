#!/usr/bin/python
########################################################################################################################
# Makes a plot of different potential energy functions.
########################################################################################################################
import math
import numpy as np
import matplotlib.pyplot as plt

k = 2743.0
d_0 = 258.9
a = 2.301
r_eq = 1.1283

def harmonic_potential(r):
    return 0.5 * k * (r - r_eq) ** 2

def morse_potential(r):
    return d_0 * (1 - math.exp(-a*(r - r_eq)))**2

x = np.arange(0.95, 1.66, 0.001)
y_harmonic = harmonic_potential(x)
y_morse = [morse_potential(i) for i in x.astype(np.float)]

fig, ax = plt.subplots()
ax.plot(x, y_harmonic, x, y_morse)
ax.legend(['Harmonic', 'Morse'])

ax.set(xlabel="r (A)", ylabel="Energy")
fig.savefig("Potential_Plot.png")
plt.show()