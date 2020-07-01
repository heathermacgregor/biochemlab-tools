#!/usr/bin/env python3
#
# Newton-Rhapson minimizer of the Rosenbrock function

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Make data.
s = 0.05  # Try s=1, 0.25, 0.1, or 0.05
X_n = np.arange(-2, 2. + s, s)  # Could use linspace instead if dividing
Y_n = np.arange(-2, 3. + s, s)  # evenly instead of stepping...

# Create the mesh grid(s) for all X/Y combos.
X_nY_n = np.meshgrid(X_n, Y_n)
# print(len(X_nY_n[0]))
Z_n = []
G_n = []
for i in range(len(X_nY_n[0])):
    x = X_nY_n[0]
    y = X_nY_n[1]
    Z_n.append(rosenbrock(x[i], y[i]))
    G_n.append(gradient(x[i], y[i]))
Z_n = np.transpose(Z_n)
G_n = np.transpose(G_n)
# X_nY_n = np.transpose(X_nY_n)

fig, ax = plt.subplots(figsize=(10, 6))
print(Z_n[0])
print(G_n[0])
print(len(G_n))
print(len(X_nY_n))

ax.contour(X_nY_n, Z_n, levels=np.logspace(0, 5, 35), norm=LogNorm(), cmap=plt.cm.jet)
ax.quiver(X_nY_n, X_nY_n - G_n, alpha=0.5)
ax.plot(np.array([1, 1]), 'r*', markersize=18)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

ax.set_xlim((np.min(x), np.max(x)))
ax.set_ylim((np.min(y), np.max(y)))

plt.show()
