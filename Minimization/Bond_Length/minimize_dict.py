import numpy as np
import sys

par_file = open(sys.argv[1])
r_ini = float(sys.argv[2])

parameters = par_file.readlines()

for parameter in parameters:
    fields = parameter.split()
    mol = (fields[0])  # Name of the molecule or bond
    k = float(fields[1])  # Harmonic force constant
    d_0 = float(fields[2])  # Dissociation energy
    a = float(fields[3])  # Morse width parameter
    r_eq = float(fields[4])  # Potential minimum distance

fcc = float(input(str("Input force convergence criterion (0.01): ")))
max_iter = int(input(str("Input maximum number of iterations (100): ")))


def kratzer_potential(r):
    return d_0 * ((r - r_eq) / r) ** 2


def force_kratzer(r):
    return d_0 * (2 * r_eq * (r - r_eq)) / (r ** 3)


def hessian_kratzer(r):
    h = d_0 * (6 * r_eq ** 2 - 4 * r * r_eq) / (r ** 4)
    return h, 1 / h


def harmonic_potential(r):
    return 0.5 * k * (r - r_eq) ** 2


def force_harmonic(r):  # First derivative
    return k * (r - r_eq)


def hessian_harmonic(r):  # Second derivative
    h = np.array(k)
    return h, 1 / h


def morse_potential(r):
    return d_0 * (1 - np.exp(-a * (r - r_eq))) ** 2


def force_morse(r):  # First derivative
    return 2 * a * d_0 * np.exp(-a * (r - r_eq)) * (1 - np.exp(-a * (r - r_eq)))


def hessian_morse(r):  # Second derivative
    h = d_0 * (2 * a ** 2 * np.exp(-2 * a * (r - r_eq)) - 2 * a ** 2 * np.exp(-a * (r - r_eq)) * (
            1 - np.exp(-a * (r - r_eq))))
    ih = 1 / h
    return h, ih


def conjugate_grad(H, b, x=None):
    n = np.array(b).size
    if not np.array(x):
        x = np.ones(n)
    r = np.array(H) * x - b
    p = -r
    rk_norm = r * r
    for k in range(2 * n):
        Hp = H * p
        alpha = rk_norm / (p * Hp)
        x += alpha * p
        r += alpha * Hp
        rk_new_norm = r * r
        beta = rk_new_norm / rk_norm
        rk_norm = rk_new_norm
        if rk_new_norm.all() < 10 ** -5:
            #print(k)
            break
        p = beta * p - r
    return x[0]


def hessian_bfgs(r0, force, iH0):  # Guess the Hessian
    p = -1 * force(r0) * iH0
    # alpha = 0.001
    alpha = line_search(morse_potential, force, r0, p)
    s = alpha * p
    r = r0 + s
    y = force(r) - force(r0)
    num = y * y
    denom = y * s
    if denom != 0:
        H = num / denom
        iH = 1 / H
    else:
        H, iH = 0, 0
    return H, iH


def line_search(f, g, x, p):
    a, b = 0.1, 0.1
    alpha = 0.001
    while f(x + alpha * p) > f(x) + a * alpha * np.dot(g(x), p):
        alpha *= b
    return alpha


def newton_rhapson(r, pot, force, hessian):
    x = [0]
    y = [r]
    E = pot(r)
    F = force(r)
    H, iH = hessian(r)
    print("STEP  DISTANCE   ENERGY       FORCE ")
    print("%3d  %8.4f    %5.3f   %8.3f " % (0, r, E, F))
    for i in range(1, max_iter):
        r -= F * iH
        E = pot(r)
        F = force(r)
        H, iH = hessian(r)
        print("%3d  %8.4f     %5.3f   %8.3f " % (i, r, E, F))
        x.append(i)
        y.append(r)
        if abs(F) < fcc or r < 0:
            break
    return x, y


def bfgs(r, pot, force, hessian):
    x = [0]
    y = [r]
    E = pot(r)
    F = force(r)
    H = hessian(r)
    #r = conjugate_grad(H, F, r)  # F * iH
    #r = r.astype(float)
    y.append(r)
    E = pot(r)
    F = force(r)
    H, iH = hessian(r)
    print("STEP  DISTANCE   ENERGY   FORCE ")
    print("%3d  %8.4f    %5.3f   %8.3f " % (1, r, E, F))
    for i in range(2, max_iter):
        r += F * iH
        E = pot(r).astype(float)
        F = force(r).astype(float)
        print(r)
        print(E)
        print(F)
        x.append(i)
        y.append(r)
        print("%3d  %8.4f     %5.3f   %8.3f " % (i, r, E, F))
        if abs(F) < fcc or r < 0:
            break
        H, iH = hessian_bfgs(r, force, iH)
    return x, y
