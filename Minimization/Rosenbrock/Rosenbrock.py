#!/usr/bin/env python3
########################################################################################################################
# Functions for Newton-Rhapson minimizer of the Rosenbrock function (Rosenbrock_NR.py) #################################
########################################################################################################################

import sympy as sym
import numpy as np
from numpy.linalg import inv                                                                           # Inverse hessian

# Rosenbrock Function:
# Weisstein, Eric W. "Rosenbrock Function." From MathWorld--A Wolfram Web Resource.
# https://mathworld.wolfram.com/RosenbrockFunction.html
def rosenbrock(x, y):
    return (1 - x) ** 2 + 100 * (y - x ** 2) ** 2


def gradient(a, b):
    x, y = sym.symbols('x y')
    f = rosenbrock(x, y)
    if a == 'x' and b == 'y':                                                                       # Solve analytically
        df_dx, df_dy = sym.diff(f, x), sym.diff(f, y)
    else:                                                                                            # Solve numerically
        df_dx, df_dy = sym.diff(f, x).subs(x, a).subs(y, b), sym.diff(f, y).subs(x, a).subs(y, b)
    return np.array((df_dx, df_dy))


def hessian(a, b):
    x, y = sym.symbols('x y')
    h = []
    for i in gradient('x', 'y'):
        if a == 'x' and b == 'y':                                                                   # Solve analytically
            j = [sym.diff(i, x), sym.diff(i, y)]
        else:                                                                                        # Solve numerically
            j = [sym.diff(i, x).subs(x, a).subs(y, b), sym.diff(i, y).subs(x, a).subs(y, b)]
        h.append(j)
    h = np.array(h).astype(float)                                                                              # Hessian
    ih = inv(h)                                                                                        # Inverse hessian
    return h, ih
