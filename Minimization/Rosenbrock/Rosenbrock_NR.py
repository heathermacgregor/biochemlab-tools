#!/usr/bin/env python3
########################################################################################################################
# Newton-Rhapson minimizer of the Rosenbrock function ##################################################################
########################################################################################################################

import os
import sys
import numpy as np
from Rosenbrock import rosenbrock, gradient, hessian

parameter_file = open(sys.argv[1])
for line in parameter_file.readlines():
    parameter = line.split(', ')
    x_0 = float(parameter[0])                                                                               # Initial x
    y_0 = float(parameter[1])                                                                               # Initial y
    step_size = float(parameter[2])                                                                         # Step size
    smaller_step_size = step_size * 0.01

try:
    os.system("rm Rosenbrock_NR.log")
    print("Removed old logfile.")
except IOError:
    print("No old logfile found.")

print('====================================================================================================================================')
print('STEP    X        Y        Z        GRADIENT                                 HESSIAN                                                 ')
print('====================================================================================================================================')

logfile = open("Rosenbrock_NR.log", "a+")
def newton_rhapson(x, y, f):
    for i in range(1, 10 ** 6):                                                               # Stop looking eventually
        x_n, y_n = x[-1], y[-1]
        if i > 20:
            if np.std([f[-10:]]) < 10 ** (-15):
                if f[-1] <= 1 * 10 ** (-15):                       # Convergence criteria met when close to the minimum
                    break
        g_n = gradient(x_n, y_n)
        h_n, ih_n = hessian(x_n, y_n)
        if i > 1000 and f_n < 1 * 10 ** (-14) or np.std([f[-100:]]) < 10 ** (-14):
            step_n = smaller_step_size * np.dot(ih_n, g_n)               # Take smaller steps when close to the minimum
        else:
            step_n = step_size * np.dot(ih_n, g_n)
        xy_n = np.array((x_n, y_n)) - step_n
        x.append(xy_n.item(0))
        y.append(xy_n.item(1))
        f_n = rosenbrock(x_n, y_n)
        f = np.append(f, f_n).astype(float)
        log = '%15.14f,  %15.14f,  %15.14f, %15.14f, %15.14f, %15.14f' % (x_n, y_n, f_n, x_n - g_n[0], x_n - g_n[1], rosenbrock(g_n[0], g_n[1]))
        logfile.write(log)
        logfile.write('\n')                                                              # Write results to a text file
        console = '%7d %7.6f %7.6f %7.6f %s     %s      %s      ' % (i, x_n, y_n, f_n, g_n, ih_n[0], ih_n[1])
        print(console)                                                                   # Show progress in the console


x_1, y_1 = [x_0], [y_0]
f_1 = np.array([rosenbrock(x_0, y_0)])

newton_rhapson(x_1, y_1, f_1)
logfile.close()
