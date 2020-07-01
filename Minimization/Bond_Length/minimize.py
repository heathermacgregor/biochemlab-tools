#!/usr/bin/python
#######################################################################################################################
# Finds the minimum-energy C-O bond length using a method and potential function of choice.
#######################################################################################################################
import sys
import matplotlib.pyplot as plt
import numpy as np
import time
import minimize_dict as m

if len(sys.argv) != 3:
    sys.exit("Usage: minimize.py param_file initial_r")

r_ini = float(sys.argv[2])


def run(pot, method, m, f, g, h, x0):
    print("\n Optimizing using " + str(method) + " Method and " + str(pot) + " Potential:")
    start = time.time()
    x, y = m(x0, f, g, h)
    end = time.time()
    t = end - start
    return x, y, t


choose_potential = int(input(str("Choose a potential: \n " +
                                 "\t (1) Harmonic Oscillator \n " +
                                 "\t (2) Morse \n " +
                                 "\t (3) Kratzer \n " +
                                 "\t (4) All \n ")))

def select_run(choose_potential):
    if choose_potential != 4:
        fig, ax = plt.subplots(2)
        if choose_potential == 1:
            pot, f, g, h, ff = "Harmonic", m.harmonic_potential, m.force_harmonic, m.hessian_harmonic, "K_NR-v-BFGS_"
        if choose_potential == 2:
            pot, f, g, h, ff = "Morse", m.morse_potential, m.force_morse, m.hessian_morse, "K_NR-v-BFGS_"
        if choose_potential == 3:
            pot, f, g, h, ff = "Kratzer", m.kratzer_potential, m.force_kratzer, m.hessian_kratzer, "K_NR-v-BFGS_"
        xNR, yNR, tNR = run(pot, "Newton-Rhapson", m.newton_rhapson, f, g, h, r_ini)
        xBFGS, yBFGS, tBFGS = run(pot, "BFGS Update", m.bfgs, f, g, h, r_ini)
        rs = np.arange(min(sum([yNR, yBFGS], [])), max(sum([yNR, yBFGS], [])), 0.001)
        ax[0].plot(rs, [f(i) for i in rs.astype(np.float)], 'k-',
                   yNR, [f(i) for i in yNR], 'co')
        ax[1].plot(rs, [f(i) for i in rs.astype(np.float)], 'k-',
                   yBFGS, [f(i) for i in yBFGS], 'm*')

        ax[0].set(title=str(pot + " NR"))
        ax[1].set(title=str(pot + " BFGS"))

    if choose_potential == 4:
        fig, axs = plt.subplots(2, 3)

        pot = "Harmonic"
        f, g, h = m.harmonic_potential, m.force_harmonic, m.hessian_harmonic
        xHarmNR, yHarmNR, tHarmNR = run(pot, "Newton-Rhapson", m.newton_rhapson, f, g, h, r_ini)
        xHarmBFGS, yHarmBFGS, tHarmBFGS = run(pot, "BFGS Update", m.bfgs, f, g, h, r_ini)
        print(tHarmNR)
        print(tHarmBFGS)
        rs = np.arange(min(sum([yHarmNR, yHarmBFGS], [])), max(sum([yHarmNR, yHarmBFGS], [])), 0.001)
        axs[0, 0].plot(rs, [m.harmonic_potential(i) for i in rs.astype(np.float)], 'k-',
                       yHarmNR, [m.harmonic_potential(i) for i in yHarmNR], 'co')
        axs[1, 0].plot(rs, [m.harmonic_potential(i) for i in rs.astype(np.float)], 'k-',
                       yHarmBFGS, [m.harmonic_potential(i) for i in yHarmBFGS], 'm*')
        axs[0, 0].set(title=str(pot + " NR"))
        axs[1, 0].set(title=str(pot + " BFGS"))

        pot = "Morse"
        f, g, h = m.morse_potential, m.force_morse, m.hessian_morse
        xMorseNR, yMorseNR, tMorseNR = run(pot, "Newton-Rhapson", m.newton_rhapson, f, g, h, r_ini)
        xMorseBFGS, yMorseBFGS, tMorseBFGS = run(pot, "BFGS Update", m.bfgs, f, g, h, r_ini)
        print(tMorseNR)
        print(tMorseBFGS)
        rs = np.arange(min(sum([yMorseNR, yMorseBFGS], [])), max(sum([yMorseNR, yMorseBFGS], [])), 0.001)
        axs[0, 1].plot(rs, [m.morse_potential(i) for i in rs.astype(np.float)], 'k-',
                       yMorseNR, [m.morse_potential(i) for i in yMorseNR], 'co')
        axs[1, 1].plot(rs, [m.morse_potential(i) for i in rs.astype(np.float)], 'k-',
                       yMorseBFGS, [m.morse_potential(i) for i in yMorseBFGS], 'm*')
        axs[0, 1].set(title=str(pot + " NR"))
        axs[1, 1].set(title=str(pot + " BFGS"))

        pot = "Kratzer"
        f, g, h = m.kratzer_potential, m.force_kratzer, m.hessian_kratzer
        xKratzerNR, yKratzerNR, tKratzerNR = run(pot, "Newton-Rhapson", m.newton_rhapson, f, g, h, r_ini)
        xKratzerBFGS, yKratzerBFGS, tKratzerBFGS = run(pot, "BFGS Update", m.bfgs, f, g, h, r_ini)
        print(tKratzerNR)
        print(tKratzerBFGS)
        rs = np.arange(min(sum([yKratzerNR, yKratzerBFGS], [])), max(sum([yKratzerNR, yKratzerBFGS], [])), 0.001)
        axs[0, 2].plot(rs, [m.kratzer_potential(i) for i in rs.astype(np.float)], 'k-',
                       yKratzerNR, [m.kratzer_potential(i) for i in yKratzerNR], 'co')
        axs[1, 2].plot(rs, [m.kratzer_potential(i) for i in rs.astype(np.float)], 'k-',
                       yKratzerBFGS, [m.kratzer_potential(i) for i in yKratzerBFGS], 'm*')
        axs[0, 2].set(title=str(pot + " NR"))
        axs[1, 2].set(title=str(pot + " BFGS"))

        ff = "HMK_NR-v-BFGS_"
    elif 5 < choose_potential or choose_potential < 0 or not int(choose_potential):
        print("Stop that! I can't interpret other inputs!! (╬ ಠ益ಠ)")
    for ax in fig.get_axes():
        ax.set(xlabel="r (A)", ylabel="Energy")
    plt.tight_layout()
    fig.savefig(str(ff + str(r_ini) + ".png"))
    plt.show()


select_run(choose_potential)
