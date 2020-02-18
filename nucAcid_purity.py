#!/usr/bin/env python

import os, sys, re
import math
import numpy as np
import matplotlib.pyplot as plt
import itertools

# Checks the purity of UV-Vis spectral data (230 - 320 nm) of a DNA sample &
# produces a plot of the entire spectrum.
# Input files are assumed to be .txt output files from Shimadzu UV spectrometers,
# which contain two lines that do not include numerical data:
# "RawData"
# "Wavelength nm."	"Abs."
# If your data files have a different format, modify the code appropriately.

spectrum = input("ENTER TXT file containing spectral data: ")
sample_name = input("ENTER sample name: ")
d_factor = input("ENTER dilution factor: ")
Wavelength = []
Absorbance = []

with open(spectrum, 'r+') as f:
    for line in itertools.islice(f, 2, None):
        edit = line.split('\t')
        #print(edit[0])
        Wavelength.append(float(edit[0]))
        Absorbance.append(float(edit[1]))

A280 = Absorbance[Wavelength.index(float("280.00"))]
A260 = Absorbance[Wavelength.index(float("260.00"))]
A320 = Absorbance[Wavelength.index(float("320.00"))]

# Concentration (µg/ml) = (A260 reading – A320 reading) × dilution factor × 50µg/ml
dna_conc = (A260 - A320)*(float(d_factor))*(50)

# DNA purity (A260/A280) = (A260 reading – A320 reading) ÷ (A280 reading – A320 reading)
dna_purity = (A260 - A320)/(A280 - A320)
dna_purity_01 = A260/A280
dna_purity_02 = A260/A320
print(str("----------------------------"))
print(str("A260 \t A280 \t A320"))
print(str("----------------------------"))
print(str(str(A260) + "\t" + str(A280) + "\t" + str(A320)))
print(str("----------------------------"))
print(str("\n"))

print(str("Sample concentration: ") + str(dna_conc) + str(" µg/ml"))
print(str("RNA contamination (A260/A280): \t") + str(dna_purity_01))
print(str("Organics and/or salts contamination (A260/A320): ") + str(dna_purity_02))
print(str("Sample purity: ") + str(dna_purity))
print(str("\n"))

if 1.7 <= dna_purity <= 2.0:
    print(str("Congratulations! Your sample looks good."))
else:
    if dna_purity_01 > 1.8:
        print(str("Oh no! Your sample appears to be contaminated by RNA."))
    elif dna_purity_02 < 1.5:
        print(str("Oh no! Your sample appears to be contaminated by chaotropic salts."))
    else:
        print(str("Oh no! Your sample appears to be contaminated."))


plt.plot(Wavelength, Absorbance, 'b--')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorbance")
# plt.show()
plt.savefig(str("./" + str(sample_name) + "_UV_Spectrum.png"))
