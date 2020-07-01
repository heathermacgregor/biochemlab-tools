#!/usr/bin/env python
import numpy as np
from itertools import combinations_with_replacement as comb
from time import time
x = 1063.9
y = 1064

aa_name = ["Gly", "Ala", "Ser", "Pro",  "Val", "Thr", "Cys", "Leu or Ile", "Asn", "Asp or isoAsp", "Gln", "Lys", "Glu", "Met", "His", "Phe", "Arg", "Tyr", "Trp", "4-hydroxyPro", "5-hydroxyLys", "6-N-methylLys", "g-carboxyGlu", "selenoCys", "phosphoSer", "phosphoThr", "phosphoTyr", "s-N-methylArg", "6-N-acetylLys", "Glu g-methyl ester", "Ornithine", "Citrulline", "3-methylHis", "N,N,N-trimethylLys", "N-acetylAla", "3-sulfinoAla", "N-acetylCys", "pyroGlu", "N-acetylGly", "Met sulfoxide", "Met sulfone", "N-acetylSer", "N-acetylThr", "Kynurenine", "Tyr O-sulfate", "cystine"]

aa_mass = [57.0215, 71.0371, 87.0320, 97.0528, 99.0684, 101.0477, 103.0092, 113.0841, 114.0429, 115.0269, 128.0586, 128.0950, 129.0426, 131.0405, 137.0589, 147.0684, 156.1011,163.0633, 186.0793, 113.048, 144.089, 142.110, 173.032, 150.954, 166.998, 181.014, 243.029, 170.116, 170.105, 143.058, 114.079, 157.085, 151.074, 170.141, 113.047, 134.999, 145.019, 111.032, 99.032, 147.035, 163.030, 129.042, 143.058, 190.047, 243.020, 222.013]


aa_dict = {aa_name[i] : aa_mass[i] for i in range(len(aa_name))}
combinations = comb(aa_dict, 10)
print('Finished Combinations')



def func(tup):
    vals = [aa_dict[k] for k in tup]
    return x <= sum(vals) <= y

comb_filt = filter(func, combinations)
print('Finished Filtering')

valid_comb = {}
part = 0
a = time()

for tup in comb_filt:
    summ = 0
    part += 1
    print('Working on part {}'.format(part))
    
    for aa in tup:
        summ += aa_dict[aa]
    tup_fancy = ' '.join(tup)
    valid_comb[tup_fancy] = summ
    #print('Sum is {}'.format(summ))
    
    if part % 5 == 0:
        #print('Part {} done'.format(part))
        time_diff = time() - a
        a = time()
        print('Time since last check: {}'.format(time_diff))
        

print('Finished Summing')


f = open('sequence.txt', "w")
for val in valid_comb:
    text = '{}: {} \n'.format(str(val), str(valid_comb[val]))
    f.write(text)
    
f.close()
print('Finished Writing')
