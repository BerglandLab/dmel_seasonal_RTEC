#### calculates fisher exact pvalues using Jamie's script for simulated data
## Oct 31, 2016
#argv = ["script","/Users/heathermachado/nescent_melCA/simulated_data/datasets/data.tmp", "testout.txt"]

import math
import numpy as np
from sys import argv

#defines the factorial, can change to use sterling's approx if desired
def st(n):
    l=math.log(math.factorial(n))
    return l

# log probability for a pair takes r2, D2 and r, D
def log_p(r2, D2, r, D):
    return (st(D-D2)-st(r-r2)-st(D-D2-(r-r2)))+(st(D2)-st(r2)-st(D2-r2))-(st(D)-st(r)-st(D-r))

def fisher_exact_p_val(r1, D1, r2, D2):
    r=r1+r2
    D=D1+D2
    log_probs=[log_p(i, D2, r, D) for i in range(r2, min(r+1, D2+1))]
    probs=[math.exp(l) for l in log_probs]
    p_val=sum(probs[1:])+probs[0]*np.random.random()
    if p_val>1.0:
        print('error: ', p_val)
    return p_val




###################### reading and writing the files
filein = argv[1]
fileout = argv[2]

data1 = filein
fout = open(fileout,'w')
print >> fout, "chrom pos sFreq fFreq sDP fDP dp.total sAlt fAlt alt.total sRef fRef greater less"

count = 0
with open(data1) as f:
    for line in f:
        #count += 1
        line2 = line.split(' ')
        chrom=line2[0]
        pos=line2[1]
        r1=int(line2[2])
        r2=int(line2[3])
        D1=int(line2[4])
        D2=int(line2[5])
        p_lesser = fisher_exact_p_val(r1, D1, r2, D2)
        p_greater = fisher_exact_p_val(r2, D2, r1, D1)
        #print p_greater
        #print p_lesser
        sFreq = float(r1)/float(D1)
        fFreq = float(r2)/float(D2)
        dp_total = D1+D2
        sRef = D1-r1
        fRef = D2-r2
        alt_total = r1+r2
        
        my_array = [chrom, pos, sFreq, fFreq, D1, D2, dp_total, r1, r2, alt_total, sRef, fRef, p_greater, p_lesser]
        print >> fout, "\t".join(str(i) for i in my_array)

fout.close()