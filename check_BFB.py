import numpy as np
import math
from copy import copy,deepcopy
#import matplotlib.pyplot as plt

#------ checkBFB.py -----------------------------
#
# This script verifies that a new code (exp2) gives 
# BFB result with a reference code (exp1) by
# comparing the h fields. 
#
#------------------------------------------------

#------ INPUT by user ---------------------------

outputdir1="/aos/home/lemieux/BFBtest/refcode/output/"
outputdir2="/aos/home/lemieux/BFBtest/newcode/output/"
exp1="10"
exp2="11"
date="2002_01_01_12_00"

#------------------------------------------------

#exp1 - reference
file1=outputdir1+"h"+date+"."+exp1
h1 = np.genfromtxt(file1, dtype=None)

#exp2 - new
file2=outputdir2+"h"+date+"."+exp2
h2 = np.genfromtxt(file2, dtype=None)

hdiff=h2-h1

metric=math.sqrt(np.sum(hdiff**2))
print("The metric for the BFB check is = ")
print(metric)
print("It is BFB if metric = 0.")

#dh=1e-6
#plt.pcolor(hdiff, cmap='bwr', vmin=-dh, vmax=dh)
#plt.colorbar()
#plt.show()
