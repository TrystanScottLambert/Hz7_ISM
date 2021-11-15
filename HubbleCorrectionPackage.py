######################
#
# Hubble Corrections
#
######################

import numpy as np 
import pylab as plt

infile = 'WFC_Aperture_Corrections.dat'
f = open(infile)
pixels = np.array([1,2,3,4,5,6,7,8,9,10,20,40])

values = []
for line in f:
	val = line.split('\t')
	values.append(val[1:])
values = np.array(values)
values = values.astype(float)

for value in values:
	plt.plot(pixels,value)
