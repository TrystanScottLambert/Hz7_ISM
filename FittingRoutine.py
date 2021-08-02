#######################################
#
# Sersic Fitting Routine with errors
#
#######################################

import numpy as np 
import pylab as plt
from scipy.optimize import curve_fit 

# use f_min
# covariance not necessarily a good representation of the uncertainty. 
# for uncertainties using the lazy monte carlo method for uncertainties. (LOOK INTO THIS)


x = np.array([-1.768, -1.56 , -1.352, -1.144, -0.936, -0.728, -0.52 , -0.312,-0.104,  0.104,  0.312,  0.52 ,  0.728,  0.936,  1.144,  1.352,1.56 ,  1.768])
y = np.array([0.05871114, 0.05070998, 0.07717734, 0.08549176, 0.1420571 ,0.26898651, 0.51840227, 0.95852468, 1.77334112, 1.77334112,0.95852468, 0.51840227, 0.26898651, 0.1420571 , 0.08549176,0.07717734, 0.05070998, 0.05871114])
yerr = np.array([0.04661739, 0.04947074, 0.05288143, 0.0578848 , 0.06364147,0.07204962, 0.08585704, 0.11047333, 0.19343844, 0.19343844,0.11047333, 0.08585704, 0.07204962, 0.06364147, 0.0578848 ,0.05288143, 0.04947074, 0.04661739])

x,y,yerr = x[9:],y[9:],yerr[9:]

def sersic(R, Ie, Re, m):
	bm = 2.0*m - 0.324
	return Ie * np.exp(-bm * ((R/Re)**(1.0/m) - 1.0))

def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

popt, pcov = curve_fit(sersic,x,y,sigma=yerr)
perr= np.sqrt(np.diag(pcov))

nstd = 1. # to draw 5-sigma intervals
popt_up = popt + nstd * perr
popt_dw = popt - nstd * perr

fit = sersic(x, *popt)
fit_up = sersic(x, *popt_up)
fit_dw = sersic(x, *popt_dw)

print('fit parameter 1-sigma error')
print('———————————–')
for i in range(len(popt)):
	print(str(popt[i])+' +- '+str(perr[i]))

xplot = np.linspace(0,2,1000)
plt.errorbar(x, y, yerr=yerr, hold=True, ecolor='k', fmt='none', label='data')
plt.plot(xplot,sersic(xplot,*popt),lw=2)
plt.fill_between(x, fit_up, fit_dw, alpha=.25, label='5-sigma interval')
plt.show()