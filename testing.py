import numpy as np
from pylab import *
from scipy.optimize import curve_fit

def func(x, a, b, c):

	return a * x *x + b*x + c

# test data and error
x = np.linspace(-10, 10, 100)
y0 = -0.07 * x * x + 0.5 * x + 2.
noise = np.random.normal(0.0, 1.0, len(x))
y = y0 + noise

# curve fit [with only y-error]
popt, pcov = curve_fit(func, x, y, sigma=1./(noise*noise))
perr = np.sqrt(np.diag(pcov))

#print fit parameters and 1-sigma estimates
print('fit parameter 1-sigma error')
print('———————————–')
for i in range(len(popt)):
	print(str(popt[i])+' +- '+str(perr[i]))

# prepare confidence level curves
nstd = 5. # to draw 5-sigma intervals
popt_up = popt + nstd * perr
popt_dw = popt - nstd * perr

fit = func(x, *popt)
fit_up = func(x, *popt_up)
fit_dw = func(x, *popt_dw)

#plot
fig, ax = plt.subplots(1)
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['font.size']= 20
errorbar(x, y0, yerr=noise, xerr=0, hold=True, ecolor='k', fmt='none', label='data')

xlabel('x', fontsize=18)
ylabel('y', fontsize=18)
title('fit with only Y-error', fontsize=18)
plot(x, fit, 'r', lw=2, label='best fit curve')
plot(x, y0, 'k-', lw=2, label='True curve')
ax.fill_between(x, fit_up, fit_dw, alpha=.25, label='5-sigma interval')
legend(loc='lower right',fontsize=18)
show()