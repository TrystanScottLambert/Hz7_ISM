########################################
#
# 2D Sersic fitting to moment map 
# Trystan Lambert 
#
########################################

import numpy as np 
import pylab as plt
from astropy.io import fits 
from astropy import modeling
import scipy.optimize as opt


infile = 'data/moment0.fits'
moment0 = fits.open(infile)
data = moment0[0].data[0][0]   # open file and strip data into the file
header = moment0[0].header
PixSize = header['CDELT2']*3600
BMAJ = header['BMAJ']*3600
BMIN = header['BMIN']*3600
AreaBeam = 2*(np.pi*BMAJ*BMIN/(8.0*np.log(2.0))) #arcsec^2
AreaBeamPixels = AreaBeam / (PixSize**2)
data = data/AreaBeamPixels



def sersicProfile2D((x,y),amplitude,r_eff,n,x0,y0,ellip,theta):
	sersicModel = modeling.models.Sersic2D(amplitude = amplitude, r_eff = r_eff, 
		n=n, x_0=x0, y_0=y0,ellip=ellip, theta=theta)
	val = sersicModel(x,y)
	return val.ravel()

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


x,y = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))

initial_guess = (5.705969e-06,25,1,155,139,1,1)
initial_guess = (5.705969e-06,25,0.6,155,139,.5,-1)
#initial_guess = (5.705969e-06,155,139,20,20,1,0)

popt,pcov = opt.curve_fit(twoD_Gaussian,(x,y),data.ravel(),p0=initial_guess)
popt,pcov = opt.curve_fit(sersicProfile2D,(x,y),data.ravel(),p0=initial_guess)


data_fitted=sersicProfile2D((x,y),*popt)
#data_fitted = twoD_Gaussian((x,y),*popt)
plottingData = data_fitted.reshape(data.shape[0],data.shape[1])

plt.imshow(data)
plt.title('Moment0 Map with 2D Sersic Countours')
plt.contour(x,y,data_fitted.reshape(data.shape[0],data.shape[1]))
plt.show()

plt.subplot(131)
plt.title('Moment 0 Map')
plt.imshow(data,vmin=-2.14e-6,vmax=5.8e-6)
plt.subplot(132)
plt.title('Sersic Fit')
plt.imshow(plottingData,vmin=-2.14e-6,vmax=5.8e-6)
plt.subplot(133)
plt.title('Moment 0 Map - Sersic Fit')
plt.imshow(data-plottingData,vmin=-2.14e-6,vmax=5.8e-6)
plt.show()
