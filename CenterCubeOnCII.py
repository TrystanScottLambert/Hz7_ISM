#######################################################
#
# Script to center HZ7_Combined.fits on [CII] line
#
#######################################################

import numpy as np 
from astropy.io import fits 
from astropy import units as u 
from spectral_cube import SpectralCube
from astropy.wcs import WCS
import pylab as plt

ciiRestFreq = 1897420620253.1646
redshift = 5.2532
centerRadioFreq = 2.517843e5 #km/s
ciiRedshifted = ciiRestFreq*(1./(1+redshift))

infile = 'data/HZ7_Combined.fits'
hdu = fits.open(infile)
currentRestFreq = hdu[0].header['RESTFRQ'] 


freqDummy = np.arange(0,hdu[0].data.shape[1])
decDummy = np.arange(0,hdu[0].data.shape[1])
raDummy = np.arange(0,hdu[0].data.shape[1])


#fits.writeto('data/HZ7_Centered.fits',hdu[0].data,hdu[0].header,overwrite=True)

meanVels=[]
tolerance = 0.00001
meanfreq = np.mean([ciiRedshifted,currentRestFreq])
meanVel=10

def calculateVelforFreq(frequency):
	cube = SpectralCube.read(hdu)
	cube_vel = cube.with_spectral_unit(u.km/u.s, rest_value = frequency*u.Hz, velocity_convention='optical')
	infile = 'delete_this_is_a_test.fits'
	cube_vel.write(infile,overwrite=True)
	hdu_vel = fits.open(infile)
	wcs = WCS(hdu_vel[0].header,naxis=3)
	_,_,currentVelRange = wcs.pixel_to_world_values(raDummy,decDummy,freqDummy)
	currentVel = np.mean(currentVelRange[52:72])  # fitted guassian +- standard deviation
	return currentVel

left = ciiRedshifted
right = currentRestFreq
val = 10
while np.abs(val) > tolerance:
	midpoint = np.mean([left,right])
	val = calculateVelforFreq(midpoint)
	if val > 0:
		right = midpoint
	elif val < 0:
		left = midpoint


optimumFrequency = midpoint

hdu[0].header['RESTFRQ'] = optimumFrequency
optimumredshift =  (-optimumFrequency+ciiRestFreq)/optimumFrequency  
fits.writeto('data/HZ7_Centered.fits',hdu[0].data,hdu[0].header,overwrite=True)


