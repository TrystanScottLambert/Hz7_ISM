################################
#
# Circular Class 
#
################################

import numpy as np 
import pylab as plt
from photutils.aperture import CircularAperture 
from photutils.aperture import aperture_photometry
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.rcParams.update({'font.size': 20})   #Always run this line just to have tick mark sizes looking good. 
from acstools import acszpt 

# function that will remove non needed dimensions of an array. In particular the stokes and frequency dimesntions of a radio cube
def reduceArrayTo2D(array):
	if len(array.shape) > 2:
		new_array = array[0]
		while len(new_array.shape) > 2:
			new_array = new_array[0]
		return new_array
	else:
		return array

def convertInstrumentFluxToPhysical(flux_array,filt,correction,date='2015-05-07'): # correction from https://www.stsci.edu/hst/instrumentation/acs/data-analysis/aperture-corrections
	flux_at_infinity = flux_array / correction
	q = acszpt.Query(date=date,detector='WFC',filt=filt)
	filter_zpt = q.fetch()

	F_lambda = flux_at_infinity * filter_zpt['PHOTFLAM']
	m_st = -2.5*np.log10(flux_at_infinity) + filter_zpt['STmag'][0].value 
	m_ab = -2.5*np.log10(flux_at_infinity) + filter_zpt['ABmag'][0].value 
	m_vega = -2.5*np.log10(flux_at_infinity) + filter_zpt['VEGAmag'][0].value 
	return m_st,m_ab,m_vega

def getFilterFromName(filename):
	filt = filename.split('_f')[1].split('w_')[0]
	return 'F'+filt+'W'


class Image:
	def __init__(self,filename):
		self.hdu = fits.open(filename)
		self.data = reduceArrayTo2D(self.hdu[0].data)
		self.header = self.hdu[0].header
		self.wcs = WCS(self.header,naxis=2)
		self.arcsecsInAPixel = self.header['CDELT2'] * 3600

	# function to work out a circular distribution of circles based on ra,dec and radius
	def setCenters(self,center,radius):
		centers = [center]
		for i in range(6):
			centers.append((2*radius*np.cos((i*np.pi)/3)+center[0],2*radius*np.sin((i*np.pi)/3)+center[1]))
		return centers

	def getCircularApetures(self,center,radius):
		positions = self.setCenters(center,radius)
		apetures = [CircularAperture(pos,r=radius) for pos in positions]
		return apetures

	def getCirclePlottingInformation(positions,radius):
		x = np.array([pos[0] for pos in positions])
		y = np.array([pos[1] for pos in positions])
		r = np.ones(len(x))*radius
		return x,y,r

	def calculateFluxes(self,center,pixelRadius):
		pixelCenter = self.wcs.world_to_pixel_values(center[0],center[1])
		pixelCenter = (float(pixelCenter[0]),float(pixelCenter[1]))
		apertures = self.getCircularApetures(pixelCenter,pixelRadius)
		error = np.sqrt(np.abs(self.data))
		
		fluxes = []
		fluxErrors = []
		for aperture in apertures:
			fluxErrors.append(list(aperture_photometry(self.data,aperture,error=error)['aperture_sum_err'])[0])
			fluxes.append(list(aperture_photometry(self.data,aperture,error=error)['aperture_sum'])[0])
		
		self.fluxes = np.array(fluxes)
		self.fluxErrors = np.array(fluxErrors)
		self.pixelCenter = pixelCenter

if __name__ == '__main__':
	center = (149.87692083333334,2.1340558333333335)
	radio = Image('data/HZ7_Collapsed.fits')
	#radio = Image('/home/trystan/Downloads/HZ7_CONTINUUM_250.fits')
	infile = 'data/ASCImages/deresolved/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved_deresolved.fits'
	optical = Image(infile)
	radio.calculateFluxes(center,4)
	optical.calculateFluxes(center,4)
	optical_fluxes_corrected = convertInstrumentFluxToPhysical(optical.fluxes,filt='F435W',correction=0.7)
	for i in range(len(radio.fluxes)):
		plt.errorbar(optical.fluxes[i],radio.fluxes[i],xerr=optical.fluxErrors[i],yerr=radio.fluxErrors[i],label=f'{i+1}',fmt='o',capsize=10)

	plt.grid()
	plt.ylabel('[CII]')
	plt.xlabel('UV')
	plt.legend()
	plt.show()