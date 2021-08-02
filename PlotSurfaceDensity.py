#############################
#                           #
# Plot Surface Density Plot #
#                           #
#############################

import numpy as np 
import pylab as plt
from astropy.io import fits 
from astropy.wcs import WCS
from astropy import modeling 
from astropy.convolution import convolve, Box1DKernel 
from astropy.modeling.models import Sersic1D
from RegscorePy import * 
from scipy.optimize import curve_fit 
from astropy.visualization import simple_norm

class RadioCube:
	def __init__(self,infile):
		self.hdu = fits.open(infile)
		self.data = self.hdu[0].data[0][0]
		self.header = self.hdu[0].header 
		self.wcs = WCS(self.header)
		self.degreePerPixel = np.abs(float(self.header['CDELT1']))
		self.arcsecondPerPixel = self.degreePerPixel*3600 #arcseconds
		self.BMAJ = self.header['BMAJ']
		self.BMIN = self.header['BMIN']
		self.beamArea = 2*(np.pi*self.BMAJ*self.BMIN/(8.0*np.log(2.0))) # arcsec^2 
		self.beamAreaInPixels = self.beamArea/(self.arcsecondPerPixel**2)
		self.beamCorrectedData = self.data/self.beamAreaInPixels
		self.arcsecondAreaPerPixel = self.arcsecondPerPixel**2


	def generateAnuli(self,maximumArcsecond,binwidthInArcseconds):
		binwidthInPixels = round(binwidthInArcseconds/self.arcsecondPerPixel)
		maximumPixel = np.ceil(maximumArcsecond/self.arcsecondPerPixel)
		self.anuli = np.arange(0,maximumPixel+binwidthInPixels,binwidthInPixels)

	def getSumAreaUncertainty(self,centerPix,maximumArcsecond,binwidthInArcseconds,RMS):
		self.generateAnuli(maximumArcsecond,binwidthInArcseconds)
		self.sums,self.areas,self.uncertainties = self.sumRegions(centerPix,RMS)

	def sumRegion(self,center,radiusMin,radiusMax,RMS,error=True):
		xs = []
		ys = []
		for i in range(len(self.data)):
			for j in range(len(self.data[0])):
				xs.append(i)
				ys.append(j)
		xs = np.array(xs)
		ys = np.array(ys)
		rs = np.sqrt((center[0]-xs)**2+(center[1]-ys)**2)

		cut = np.where((rs<=radiusMax) & (rs>=radiusMin))[0]
		xs = xs[cut]
		ys = ys[cut]
		rs = rs[cut]

		val = []
		for i in range(len(xs)):
			val.append(self.beamCorrectedData[ys[i]][xs[i]])

		sumVal = np.nansum(val)
		annuliArea = len(xs)*self.arcsecondAreaPerPixel

		if error==False:
			return sumVal, len(xs)*arcsecArea
		elif error == True:
			return sumVal,annuliArea,RMS*(np.sqrt(annuliArea/self.beamArea))/annuliArea

	def sumRegions(self,center,RMS):
		sums=[]
		areas=[]
		uncertainties=[]
		for i in range(len(self.anuli)-1):
			currentSum,currentArea,currentUncertainty = self.sumRegion(center,self.anuli[i],self.anuli[i+1],RMS,error=True)
			sums.append(currentSum)
			areas.append(currentArea)
			uncertainties.append(currentUncertainty) 
		return np.array(sums),np.array(areas),np.array(uncertainties)


if __name__ == '__main__':
	infile = 'data/HZ7_Collapsed.fits'
	rc = RadioCube(infile)
	rc.getSumAreaUncertainty((155,139),2,0.2,6e-5)
