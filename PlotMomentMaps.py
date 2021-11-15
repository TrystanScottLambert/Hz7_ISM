###############################################
#                                             #
# Script to Plot moment maps of Hz7 for paper #
#                                             #
###############################################

import numpy as np 
import pylab as plt
import matplotlib
from matplotlib.patches import Ellipse
from astropy.io import fits 
from astropy.wcs import WCS
import plotBeam as pb
matplotlib.rcParams.update({'font.size': 20})   #Always run this line just to have tick mark sizes looking good. 
rms_val = 2.3e-2
levels = np.arange(3,6)*rms_val


def plotFits(fitsFileName,contourFitsFileName,showContours = True,**kwargs):
	fitsData,fitsHeader, fitsWCS = openFitsFile(fitsFileName)
	contourData,_,_ = openFitsFile(contourFitsFileName)
	fig = plt.figure()
	ax = fig.add_subplot()#projection=fitsWCS)
	image = ax.imshow(fitsData,cmap='coolwarm',**kwargs)

	if showContours == True:
		ax.contour(contourData, cmap='Greys_r', alpha=0.5, levels=levels)

	ax.text(155.15,138.05,'1',fontsize=15)
	ax.text(133.11,146.72,'2',fontsize=15)
	ax.text(184.93,120.25,'3',fontsize=15)
	ax.text(182.07,141.55,'4',fontsize=15)
	ax.text(171.74,156.30,'5',fontsize=15)
	plt.colorbar(image,label=fitsHeader['BUNIT'])
	pb.drawBeam(fitsHeader,ax)
	prettifyPlot('RA','Dec')
	plt.xlim(100,200)
	plt.ylim(170,100)
	plt.show()

def openFitsFile(fitsFileName):
	hdu = fits.open(fitsFileName)
	fitsData = hdu[0].data[0][0] 
	fitsHeader = hdu[0].header 
	fitsWCS = WCS(fitsHeader,naxis=2)
	return fitsData,fitsHeader,fitsWCS

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)

if __name__ == '__main__':
	plotFits('data/HZ7_Collapsed.fits','data/HZ7_Collapsed.fits')
	plotFits('data/HZ7_mom_mask2.weighted_coord.fits','data/HZ7_Collapsed.fits',vmin=-100,vmax=100)
	plotFits('data/HZ7_mom_mask2.weighted_dispersion_coord.fits','data/HZ7_Collapsed.fits')