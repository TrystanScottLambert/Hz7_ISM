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
matplotlib.rcParams.update({'font.size': 2})   #Always run this line just to have tick mark sizes looking good. 
rms_val = 2.3e-2
levels = np.arange(3,6)*rms_val
import matplotlib as mpl  

#making the plots look nice
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 

def plotFits(fitsFileName, contourFitsFileName, outFileName, showContours = True, **kwargs):
	fitsData,fitsHeader, fitsWCS = openFitsFile(fitsFileName)
	contourData,_,_ = openFitsFile(contourFitsFileName)
	fig = plt.figure(figsize=(3.54,3.54),dpi=600)
	ax = fig.add_subplot(projection=fitsWCS)
	image = ax.imshow(fitsData,**kwargs)

	if showContours == True:
		ax.contour(contourData, cmap='Greys_r', alpha=0.5, levels=levels)

	ax.text(155.15,138.05,'1',fontsize=10)
	ax.text(133.11,146.72,'2',fontsize=10)
	ax.text(184.93,120.25,'3',fontsize=10)
	ax.text(182.07,141.55,'4',fontsize=10)
	ax.text(171.74,156.30,'5',fontsize=10)
	#  can use label = f'[fitsHeader["BUNIT"]]' if you want
	#  mom0 label = r'[Jy beam$^{-1}$ km s$^{-1}$]'
	#  mom1 label = r'Velocity offset [km s$^{-1}$]'
	#  mom2 label = r'Velocity dispersion [km s$^{-1}$]'
	plt.colorbar(image,label=r'[Jy beam$^{-1}$ km s$^{-1}$]',fraction=0.035)

	pb.drawBeam(fitsHeader,ax)
	prettifyPlot('Right ascension', 'Declination')
	plt.xlim(100,200)
	plt.ylim(100,170)
	plt.savefig(outFileName,bbox_inches='tight',transparent=False)
	plt.show()

def openFitsFile(fitsFileName):
	hdu = fits.open(fitsFileName)
	fitsData = hdu[0].data[0][0] 
	fitsHeader = hdu[0].header 
	fitsWCS = WCS(fitsHeader,naxis=2)
	return fitsData,fitsHeader,fitsWCS

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=12)
	plt.ylabel(ylabel,fontsize=12)
	#plt.tick_params(axis='both',labelsize=2)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=4, direction='in') #If you want tick marks on the opposite side also, add right=True
	#plt.tick_params(which='minor', length=4)

if __name__ == '__main__':
	plotFits('data/HZ7_Collapsed.fits','data/HZ7_Collapsed.fits','moment0.png',cmap='inferno')
	plotFits('data/HZ7_mom_mask2.weighted_coord.fits','data/HZ7_Collapsed.fits','moment1.png',vmin=-100,vmax=100,cmap='coolwarm')
	plotFits('data/HZ7_mom_mask2.weighted_dispersion_coord.fits','data/HZ7_Collapsed.fits','moment2.png',cmap='inferno')
	plotFits('data/HZ7_mom8_mask2.fits','data/HZ7_Collapsed.fits','moment8.png',cmap='inferno')

	# Jorges Corrections 
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits')
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_mom_mask2.weighted_coord.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits',vmin=-100,vmax=100)
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_mom_mask2.weighted_dispersion_coord.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits')