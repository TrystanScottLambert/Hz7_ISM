###############################################
#                                             #
# Script to Plot moment maps of Hz7 for paper #
#                                             #
###############################################

import numpy as np 
import pylab as plt
from astropy.io import fits 
from astropy.wcs import WCS
import plotBeam as pb
from calc_channels import calc_rms
import matplotlib as mpl  

#making the plots look nice
mpl.rcParams.update({'font.size': 2})   #Always run this line just to have tick mark sizes looking good. 
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 


CBAR_LABELS = {
	0: r'[Jy beam$^{-1}$ km s$^{-1}$]',
	1: r'Velocity offset [km s$^{-1}$]',
	2: r'Velocity dispersion [km s$^{-1}$]'
}

class MomentMap:
	def __init__(self, infile: str, order: int) -> None:
		self.infile = infile
		self.order = order
	
	def plot_moment(self, outfile: str, contour_name: str = None, **kwargs) -> None:
		fits_data,fits_header, fits_wcs = openFitsFile(self.infile)
		
		fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
		ax = fig.add_subplot(projection = fits_wcs)
		image = ax.imshow(fits_data, **kwargs)

		if contour_name is not None:
			contour_data, _, _  = openFitsFile(contour_name)
			rms_val = calc_rms(contour_data, 20)
			levels = np.arange(3,6)*rms_val
			ax.contour(contour_data, cmap='Greys_r', levels = levels)
		
		plt.colorbar(image, label = CBAR_LABELS[self.order], fraction = 0.035)
		pb.drawBeam(fits_header,ax)
		prettifyPlot('Right ascension', 'Declination')
		plt.xlim(100,200)
		plt.ylim(100,170)
		plt.savefig(outfile, bbox_inches = 'tight', transparent = False)
		plt.clf()
			
def openFitsFile(fitsFileName):
	hdu = fits.open(fitsFileName)
	fitsData = hdu[0].data
	fitsHeader = hdu[0].header 
	fitsWCS = WCS(fitsHeader,naxis=2)
	return fitsData,fitsHeader,fitsWCS

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=12)
	plt.ylabel(ylabel,fontsize=12)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') 
	plt.tick_params(which='major', length=4, direction='in')

if __name__ == '__main__':
	moment_0 = MomentMap('data/HZ7_integrated.fits', order = 0)
	moment_1 = MomentMap('data/HZ7_Masked_weighted_coord.fits', order = 1)
	moment_2 = MomentMap('data/HZ7_Masked_weighted_dispersion_coord.fits', order = 2)

	moment_0.plot_moment('moment0.png', 'data/HZ7_integrated.fits', cmap = 'inferno')
	moment_1.plot_moment('moment1.png', 'data/HZ7_integrated.fits', vmin = -100, vmax = 100, cmap = 'coolwarm')
	moment_2.plot_moment('moment2.png', 'data/HZ7_integrated.fits', cmap = 'inferno')

	'''
	ax.text(155.15,138.05,'1',fontsize=10)
	ax.text(133.11,146.72,'2',fontsize=10)
	ax.text(184.93,120.25,'3',fontsize=10)
	ax.text(182.07,141.55,'4',fontsize=10)
	ax.text(171.74,156.30,'5',fontsize=10)
	'''
	
	#plotFits('data/HZ7_integrated.fits','data/HZ7_integrated.fits','moment0.png',cmap='inferno')
	#plotFits('data/HZ7_Masked_weighted_coord.fits','data/HZ7_integrated.fits','moment1.png',vmin=-100,vmax=100,cmap='coolwarm')
	#plotFits('data/HZ7_Masked_weighted_dispersion_coord.fits','data/HZ7_integrated.fits','moment2.png',cmap='inferno')
	#plotFits('data/HZ7_Masked_fits','data/HZ7_integrated.fits','moment8.png',cmap='inferno')

	# Jorges Corrections 
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits')
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_mom_mask2.weighted_coord.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits',vmin=-100,vmax=100)
	#plotFits('data/Jorge_cut_HZ7/data/HZ7_mom_mask2.weighted_dispersion_coord.fits','data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits')