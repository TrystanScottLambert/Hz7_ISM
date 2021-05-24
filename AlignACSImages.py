##########################
#
# Alligning ACS images 
#
##########################

import numpy as np 
import pylab as plt
import astroalign as aa 
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval 

infile1 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_acs_wfc_f475w_sci.fits'
infile2 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_acs_wfc_f606w_sci.fits'
infile3 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_acs_wfc_f814w_sci.fits'
infile4 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_acs_wfc_total_sci.fits'

f475 = fits.open(infile1) 
f606 = fits.open(infile2) 
f814 = fits.open(infile3)
ftotal = fits.open(infile4) 

f475Data = f475[1].data
f606Data = f606[1].data 
f814Data = f814[1].data
ftotalData = ftotal[1].data



def locateStarsInImage(imageArray):
	mean,median,std = sigma_clipped_stats(imageArray,sigma=3)
	daofind = DAOStarFinder(fwhm=3,threshold=7*std)
	sources = daofind(imageArray-median)
	return sources

def testImage(imageArray):
	interval = ZScaleInterval()
	limits=interval.get_limits(imageArray)
	sources = locateStarsInImage(imageArray)
	positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
	apertures = CircularAperture(positions, r=4.)
	norm = ImageNormalize(stretch=SqrtStretch())
	plt.imshow(imageArray, cmap='Greys', origin='lower', norm=norm,interpolation='nearest',vmin=limits[0],vmax=limits[1])
	apertures.plot(color='red', lw=1.5, alpha=0.5)
	plt.show()


testImage(ftotalData[1800:3000,1000:3000])
testImage(f475Data[1800:3000,1000:3000])
testImage(f606Data[1800:3000,1000:3000])
testImage(f814Data[1800:3000,1000:3000])






#transf, (s_list,t_list) = aa.find_transform(source,target)



