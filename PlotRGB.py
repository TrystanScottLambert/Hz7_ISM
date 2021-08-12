##########################################
#                                        #
# Program to make an rgb image in python # 
#                                        #
##########################################

import numpy as np 
import pylab as plt 
from astropy.io import fits 
import astropy.units as u 
from astropy.wcs import WCS

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)

infiles= ['data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected.fits']

imageData = [fits.open(infile)[0].data for infile in infiles]
imageHeader = fits.open(infiles[0])[0].header

maskedCIIDetection = fits.open('data/HZ7_mom_mask2.integrated.fits')
rms_val = 2.3e-2
levels = np.arange(4,8)*rms_val
maskedCIIData = maskedCIIDetection[0].data
ciiWCS = WCS(maskedCIIDetection[0].header,naxis=2)


hubbleWCS = WCS(imageHeader)
fig = plt.figure()
ax = fig.add_subplot(111,projection=hubbleWCS)
plt.imshow(np.dstack([imageData[0],imageData[1],imageData[2]]),vmin=0.35,vmax=0.45)
ax.contour(maskedCIIData[0][0], cmap='Greys_r', alpha=0.5, levels=levels,transform = ax.get_transform(ciiWCS))
prettifyPlot('RA','Dec')
plt.show()