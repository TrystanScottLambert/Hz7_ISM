#############################################################
#
# Script to match resoltions of HST images and ALMA image 
#
############################################################

import numpy as np 
import pylab as plt 
from astropy.io import fits 
from spectral_cube import SpectralCube 
from astropy.wcs import WCS 
from reproject import reproject_interp

ALMAFitsFile = fits.open('data/HZ7_Collapsed.fits') 
ALMASpectralCube = SpectralCube.read(ALMAFitsFile)
ALMAWCS = WCS(ALMASpectralCube.header,naxis=2)


ASCFiles = ['hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved.fits',
				'hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected_convolved.fits',
				'hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected_convolved.fits',
				'hst_13641_07_wfc3_ir_total_sci_gaia_corrected_convolved.fits']

ASCDirectory = 'data/ASCImages/convolved/'


for file in ASCFiles:
	infile = ASCDirectory + file
	ASCFitsFile = fits.open(infile)
	print('IR Resolution (dx,dy) = ', ASCFitsFile[0].header['CD1_1'], ASCFitsFile[0].header['CD2_2'])
	print('HI Resolution (dx,dy) = ', ALMAFitsFile[0].header['CDELT1'], ALMASpectralCube[0].hdu.header['CDELT2'])

	ASCheader = ASCFitsFile[0].header 
	ASCwcs = WCS(ASCheader)
	ASChdu = ASCFitsFile[0].data

	rescaled_ASC_data, _ = reproject_interp(ASCFitsFile,ALMAWCS,shape_out = ASCFitsFile[0].shape)
	rescaled_ASC_imagehdu = fits.PrimaryHDU(data = rescaled_ASC_data, header = ALMASpectralCube[0].hdu.header)

	rescaled_ASC_imagehdu.writeto('data/ASCImages/deresolved/' + file.split('.fits')[0] + '_deresolved.fits')

	image_nan_locs = np.isnan(ASChdu)
	rescaled_ASC_data_nonans = ASChdu
	rescaled_ASC_data_nonans[image_nan_locs] = 0
