################################################################
#
# Convolving images into one another using the ANA algorithm 
#
################################################################

import numpy as np 
import pylab as plt
from astropy.io import fits 
from astropy.convolution import Gaussian2DKernel 
from astropy.convolution import convolve, convolve_fft
from astropy.wcs import WCS 
from tqdm import tqdm

class HSTImage:
	def __init__(self,filename):
		self.hdu = fits.open(filename)
		self.data = self.hdu[0].data 
		self.header = self.hdu[0].header 
		self.wcs = WCS(self.header,naxis=2)
		self.pixScaleDeg = self.header['CD2_2']
		self.pixScaleArcSec = self.pixScaleDeg*3600
		self.imageSize = self.data.shape

	def writeConvolvedData(self,convolvedDataArray,newFileName):
		self.hdu[0].data = convolvedDataArray
		self.hdu.writeto(newFileName,overwrite=True)

class RadioImage(HSTImage):
	def __init__(self,filename):
		self.hdu = fits.open(filename)
		self.data = self.hdu[0].data[0][0]
		self.header = self.hdu[0].header 
		self.wcs = WCS(self.header,naxis=2)
		self.pixScaleDeg = self.header['CDELT2']
		self.pixScaleArcSec = self.pixScaleDeg*3600
		self.imageSize = self.data.shape 
		self.bmajDeg = self.header['BMAJ']
		self.bminDeg = self.header['BMIN']
		self.bpa = self.header['BPA']
		self.bmajArcSec = self.bmajDeg*3600
		self.bminArcSec = self.bminDeg*3600

def convolveImages(HSTImageObject,radioImageObject):
	FWHMconst = 2.355
	fwhm_HST_arcsec = 0.18 # how do we get this

	sigma_HST_pix_in_ALMA = (fwhm_HST_arcsec/FWHMconst)/(radioImageObject.pixScaleArcSec)
	sigma_HST_pix_in_HST = (fwhm_HST_arcsec/FWHMconst)/(HSTImageObject.pixScaleArcSec)

	sigma_maj_ALMA_pix_in_HST = (radioImageObject.bmajArcSec/FWHMconst)/HSTImageObject.pixScaleArcSec
	sigma_min_ALMA_pix_in_HST = (radioImageObject.bminArcSec/FWHMconst)/HSTImageObject.pixScaleArcSec
	theta_ALMA_degree = radioImageObject.bpa

	sigma_maj_ALMA_pix_in_ALMA = (radioImageObject.bmajArcSec/FWHMconst)/radioImageObject.pixScaleArcSec
	sigma_min_ALMA_pix_in_ALMA = (radioImageObject.bmajArcSec/FWHMconst)/radioImageObject.pixScaleArcSec

	beam_ALMA_in_HST = Gaussian2DKernel(x_stddev = sigma_maj_ALMA_pix_in_HST, y_stddev=sigma_min_ALMA_pix_in_HST, theta = ((theta_ALMA_degree +90) * np.pi)/180,x_size = HSTImageObject.imageSize[0],y_size = HSTImageObject.imageSize[1])
	beam_ALMA_in_ALMA= Gaussian2DKernel(x_stddev = sigma_maj_ALMA_pix_in_ALMA, y_stddev=sigma_min_ALMA_pix_in_ALMA, theta = ((theta_ALMA_degree +90) * np.pi)/180,x_size = radioImageObject.imageSize[0],y_size = radioImageObject.imageSize[1])
	
	psf_HST_in_ALMA = Gaussian2DKernel(x_stddev = sigma_HST_pix_in_ALMA, y_stddev=sigma_HST_pix_in_ALMA, theta = 0, x_size = radioImageObject.imageSize[0], y_size = radioImageObject.imageSize[1])
	psf_HST_in_HST = Gaussian2DKernel(x_stddev = sigma_HST_pix_in_HST, y_stddev=sigma_HST_pix_in_HST, theta = 0, x_size = HSTImageObject.imageSize[0], y_size = HSTImageObject.imageSize[1] )

	convolvedUV = convolve_fft(HSTImageObject.data,beam_ALMA_in_HST,allow_huge=True)
	convolvedALMA = convolve_fft(radioImageObject.data,psf_HST_in_ALMA,allow_huge=True)

	return convolvedUV, convolvedALMA

def generateConvolvedImages(opticalFileName,radioFileName,newOpticalFileName,newRadioFileName):
	Optical = HSTImage(opticalFileName)
	Radio = RadioImage(radioFileName)
	convolvedOptical, convolvedRadio = convolveImages(Optical,Radio)
	Optical.writeConvolvedData(convolvedOptical,newOpticalFileName)
	Radio.writeConvolvedData(convolvedRadio,newRadioFileName)

if __name__ == '__main__':
	radioFile = 'data/HZ7_Collapsed.fits'
	infile1 = 'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits'
	infile2 = 'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected.fits'  
	infile3 = 'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected.fits' 
	infile4 = 'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_total_sci_gaia_corrected.fits'
	infiles = [infile1,infile2,infile3,infile4]

	for infile in tqdm(infiles):
		generateConvolvedImages(infile,radioFile,infile.split('gaiacorrected/')[0]+'convolved/'+infile.split('gaiacorrected/')[1].split('.fits')[0]+'_convolved.fits','TEST_RADIO_CONVOLVED.fits')
