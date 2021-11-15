################################################################
#
# Convolving images into one another using the ANA algorithm 
#
################################################################

# From Roberto
# If they are Gaussian PSFs, you take sigma_best and sigma_worst and you convolve the image with PSF sigma_best with a Gaussian kernel of width sigma = (sigma_worst**2-sigma_best**2)**0.5

import numpy as np 
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
	fwhmHSTArcSec = 0.18 

	sigmaHSTPixInALMA = (fwhmHSTArcSec/FWHMconst)/(radioImageObject.pixScaleArcSec)
	sigmaHSTPixInHST = (fwhmHSTArcSec/FWHMconst)/(HSTImageObject.pixScaleArcSec)

	sigmaMajALMAPixInHST = (radioImageObject.bmajArcSec/FWHMconst)/HSTImageObject.pixScaleArcSec
	sigmaMinALMAPixInHST = (radioImageObject.bminArcSec/FWHMconst)/HSTImageObject.pixScaleArcSec
	thetaALMADegree = radioImageObject.bpa

	sigmaMajALMAPixInALMA = (radioImageObject.bmajArcSec/FWHMconst)/radioImageObject.pixScaleArcSec
	sigmaMinALMAPixInALMA = (radioImageObject.bmajArcSec/FWHMconst)/radioImageObject.pixScaleArcSec

	beamALMAInHST = Gaussian2DKernel(x_stddev = sigmaMajALMAPixInHST, y_stddev=sigmaMinALMAPixInHST, theta = ((thetaALMADegree +90) * np.pi)/180,x_size = HSTImageObject.imageSize[0],y_size = HSTImageObject.imageSize[1])
	beamALMAInALMA= Gaussian2DKernel(x_stddev = sigmaMajALMAPixInALMA, y_stddev=sigmaMinALMAPixInALMA, theta = ((thetaALMADegree +90) * np.pi)/180,x_size = radioImageObject.imageSize[0],y_size = radioImageObject.imageSize[1])
	
	psfHSTInALMA = Gaussian2DKernel(x_stddev = sigmaHSTPixInALMA, y_stddev=sigmaHSTPixInALMA, theta = 0, x_size = radioImageObject.imageSize[0], y_size = radioImageObject.imageSize[1])
	psfHSTInHST = Gaussian2DKernel(x_stddev = sigmaHSTPixInHST, y_stddev=sigmaHSTPixInHST, theta = 0, x_size = HSTImageObject.imageSize[0], y_size = HSTImageObject.imageSize[1] )

	convolvedUV = convolve_fft(HSTImageObject.data,beamALMAInHST,allow_huge=True)
	convolvedALMA = convolve_fft(radioImageObject.data,psfHSTInALMA,allow_huge=True)

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