##########################
#
# Alligning ACS images 
#
##########################

# get the gaia data for your system from https://gea.esac.esa.int/archive/

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
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import os

##########################
#### Input Parameters ####
##########################
sigmaParameter = 3 # the standard deviation from the mean when working out sigma clipped stats (i.e. excluding extreme values)
fwhmParameter = 3 # the FWHM value to accept a source as a detection (and not a bad pixel)
thresholdParameter = 5  # the number of standard deviations to accept a source as a detection 
gaiaFile = 'data/gaiaSources.txt' # location and name of your gaia sources (https://gea.esac.esa.int/archive/). Save as csv. 
searchRadius = 10 # the arcsecond to determine a match between the gaia catalog and your image sources
##########################

#####################
#### Input Files ####
#####################

#files you want to correct to gaia
infile1 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits'
infile2 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits'
infile3 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits'
infile4 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_total_sci.fits'
infiles = [infile1,infile2,infile3,infile4]

#reading in the gaia data 
raGaia, decGaia = np.loadtxt(gaiaFile,unpack=True,delimiter=',',skiprows=1)

###################
#### Functions ####
###################

#make sure to choose the correct hdu number (for asc this is hdu[1])
def locateStarPositionsInImage(fitsFile):
	imageData = fitsFile.data
	imageWCS = WCS(fitsFile.header)

	xPositions, yPositions = convertSourcesToArrays(locateStarsInImage(imageData))
	raPositions, decPositions = imageWCS.pixel_to_world_values(xPositions,yPositions)

	return raPositions, decPositions


def locateStarsInImage(imageArray):
	mean,median,std = sigma_clipped_stats(imageArray,sigma=sigmaParameter)
	daofind = DAOStarFinder(fwhm=fwhmParameter,threshold=thresholdParameter*std)
	sources = daofind(imageArray-median)

	return sources


def convertSourcesToArrays(sourcesOutput):
	xPositions = np.array(list(sourcesOutput['xcentroid']))
	yPositions = np.array(list(sourcesOutput['ycentroid']))

	return xPositions,yPositions

# A plotting function to help choose the correct values for the input Parameters
def plotTestImage(imageArray):
	interval = ZScaleInterval()
	limits = interval.get_limits(imageArray)
	sources = locateStarsInImage(imageArray)
	positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
	apertures = CircularAperture(positions, r=4.)
	norm = ImageNormalize(stretch=SqrtStretch())
	plt.imshow(imageArray, cmap='Greys', origin='lower', norm=norm,interpolation='nearest',vmin=limits[0],vmax=limits[1])
	apertures.plot(color='red', lw=1.5, alpha=0.5)
	plt.show()

def calculatePositionalDifferences(gaiaRA,gaiaDec,sourceRA,sourceDec,plot=False):
	acceptableGaiaMatches, acceptableSourceMatches = crossMatchToGaia(gaiaRA,gaiaDec,sourceRA,sourceDec)

	raDifferences = np.abs(gaiaRA[acceptableGaiaMatches] - sourceRA[acceptableSourceMatches])
	decDifferences = np.abs(gaiaDec[acceptableGaiaMatches] - sourceDec[acceptableSourceMatches])

	raOffset = np.mean(raDifferences)
	decOffset = np.mean(decDifferences)
	raOffsetError = np.std(raDifferences)
	decOffsetError = np.std(decDifferences)

	if plot==True:
		plotHistogramDifferences(raDifferences,decDifferences)

	return raOffset, decOffset, raOffsetError, decOffsetError

def crossMatchToGaia(gaiaRA,gaiaDec,sourceRA,sourceDec):
	gaiaCatalog = SkyCoord(ra=gaiaRA*u.deg,dec=gaiaDec*u.deg,frame='icrs')
	sourceCatalog = SkyCoord(ra=sourceRA*u.deg,dec=sourceDec*u.deg,frame='icrs')

	idxMatches,Distances,_ = gaiaCatalog.match_to_catalog_sky(sourceCatalog)
	separationContraint = Distances < searchRadius*u.arcsec

	acceptableGaiaMatches = separationContraint
	acceptableSourceMatches = idxMatches[separationContraint]

	print(f'{len(acceptableSourceMatches)} confirmed matches...')

	return acceptableGaiaMatches, acceptableSourceMatches 

def plotHistogramDifferences(raDifferences,decDifferences):
	plt.subplot(211)
	plt.hist(raDifferences,bins=20,lw=3,histtype='step')
	prettifyPlot(r'$\Delta \alpha$','Counts')

	plt.subplot(212)
	plt.hist(decDifferences,bins=20,lw=3,histtype='step')
	prettifyPlot(r'$\Delta \delta$','Counts')
	plt.show()


def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)


def correctFitsImage(fitsFile,raOffset,decOffset,ouputName):
	fitsFile.header['CRVAL1'] = fitsFile.header['CRVAL1'] + raOffset
	fitsFile.header['CRVAL2'] = fitsFile.header['CRVAL2'] + decOffset
	fits.writeto(ouputName,fitsFile.data,fitsFile.header,overwrite=True)

#################
#### Example ####
#################

# Before fully running the script you should make sure that the source finder is finding appropriate sources in your ASC image.
# use the plotTestImage function to see the sources which have been identified in the image BASED ON THE INPUT PARAMATERS

# e.g.
# f105w = fits.open(infile1)
# plotTestImage(f105w[1].data)

# you an align the images in this fasion (I'm just doing all of them at once but you can do one at a time via terminal if you wish)
f = open('data/ASCImages/gaiacorrected/Offsets.txt','w')  # file to write the offsets
f.write('#File RA_offset Dec_offset RA_offset_err Dec_offset_err \n')
for infile in infiles:
	outfile = 'data/ASCImages/gaiacorrected/' + infile.split('/')[-1].split('.fit')[0] + '_gaia_corrected.fits'  # choose your outputfile name
	fitsFileHDU = fits.open(infile)  # opening your file 
	fitsFile = fitsFileHDU[1] # for the case of the ACS image this [1] but use the correct value. 

	raSources, decSources = locateStarPositionsInImage(fitsFile) 
	raOffset,decOffset,raOffsetError,decOffsetError = calculatePositionalDifferences(raGaia,decGaia,raSources,decSources,plot=False) # use plot=True to get a histogram plot of the offsets (though this is usually small numbers)
	correctFitsImage(fitsFile,raOffset,decOffset,outfile) 

	f.write(f'{infile} {raOffset} {decOffset} {raOffsetError} {decOffsetError} \n')

	# optional # 
	# you can (and probably should) do a quick scatter plot to make sure that the number of matches and the types of matches are looking reasonable
	# plt.scatter(raGaia, decGaia, s = 100)
	# plt.scatter(raSource, decSource, s = 50)
	# plt.show()
f.close()

# Be sure to double check alignment in ds9