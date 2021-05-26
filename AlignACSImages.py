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


astrometry_apikey = 'ywayokgixwuntvsd'
Hz7_RA = 1.498769860000E+02
Hz7_Dec =  2.134113000000E+00

infile1 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits'
infile2 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits'
infile3 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits'
infile4 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_total_sci.fits'
infiles = [infile1,infile2,infile3,infile4]

f475 = fits.open(infile1) 
f606 = fits.open(infile2) 
f814 = fits.open(infile3)
ftotal = fits.open(infile4) 

f475Data = f475[1].data
f606Data = f606[1].data 
f814Data = f814[1].data
ftotalData = ftotal[1].data


#make sure to choose the correct hdu number (for acs this is hdu[1])
def locateStarPositionsInImage(fitsFile):
	imageData = fitsFile.data
	imageWCS = WCS(fitsFile.header)

	xPositions, yPositions = convertSourcesToArrays(locateStarsInImage(imageData))
	raPositions, decPositions = imageWCS.pixel_to_world_values(xPositions,yPositions)

	return raPositions, decPositions


def locateStarsInImage(imageArray):
	mean,median,std = sigma_clipped_stats(imageArray,sigma=3)
	daofind = DAOStarFinder(fwhm=3,threshold=5*std)
	sources = daofind(imageArray-median)

	return sources


def convertSourcesToArrays(sourcesOutput):
	xPositions = np.array(list(sourcesOutput['xcentroid']))
	yPositions = np.array(list(sourcesOutput['ycentroid']))

	return xPositions,yPositions


def plotTestImage(imageArray):
	imageArray = imageArray[1800:4000,1000:4000]
	interval = ZScaleInterval()
	limits = interval.get_limits(imageArray)
	sources = locateStarsInImage(imageArray)
	positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
	apertures = CircularAperture(positions, r=4.)
	norm = ImageNormalize(stretch=SqrtStretch())
	plt.imshow(imageArray, cmap='Greys', origin='lower', norm=norm,interpolation='nearest',vmin=limits[0],vmax=limits[1])
	apertures.plot(color='red', lw=1.5, alpha=0.5)
	plt.show()


def crossMatchPositionalArrays(xpositionArray1,ypositionArray1,xpositionArray2,ypositionArray2,pixelTolerance):
	idxMatches = []
	matchDifferences = []
	for i in range(len(xpositionArray2)):
		differences = np.sqrt((xpositionArray2[i]-xpositionArray1)**2+(ypositionArray2[i]-ypositionArray1)**2)
		closestMatch = np.where(differences<pixelTolerance)[0]

		if len(closestMatch) != 0 :
			closestMatch = closestMatch[0]
			difference = differences[closestMatch]

			if closestMatch in idxMatches:
				duplicatePosition = np.where(np.array(idxMatches)==closestMatch)[0]

				if matchDifferences[duplicatePosition[0]] > difference:
					matchDifferences[duplicatePosition[0]] = closestMatch
			else:
				idxMatches.append(closestMatch)
			matchDifferences.append(difference)

	return idxMatches,matchDifferences


def crossMatchToGaia(gaiaRA,gaiaDec,sourceRA,sourceDec,onSkyLimit):
	gaiaCatalog = SkyCoord(ra=gaiaRA*u.deg,dec=gaiaDec*u.deg,frame='icrs')
	sourceCatalog = SkyCoord(ra=sourceRA*u.deg,dec=sourceDec*u.deg,frame='icrs')

	idxMatches,Distances,_ = gaiaCatalog.match_to_catalog_sky(sourceCatalog)
	separationContraint = Distances < onSkyLimit*u.arcsec

	acceptableGaiaMatches = separationContraint
	acceptableSourceMatches = idxMatches[separationContraint]

	print(f'{len(acceptableSourceMatches)} confirmed matches...')

	return acceptableGaiaMatches, acceptableSourceMatches 


def calculatePositionalDifferences(gaiaRA,gaiaDec,sourceRA,sourceDec,onSkyLimit,plot=False):
	acceptableGaiaMatches, acceptableSourceMatches = crossMatchToGaia(gaiaRA,gaiaDec,sourceRA,sourceDec,onSkyLimit)

	raDifferences = np.abs(gaiaRA[acceptableGaiaMatches] - sourceRA[acceptableSourceMatches])
	decDifferences = np.abs(gaiaDec[acceptableGaiaMatches] - sourceDec[acceptableSourceMatches])

	raOffset = np.mean(raDifferences)
	decOffset = np.mean(decDifferences)
	raOffsetError = np.std(raDifferences)
	decOffsetError = np.std(decDifferences)

	if plot==True:
		plotHistogramDifferences(raDifferences,decDifferences)

	return raOffset, decOffset, raOffsetError, decOffsetError


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




#plotTestImage(ftotalData[1800:3000,1000:3000])
#plotTestImage(f475Data[1800:3000,1000:3000])
#plotTestImage(f606Data[1800:3000,1000:3000])
#plotTestImage(f814Data[1800:3000,1000:3000])


#reading in the gaia data 
raGaia, decGaia = np.loadtxt('data/gaiaSources.txt',unpack=True,delimiter=',',skiprows=1)

for infile in infiles:
	print(infile)
	fitsFileHDU = fits.open(infile)
	fitsFile = fitsFileHDU[1] # for the case of the ACS images
	raSources, decSources = locateStarPositionsInImage(fitsFile)
	raOffset,decOffset,raOffsetError,decOffsetError = calculatePositionalDifferences(raGaia,decGaia,raSources,decSources,10,plot=True)
	correctFitsImage(fitsFile,raOffset,decOffset,infile.split('/')[-1].split('.fit')[0]+str('_gaia_corrected.fits'))

