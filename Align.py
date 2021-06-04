##############################################
#                                            #
# Program to Align HST images to one another #
#                                            #
##############################################

import numpy as np 
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval 
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.io import fits 
from astropy.wcs import WCS
from scipy.stats import norm
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from astroquery.gaia import Gaia
import astropy.units as u 
from astropy.coordinates import SkyCoord

def searchGaiaArchives(ra,dec,height,width):
	coord = SkyCoord(ra = ra*u.deg, dec=dec*u.deg)
	width = u.Quantity(width*u.deg)
	height = u.Quantity(height*u.deg)
	results = Gaia.query_object_async(coordinate=coord, width=width, height=height)
	return np.array(list(results['ra'])), np.array(list(results['dec']))

def getFitsCenter(hdu):
	localWCS = WCS(hdu.header)
	y,x = hdu.shape
	raCenter, decCenter = localWCS.pixel_to_world_values(x/2,y/2)
	return float(raCenter), float(decCenter)

def getFitsDimensions(hdu):
	localWCS = WCS(hdu.header)
	pixHeight, pixWidth = hdu.data.shape
	xPixels = np.array([0,0,pixWidth,pixWidth])
	yPixels = np.array([0,pixHeight,0,pixHeight])
	raPixels,decPixels = localWCS.pixel_to_world_values(xPixels,yPixels)
	return np.abs(raPixels[0] - raPixels[2]), np.abs(decPixels[0]- decPixels[1])

def locateImageHDUValue(hdulist):
	for i in range(len(hdulist)):
		try:
			hdulist[i].header['CRVAL1'] = hdulist[i].header['CRVAL1']
			value = i
		except:
			pass
	return value

def splitTupples(listOfTupples):
	firstElements = [tupple[0] for tupple in listOfTupples]
	secondElements = [tupple[1] for tupple in listOfTupples]
	return firstElements, secondElements

def makeTupples(firstElementsArray, secondElementsArray):
	tupples = np.array([(firstElementsArray[i], secondElementsArray[i]) for i in range(len(firstElementsArray))])
	return tupples

def sumColumnsAndRows(array2D):
	columnSum = np.sum(array2D,axis=0) # axis 0 gives columns  
	rowSum = np.sum(array2D,axis=1)
	return columnSum, rowSum


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gaussFit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt

def cutStar(xPosition,yPosition,imageData):
	xPositionRounded, yPositionRounded = int(round(xPosition)), int(round(yPosition))
	padding = 50
	imageDataPostageStamp = imageData[yPositionRounded-padding:yPositionRounded+padding,xPositionRounded-padding:xPositionRounded+padding]


	localYPosition, localXPosition = getStarPosition(imageDataPostageStamp)
	accurateXPosition, accurateYPosition = xPositionRounded-padding + localYPosition,yPositionRounded-padding + localXPosition

	return imageDataPostageStamp, accurateYPosition, accurateXPosition

def getStarPosition(star):
	columnSum, rowSum = sumColumnsAndRows(star)
	xValues = np.arange(len(columnSum))
	columnPopt = gaussFit(xValues,columnSum)
	rowPopt = gaussFit(xValues,rowSum)
	return columnPopt[2], rowPopt[2]


def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)



class AlignmentDuo:
	def __init__(self, HSTImage1, HSTImage2):
		self.hdulist1 = fits.open(HSTImage1)
		self.hdulist2 = fits.open(HSTImage2)
		self.HSTImage1Data = self.hdulist1[1].data
		self.HSTImage2Data = self.hdulist2[1].data
		self.HSTImage1Header = self.hdulist1[1].header
		self.HSTImage2Header = self.hdulist2[1].header
		self.HSTImage1WCS = WCS(self.HSTImage1Header,naxis=2)
		self.HSTImage2WCS = WCS(self.HSTImage2Header,naxis=2)
		self.plotSideBySide()
		self.generateWorldCoordinates()
		self.calculateXYOffset()

	def plotSideBySide(self):
		self.coords1 = []
		self.coords2 = []
		fig1 = plt.figure()
		interval = ZScaleInterval()

		ax1 = fig1.add_subplot()
		limits = interval.get_limits(self.HSTImage1Data)	
		ax1.imshow(self.HSTImage1Data, cmap='Greys', origin='lower',interpolation='nearest',vmin=limits[0],vmax=limits[1])

		def onclick1(event):
			if event.dblclick:
				print(f'x = {event.xdata}, y = {event.ydata}')
				star,columnPix, rowPix = cutStar(event.xdata,event.ydata,self.HSTImage1Data)
				if np.abs(rowPix-event.xdata) < 15 and np.abs(columnPix-event.ydata) < 15:
					plt.scatter(rowPix, columnPix,facecolors='none',edgecolors='r',s=50)
					fig1.canvas.draw()
					self.coords1.append((rowPix,columnPix))
		cid = fig1.canvas.mpl_connect('button_press_event', onclick1)

		fig2 = plt.figure()
		ax2 = fig2.add_subplot()
		limits = interval.get_limits(self.HSTImage2Data)
		ax2.imshow(self.HSTImage2Data, cmap='Greys', origin='lower',interpolation='nearest',vmin=limits[0],vmax=limits[1])


		def onclick2(event):
			if event.dblclick:
				print(f'x = {event.xdata}, y = {event.ydata}')
				star,columnPix, rowPix = cutStar(event.xdata,event.ydata,self.HSTImage2Data)
				if np.abs(rowPix-event.xdata) < 15 and np.abs(columnPix-event.ydata) < 15:
					plt.scatter(rowPix, columnPix,facecolors='none',edgecolors='b',s=50)
					fig2.canvas.draw()
					self.coords2.append((rowPix,columnPix))
		cid = fig2.canvas.mpl_connect('button_press_event', onclick2)
		plt.show()

	def generateWorldCoordinates(self):
		x1,y1 = splitTupples(self.coords1)
		x2,y2 = splitTupples(self.coords2)
		self.worldCoordinates1 = self.HSTImage1WCS.pixel_to_world_values(x1,y1)
		self.worldCoordinates2 = self.HSTImage2WCS.pixel_to_world_values(x2,y2)

	def calculateXYOffset(self):
		raOffsets = self.worldCoordinates1[0] - self.worldCoordinates2[0]
		decOffsets = self.worldCoordinates1[1] - self.worldCoordinates2[1]
		self.raOffset = np.mean(raOffsets)
		self.decOffset = np.mean(decOffsets)
		self.raSTD = np.std(np.abs(raOffsets))
		self.decSTD = np.std(np.abs(decOffsets))

	def applySimpleCorrection(self,outputFileName):
		self.HSTImage2Header['CRVAL1'] = self.HSTImage2Header['CRVAL1'] + self.raOffset
		self.HSTImage2Header['CRVAL2'] = self.HSTImage2Header['CRVAL2'] + self.decOffset
		self.hdulist2.writeto(outputFileName,overwrite=True)



class GaiaAlignment:
	def __init__(self, inputFile):
		self.hdulist = fits.open(inputFile)
		hduValue = locateImageHDUValue(self.hdulist)
		self.hdu = self.hdulist[hduValue]
		self.hduHeader = self.hdu.header 
		self.hduData = self.hdu.data 
		self.wcs = WCS(self.hduHeader)


		centerRa,centerDec = getFitsCenter(self.hdu)
		width,height = getFitsDimensions(self.hdu)
		self.gaiaRA, self.gaiaDec = searchGaiaArchives(centerRa,centerDec,height,width)
		self.gaiaX, self.gaiaY = self.wcs.world_to_pixel_values(self.gaiaRA,self.gaiaDec)

		self.plotImage()
		self.generateWorldCoordinates()
		self.calculateXYOffset()


	def plotImage(self):
		self.coords = []
		fig = plt.figure()
		interval = ZScaleInterval()

		ax = fig.add_subplot()
		limits = interval.get_limits(self.hduData)	
		ax.imshow(self.hduData, cmap='Greys', origin='lower',interpolation='nearest',vmin=limits[0],vmax=limits[1])
		#ax.set_xlim(0,len(self.hduData[1]))
		#ax.set_ylim(0,len(self.hduData[0]))
		ax.scatter(self.gaiaX,self.gaiaY,marker='*',s=100,facecolors='none',edgecolors='cyan')

		def onclick(event):
			if event.dblclick:
				star,columnPix, rowPix = cutStar(event.xdata,event.ydata,self.hduData)
				if np.abs(rowPix-event.xdata) < 15 and np.abs(columnPix-event.ydata) < 15:
					print(f'x = {event.xdata}, y = {event.ydata}')
					plt.scatter(rowPix, columnPix,facecolors='none',edgecolors='r',s=50)
					fig.canvas.draw()
					self.coords.append((rowPix,columnPix))
		cid = fig.canvas.mpl_connect('button_press_event', onclick)
		plt.show()

	def generateWorldCoordinates(self):
		x ,y = splitTupples(self.coords)
		self.starWorldCoordinates = self.wcs.pixel_to_world_values(x,y)

	def calculateXYOffset(self):
		counter = 0
		raOffsets = []
		decOffsets = []
		print(f'\t Gaia_Star_RA, Gaia_Star_Dec, HST_Star_RA, HST_Star_Dec, RA_Offset, Dec_Offset')
		for i in range(len(self.coords)):
			distances = np.sqrt((self.coords[i][0] - self.gaiaX)**2 + (self.coords[i][1] - self.gaiaY)**2)
			idxMatch = np.where(distances == np.min(distances))[0]

			matchingGaiaRA, matchingGaiaDec = self.gaiaRA[idxMatch], self.gaiaDec[idxMatch]
			currentRAOffset, currentDecOffset = matchingGaiaRA - self.starWorldCoordinates[0][i], matchingGaiaDec - self.starWorldCoordinates[1][i]
			raOffsets.append(currentRAOffset)
			decOffsets.append(currentDecOffset)

			print(f'Gaia Star {i+1}: {matchingGaiaRA}, {matchingGaiaDec}, {self.starWorldCoordinates[0][i]}, {self.starWorldCoordinates[1][i]}, {currentRAOffset}, {currentDecOffset}')

		self.raOffset = np.mean(raOffsets)
		self.decOffset = np.mean(decOffsets)	
		self.raSTD = np.std(np.abs(raOffsets))
		self.decSTD = np.std(np.abs(decOffsets))

		print(f' \n\n TOTAL OFFSET: \n RA: {self.raOffset} +- {self.raSTD} \n Dec: {self.decOffset} +- {self.decSTD}')

	def applySimpleCorrection(self,outputFileName):
		self.hduHeader['CRVAL1'] = self.hduHeader['CRVAL1'] + self.raOffset
		self.hduHeader['CRVAL2'] = self.hduHeader['CRVAL2'] + self.decOffset
		self.hdulist.writeto(outputFileName,overwrite=True)


if __name__ == '__main__':
	infile1 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits'
	infile2 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits'
	infile3 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits'
	infile4 = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_total_sci.fits'
	infiles = [infile1,infile2,infile3,infile4]

	for infile in infiles:
		ga = GaiaAlignment(infile)
		ga.applySimpleCorrection('data/ASCImages/gaiacorrected/' + infile.split('/')[-1].split('.fit')[0] + '_gaia_corrected.fits')



###############################################################################################################################################
gaiaFile = 'data/gaiaSources.txt' # location and name of your gaia sources (https://gea.esac.esa.int/archive/). Save as csv. 
raGaia, decGaia = np.loadtxt(gaiaFile,unpack=True,delimiter=',',skiprows=1)


hdu = fits.open('/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits')
hduData = hdu[1].data
hduHeader = hdu[1].header
wcs = WCS(hduHeader)

xGaia, yGaia = wcs.world_to_pixel_values(raGaia,decGaia)

fig = plt.figure()
ax = fig.add_subplot()
interval = ZScaleInterval()
limits = interval.get_limits(hduData)
ax.imshow(hduData, cmap='Greys', origin='lower',interpolation='nearest',vmin=limits[0],vmax=limits[1])
#ax.scatter(xGaia,yGaia, marker='*', s=100, facecolors='none', edgecolors='cyan')
ax.set_xlim(0,len(hduData[1]))
ax.set_ylim(0,len(hduData[0]))

coords = []
stars = []
def onclick(event):
	if event.dblclick:
		#print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
	    #  (event.button, event.x, event.y, event.xdata, event.ydata))
		star,columnPix, rowPix = cutStar(event.xdata,event.ydata,hduData)
		print(np.abs(rowPix-event.xdata), np.abs(columnPix-event.ydata))
		if np.abs(rowPix-event.xdata) < 15 and np.abs(columnPix-event.ydata) < 15:

			plt.scatter(rowPix, columnPix,color='r',s=10)
			stars.append(star)

			fig.canvas.draw()
			coords.append((event.xdata,event.ydata))

connection_id = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()	



for star in stars:
	columnSum, rowSum = sumColumnsAndRows(star)
	xValues = np.arange(len(columnSum))
	xValuesPlotting = np.linspace(xValues[0],xValues[-1],1000)

	columnPopt = gaussFit(xValues,columnSum)
	rowPopt = gaussFit(xValues,rowSum)


	# start with a square Figure
	fig = plt.figure(figsize=(8, 8))
	gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7), left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)
	ax = fig.add_subplot(gs[1, 0])
	ax_x = fig.add_subplot(gs[0, 0], sharex=ax)
	ax_y = fig.add_subplot(gs[1, 1], sharey=ax)

	ax_x.tick_params(axis="x", labelbottom=False)
	ax_y.tick_params(axis="y", labelleft=False)

	ax.imshow(star,cmap='Greys')
	ax.axhline(rowPopt[2],color='k',ls=':')
	ax.axvline(columnPopt[2],color='k',ls=':')
	penis = getStarPosition(star)
	ax.scatter(penis[0],penis[1])

	ax_x.scatter(xValues,columnSum,color='k')
	ax_x.plot(xValuesPlotting,gauss(xValuesPlotting,*columnPopt),color='r',lw=3)
	ax_x.axvline(columnPopt[2],ls=':',color='k')

	ax_y.scatter(rowSum,xValues,color='k')
	ax_y.plot(gauss(xValuesPlotting,*rowPopt),xValuesPlotting,color='r',lw=3)
	ax_y.axhline(rowPopt[2],ls=':',color='k')


plt.show()
