#########################################################
#
# Script to Measure Flux in Moment1 map and ASC images
#
#########################################################

# To Do:
# Flux/area (or equivalent) plot. 
# FluxIR/FluxUV maybe? 

import numpy as np 
import pylab as plt
from photutils.aperture import CircularAperture 
from photutils.aperture import aperture_photometry
from astropy.io import fits
from astropy.wcs import WCS
from astropy import modeling
import matplotlib
from scipy.optimize import curve_fit
matplotlib.rcParams.update({'font.size': 20})   #Always run this line just to have tick mark sizes looking good. 

#function to work out a circular distribution of circles based on ra,dec and radius
def setCenters(center,radius):
	centers = [center]
	for i in range(6):
		centers.append((2*radius*np.cos((i*np.pi)/3)+center[0],2*radius*np.sin((i*np.pi)/3)+center[1]))
	return centers

def getCircularApetures(center,radius):
	positions = setCenters(center,radius)
	apetures = [CircularAperture(pos,r=radius) for pos in positions]
	return apetures

def getCirclePlottingInformation(positions,radius):
	x = np.array([pos[0] for pos in positions])
	y = np.array([pos[1] for pos in positions])
	r = np.ones(len(x))*radius
	return x,y,r

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)

infile = 'data/HZ7_Collapsed.fits'
hdu = fits.open(infile)
hdu = hdu[0]
data = hdu.data[0][0]
arcsecsInAPixel = hdu.header['CDELT2'] * 3600
squareArcsecsInAPixel = arcsecsInAPixel**2 

pixelRadius = 8
squarePixelsInAperture = pixelRadius*np.pi**2 
arcsecAreaInAperture = squarePixelsInAperture*squareArcsecsInAPixel
center = (149.87692083333334,2.1340558333333335)

def calculateFluxes(data,wcs):
	pixelCenter = wcs.world_to_pixel_values(center[0],center[1])
	pixelCenter = (float(pixelCenter[0]),float(pixelCenter[1]))

	apertures = getCircularApetures(pixelCenter,pixelRadius)
	error = np.sqrt(data)   # have to remove the nans 
	nanLocs = np.isnan(error)
	error[nanLocs]=0

	fluxes = []
	fluxErrors = []
	for aperture in apertures:
		print(aperture_photometry(data,aperture,error=error))
		fluxErrors.append(list(aperture_photometry(data,aperture,error=error)['aperture_sum_err'])[0])
		fluxes.append(list(aperture_photometry(data,aperture,error=error)['aperture_sum'])[0])

	return np.array(fluxes), np.array(fluxErrors),pixelCenter

fluxes,fluxErrors,pixelCenter = calculateFluxes(data,WCS(hdu.header,naxis=2))

######################################################################################

norm = plt.Normalize()
colors = plt.cm.jet(norm(fluxes))
rainbowColorMapValue = plt.get_cmap('rainbow')
colors = rainbowColorMapValue(norm(fluxes))

infile = 'data/HZ7_Collapsed.fits'

hdu = fits.open(infile)
data = hdu[0].data[0][0]
wcs = WCS(hdu[0].header,naxis=2)

# Plotting low-res moment-1 plot
fig = plt.figure()
ax = fig.add_subplot(projection=wcs)
im = ax.imshow(data,cmap='rainbow')

positions = setCenters(pixelCenter,pixelRadius)
x,y,r = getCirclePlottingInformation(positions,pixelRadius)
sc = ax.scatter(x,y,c=fluxes,cmap='rainbow',s=0)
fig.colorbar(im,ax=ax)

for i in range(len(x)):
    currentCircle = plt.Circle((x[i], y[i]), r[i], color = colors[i],fill=False,lw=3)
    ax.add_artist(currentCircle)

prettifyPlot('RA','DEC')
plt.xlim(x[0]-50,x[0]+50)
plt.ylim(y[0]-50,y[0]+50)
plt.show()

plt.subplot(211)
plt.errorbar(np.arange(len(fluxes)),fluxes,yerr=fluxErrors,fmt='or')
#plt.show()


infile = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/deresolved/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved_deresolved.fits'
hdu = fits.open(infile)
hdu = hdu[0]
data = hdu.data
wcs = WCS(hdu.header)

plt.subplot(212)
fluxes,fluxErrors,pixelCenter = calculateFluxes(data,WCS(hdu.header,naxis=2))
plt.errorbar(np.arange(len(fluxes)),fluxes/arcsecAreaInAperture,yerr=fluxErrors/arcsecAreaInAperture,fmt='o')
plt.show()



positions = setCenters(pixelCenter,pixelRadius)
x,y,r = getCirclePlottingInformation(positions,pixelRadius)


fig = plt.figure()
ax = fig.add_subplot(projection=wcs)
im = ax.imshow(data,cmap='rainbow',vmin=0.3,vmax=0.5)

for i in range(len(x)):
    currentCircle = plt.Circle((x[i], y[i]), r[i], color = colors[i],fill=False,lw=3)
    ax.text(x[i],y[i],str(i+1),size=20)
    ax.add_artist(currentCircle)

prettifyPlot('RA','DEC')
plt.xlim(x[0]-50,x[0]+50)
plt.ylim(y[0]-50,y[0]+50)
plt.show()

