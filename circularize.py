##########################################################
# 
# Script to do Circular analysis as per Fujimoto 19
#
# Trystan Lambert 
#
###########################################################

import numpy as np 
import pylab as plt
from photutils.aperture import CircularAperture 
from photutils.aperture import aperture_photometry
from astropy.io import fits
from astropy.wcs import WCS
from astropy import modeling
import matplotlib
from scipy.optimize import curve_fit
import os

matplotlib.rcParams.update({'font.size': 20})   #Always run this line just to have tick mark sizes looking good. 

def gaussianFunction(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

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

infile = 'data/HZ7_Centered.fits'
hdu = fits.open(infile)
hdu = hdu[0]

center = ('09:59:30.461','02:08:02.601')
pixelCenter = (155.888,139.095)
arcsecsInAPixel = hdu.header['CDELT2'] * 3600

#arsecondRadius = 0.1
pixelRadius = 8

data = hdu.data
wcsCube = WCS(hdu.header,naxis=3)
restFrequency = hdu.header['RESTFRQ']

apertures = getCircularApetures(pixelCenter,pixelRadius)

fluxes = []
for aperture in apertures:
	fluxes.append([list(aperture_photometry(data[0][val],aperture)['aperture_sum'])[0] for val in range(len(data[0]))])

x = np.arange(len(data[0]))
_,_,freqs = wcsCube.pixel_to_world_values(x,x,x)
vels = 300000*((-freqs + restFrequency)/restFrequency)

for i in range(7):
	f = open(f'nice_circle{i+1}.txt','w')
	for j in range(len(x)):
		f.write(f'{vels[j]} {fluxes[i][j]} \n')
	f.close()

#####################################################################################
infiles = []
for i in range(7):
	infiles.append(f'nice_circle{i+1}.txt')

x_model = np.linspace(-1000,1000,1000)
means = []
meansErrors = []
amplitudes = []
amplitudesErrors = []
fwhms = []
fwhmsErrors = []
for i in range(len(infiles)):
	x,y = np.loadtxt(infiles[i],unpack=True)
	popt,pcov = curve_fit(gaussianFunction,x,y,p0 = [1,0,200])
	perr = np.sqrt(np.diag(pcov))

	means.append(popt[1])
	amplitudes.append(popt[0])
	fwhms.append(popt[2])

	meansErrors.append(perr[1])
	amplitudesErrors.append(perr[0])
	fwhmsErrors.append(perr[2]) 
	
	plt.plot(x_model,gaussianFunction(x_model,*popt),lw=3,label=str(i+1))
	plt.step(x,y,label='data')
	plt.axhline(0,color='k',alpha=0.5,ls=':')
	plt.legend()
	plt.show()

means = np.array(means)
amplitudes = np.array(amplitudes)
fwhms = np.array(fwhms)
######################################################################################

plt.xlabel('Radio Velocity [km/s]')
plt.ylabel('Arbitary units [Jy]')
plt.show()

norm = plt.Normalize()
colors = plt.cm.jet(norm(means))
rainbowColorMapValue = plt.get_cmap('rainbow')
colors = rainbowColorMapValue(norm(means))

# reading in the moment map
infile = 'data/HZ7_Collapsed.fits'
#infile = 'HZ7_mom_includepix.weighted_coord.fits'
#infile =  'HZ7_mom_mask2.weighted_coord.fits'

hdu = fits.open(infile)
data = hdu[0].data[0][0]
wcs = WCS(hdu[0].header,naxis=2)


# Plotting low-res moment-1 plot
fig = plt.figure()
ax = fig.add_subplot(projection=wcs)
im = ax.imshow(data,cmap='rainbow')

positions = setCenters(pixelCenter,pixelRadius)
x,y,r = getCirclePlottingInformation(positions,pixelRadius)
sc = ax.scatter(x,y,c=means,cmap='rainbow',s=0)
fig.colorbar(im,ax=ax,label='km/s')

for i in range(len(x)):
    currentCircle = plt.Circle((x[i], y[i]), r[i], color = colors[i],fill=False,lw=3)
    ax.add_artist(currentCircle)

prettifyPlot('RA','DEC')
plt.show()

plt.errorbar(np.arange(len(means))+1,means,yerr=meansErrors,fmt='o') 
prettifyPlot('Ring','Means')
plt.show()

plt.errorbar(np.arange(len(amplitudes))+1,means,yerr=amplitudesErrors,fmt='o')
prettifyPlot('Ring','Amplitudes')
plt.show()

plt.errorbar(np.arange(len(fwhms))+1,fwhms,yerr=fwhmsErrors,fmt='o')
prettifyPlot('Ring','FWHMS')
plt.show()

os.system('rm nice_circle*.txt') # tear down temp files.