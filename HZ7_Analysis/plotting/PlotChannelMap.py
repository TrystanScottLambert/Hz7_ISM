########################################################
#                                                      #
# Script to Generate Channel Maps of HZ7_Combined.fits # 
#                                                      #
########################################################

import numpy as np 
import pylab as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits 
from astropy.wcs import WCS
import matplotlib
from spectral_cube import SpectralCube
import astropy.units as u 
import plotBeam as pb
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 5
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=4) 
mpl.rc('ytick', labelsize=4) 

infiles= ['data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected.fits']

imageData = [fits.open(infile)[1].data for infile in infiles]
imageHeader = fits.open(infiles[0])[1].header
imageWCS = WCS(imageHeader)
maskedCIIDetection = fits.open('data/HZ7_mom_mask2.integrated.fits')
rms_val = 2.3e-4
levels = np.arange(3,6)*rms_val
negativeLevels = np.arange(-5,-2)*rms_val
maskedCIIData = maskedCIIDetection[0].data
ciiWCS = WCS(maskedCIIDetection[0].header,naxis=2)
hubbleWCS = WCS(imageHeader)

def getSpectralValues(cube):
	wcs = WCS(cube.header,naxis=3)
	print(wcs)
	_, _, spectralValues = wcs.pixel_to_world_values(np.ones(len(cube)), np.ones(len(cube)),np.arange(len(cube)))
	return  spectralValues

def getChannels(box):
	numberFrames = box[0]**2
	if numberFrames%2 != 0:
		startingChannel = centralChannel - np.floor(numberFrames/2)
		endingChannel = centralChannel + np.ceil(numberFrames/2 + 1)
	else:
		startingChannel = centralChannel - numberFrames/2
		endingChannel = centralChannel + numberFrames/2
	return int(startingChannel), int(endingChannel)

def getPlottingChannelInfo(box):
	startingChannel, endingChannel = getChannels(box)
	leftPlots = np.arange(0,numberFrames-box[0],box[0])
	bottomPlots = np.arange(numberFrames-box[0]+2,numberFrames+1) - 1
	cornerPlot = numberFrames-box[0]
	return startingChannel, endingChannel, leftPlots, bottomPlots, cornerPlot

######################
#### Input Params ####
######################
box = (3,3)
centralChannel = 62

numberFrames = box[0]**2
startingChannel, endingChannel, leftPlots, bottomPlots, cornerPlot = getPlottingChannelInfo(box)

###################
### Input File ####
###################
infile = 'data/HZ7_Centered.fits'
hdu = fits.open(infile)
hduData = hdu[0].data[0]
hduData = hduData[startingChannel:endingChannel]
hduHeader = hdu[0].header

hduWCS = WCS(hduHeader,naxis=2)

hduBeamData = hdu[1].data   # Beam Table
beamData = hduBeamData[startingChannel:endingChannel]

restFrequency = hduHeader['RESTFRQ']
CDELT = hduHeader['CDELT1']
cube = SpectralCube.read(hdu)
cubeVel = cube.with_spectral_unit(u.km/u.s, rest_value = restFrequency*u.Hz, velocity_convention='optical')
velocities = getSpectralValues(cubeVel)
velocitiesInChannels = velocities[startingChannel:endingChannel]

fig, axs = plt.subplots(box[0], box[1],figsize=(3.54,3.54),dpi=600, sharex=True, sharey=True, subplot_kw={'projection':hduWCS})  #hubbleWCS for RGB overlay
fig.subplots_adjust(wspace=0.05, hspace=0.05)
for i in range(numberFrames):
	xlabel=r'$\alpha$'
	ylabel=r'$\delta$'
	ax = axs.flat[i]
	#im = ax.imshow(hduData[i]*1e3,vmin = -1.5, vmax= 1.5, cmap='RdYlBu_r')
	im = ax.imshow(imageData[0],vmin=0.35,vmax=0.6,transform=ax.get_transform(imageWCS),cmap='gray_r')
	ax.set_ylim(150-50,150+50)
	ax.set_xlim(142-50,142+50)
	ax.errorbar(142+30, 150-35, xerr=5*6.287/2, color='k', capsize=1.5)
	ax.text( 142+30, 146-35, '5 kpc',  horizontalalignment='center', verticalalignment='top')

	lon = ax.coords[0]
	lat = ax.coords[1]
	lon.set_major_formatter('hh:mm')
	lat.set_major_formatter('hh:mm')	
	raAxes = ax.coords[0]
	decAxes = ax.coords[1]
	raAxes.set_ticks(exclude_overlapping=True)
	if i in leftPlots:
		ax.set_ylabel(ylabel,family='serif')
		raAxes.set_ticklabel_visible(False)
		raAxes.set_axislabel('')
	elif i in bottomPlots:
		ax.set_xlabel(xlabel,family='serif')
		decAxes.set_ticklabel_visible(False)
		decAxes.set_axislabel('')
	elif i == cornerPlot:
		ax.set_ylabel(ylabel,family='serif')
		ax.set_xlabel(xlabel,family='serif')
	else:
		decAxes.set_ticklabel_visible(False)
		decAxes.set_axislabel('')
		raAxes.set_ticklabel_visible(False)
		raAxes.set_axislabel('')

	pb.drawBeamManually(beamData[i][0],beamData[i][1],beamData[i][2],CDELT,ax)
	#ax.imshow(np.dstack([imageData[0],imageData[1],imageData[2]]),vmin=0.3,vmax=0.6)
	ax.contour(hduData[i]*1e3, cmap='rainbow',linewidths=1, levels=levels*1e3,transform = ax.get_transform(ciiWCS),lw=1)
	ax.contour(hduData[i]*1e3, cmap='gray_r', linestyles=':',linewidths=0.5, levels=negativeLevels*1e3,transform = ax.get_transform(ciiWCS))
	ax.text(0.5, 0.85, str(round(1e-3*velocitiesInChannels[i]))+ r' km s$^{-1}$' , transform=ax.transAxes, fontsize=5,color='black', bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 2},ha='right')
	ax.minorticks_on() #minor ticks work for the multi-plots
	ax.tick_params(which='major', length=2, width=1, direction='in')
	ax.tick_params(which='minor', length=4)
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax.get_yticklabels(), visible=False)

#Only put on colorbar when not using the HST image
#cb_ax = fig.add_axes([0.1, 0.95, 0.8, 0.02])
#cbar = fig.colorbar(im, cax=cb_ax,orientation='horizontal',label='mJy/Beam')

plt.savefig('ChannelMaps.png',bbox_inches='tight',transparent=False)
plt.show()