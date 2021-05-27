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

infile = 'data/HZ7_Combined.fits'
hdu = fits.open(infile)
hduData = hdu[0].data[0]
hduHeader = hdu[0].header
hduWCS = WCS(hduHeader,naxis=2)
matplotlib.rcParams.update({'font.size': 15}) 

'''fig = plt.figure(figsize=(14,11))
grid = ImageGrid(fig, (0.06, 0f.09, 0.40, 0.94), nrows_ncols=(4, 3),axes_pad=0.0, cbar_mode='single', cbar_location='right',share_all=True)
for i in range(12):
	ax = grid[i]
	ax.minorticks_on()
	ax.tick_params(which='major', length=6, width=2, direction='in')
	ax.tick_params(which='minor', length=2, width=2, direction='in')
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')

	im = ax.imshow(hduData[i], cmap='RdYlBu_r', origin='lower')

cb = plt.colorbar(im,cax = grid.cbar_axes[0],pad=0.1)
plt.show()'''

#CubeRMS = np.sqrt(Cube.variance[:, 0, 0]) * 1E3
#MedRMS = np.median(CubeRMS)
#vrange = [-3 * MedRMS, 11 * MedRMS]
#vel = Cube.get_velocity()


fig, axs = plt.subplots(3, 3,figsize=(12,12), sharex=True, sharey=True, subplot_kw={'projection':hduWCS})
fig.subplots_adjust(wspace=0.05, hspace=0.05)
for i in range(9):
	xlabel=r'$\alpha$'
	ylabel=r'$\delta$'
	im = axs.flat[i].imshow(hduData[i]*1e3,cmap='RdYlBu_r')
	raAxes = axs.flat[i].coords[0]
	decAxes = axs.flat[i].coords[1]
	raAxes.set_ticks(exclude_overlapping=True)
	if i == 0 or i == 3:
		axs.flat[i].set_ylabel(ylabel,family='serif')
		raAxes.set_ticklabel_visible(False)
		raAxes.set_axislabel('')
	elif i == 7 or i ==8:
		axs.flat[i].set_xlabel(xlabel,family='serif')
		decAxes.set_ticklabel_visible(False)
		decAxes.set_axislabel('')
	elif i ==6:
		axs.flat[i].set_ylabel(ylabel,family='serif')
		axs.flat[i].set_xlabel(xlabel,family='serif')
	else:
		decAxes.set_ticklabel_visible(False)
		decAxes.set_axislabel('')
		raAxes.set_ticklabel_visible(False)
		raAxes.set_axislabel('')


	axs.flat[i].minorticks_on() #minor ticks work for the multi-plots
	axs.flat[i].tick_params(which='major', length=6, width=2, direction='in')
	axs.flat[i].tick_params(which='minor', length=4)
	axs.flat[i].yaxis.set_ticks_position('both')
	axs.flat[i].xaxis.set_ticks_position('both')
	plt.setp(axs.flat[i].get_xticklabels(), visible=False)
	plt.setp(axs.flat[i].get_yticklabels(), visible=False)

cb_ax = fig.add_axes([0.1, 0.95, 0.8, 0.02])
cbar = fig.colorbar(im, cax=cb_ax,orientation='horizontal',label='mJy/Beam')

plt.show()