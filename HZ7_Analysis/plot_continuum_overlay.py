##########################################
#                                        #
# Program to make an rgb image in python # 
# and compare it to the continuum image  #
#   									 #
##########################################

import numpy as np 
import pylab as plt 
from astropy.io import fits 
import astropy.units as u 
from astropy.wcs import WCS
from astropy.visualization import simple_norm
import matplotlib as mpl  
from calc_channels import calc_rms
import plotBeam as pb

#making the plots look nice
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 

def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=12)
	plt.ylabel(ylabel,fontsize=12)
	plt.tick_params(axis='both',labelsize=10)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=4, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=2)

infiles = ['data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected.fits',
			'data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected.fits']


imageData = [fits.open(infile)[1].data for infile in infiles]
imageHeader = fits.open(infiles[0])[1].header


continuumCIIDetection = fits.open('data/HZ7_CONTINUUM_250.fits')
continuumCIIData = continuumCIIDetection[0].data[0][0]
rmsValContinuum = calc_rms(continuumCIIData[100:300, 100:300], 20)
continuumLevels = np.arange(2,6)*rmsValContinuum
continuumCIIWCS = WCS(continuumCIIDetection[0].header,naxis=2)

CIIDetection = fits.open('data/HZ7_integrated.fits')
CIIData = CIIDetection[0].data
rmsCII = calc_rms(CIIData, 20)
ciiLevels = np.arange(3,6)*rmsCII
CIIWCS = WCS(CIIDetection[0].header,naxis=2)


# Overlapping Figure
hubbleWCS = WCS(imageHeader)
fig = plt.figure(figsize=(3.54,3.54),dpi=600)
ax = fig.add_subplot(111,projection=hubbleWCS)
ax.imshow(imageData[0],vmin=0.3,vmax=0.55,cmap='gray_r')
ax.contour(continuumCIIData, colors='b', alpha=0.8,lw=100, levels=continuumLevels,transform = ax.get_transform(continuumCIIWCS))
ax.contour(CIIData, colors='r', alpha=0.7,lw=2,levels=ciiLevels,transform = ax.get_transform(CIIWCS))
ax.set_ylim(1068-40,1068-10)
ax.set_xlim(1144-10,1144+20)
SCALE = 6.129 # kpc/" from cosmology calculator
arc_per_pix = imageHeader['CD2_2'] * 3600
error_bar_width = SCALE * arc_per_pix
ax.errorbar(1144-2, 1068-36, xerr = (5 / error_bar_width) /2, color='k', capsize=1.5) #1 kpc scale bar
ax.text( 1144-2, 1070-36, '5 kpc',  horizontalalignment='center', verticalalignment='top')
#pb.drawBeam(CIIDetection[0].header,ax) # Doens't work because pixel scales are different. So this would need to be donw manually
prettifyPlot('RA','Dec')
plt.savefig('plots/Continuum.png',bbox_inches='tight',transparent=False)
plt.show()

# Side by Side figure
fig = plt.figure(figsize=(3.54,3.54),dpi=600)
ax1=fig.add_subplot(121,projection=hubbleWCS)
ax1.imshow(imageData[0],vmin=0.3,vmax=0.55,cmap='gray_r')
ax1.contour(continuumCIIData, cmap='winter', alpha=1,lw=1.5, levels=continuumLevels,transform = ax1.get_transform(continuumCIIWCS))
ax1.set_ylim(1068-40,1068-10)
ax1.set_xlim(1144-10,1144+20)
prettifyPlot('RA','Dec')

ax2 = fig.add_subplot(122,projection=hubbleWCS)
ax2.imshow(imageData[0],vmin=0.3,vmax=0.55,cmap='gray_r')

ax2.contour(CIIData, cmap='autumn', alpha=0.8,lw=2,ls='--' ,levels=ciiLevels,transform = ax2.get_transform(CIIWCS))
ax2.set_ylim(1068-40,1068-10)
ax2.set_xlim(1144-10,1144+20)
prettifyPlot('RA',' ')
plt.savefig('plots/SideBySide.png',bbox_inches='tight',transparent=False)
plt.show()