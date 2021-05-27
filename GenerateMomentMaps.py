import casa 
import os
import astropy.io.fits as fits
import numpy as np

names = ['HZ7_mom_includepix.integrated',
			'HZ7_mom_includepix.weighted_coord',
			'HZ7_mom_includepix.weighted_dispersion_coord',
			'HZ7_mom_mask1.integrated',
			'HZ7_mom_mask1.weighted_coord',
			'HZ7_mom_mask1.weighted_dispersion_coord',
			'HZ7_mom_mask2.integrated',
			'HZ7_mom_mask2.weighted_coord',
			'HZ7_mom_mask2.weighted_dispersion_coord']

###############
#### input ####
###############
rms_val = 2.3e-2  #rms from casa viewer 
sigma = 2.5   
infile = 'data/HZ7_Combined.fits'
emissionChannels = '51~75'   #channels from casa viewer

###########################################################
#### Option 1: Using the includpix to make moment maps ####
###########################################################
casa.immoments(axis='spec', imagename=infile, moments=[0,1,2], outfile='HZ7_mom_includepix',includepix=[sigma*rms_val,100],chans=emissionChannels)

######################################################################################################
#### Option 2: Making a mask from the -1 map and choosing everything above 0.0001 Janskys as real ####
######################################################################################################
#using a mask created by collapsing the line in casaviewer. Here I manually select that everything above 0.1 mJy is real. 
casa.immoments(axis='spec',imagename=infile,moments=[-1],outfile='HZ7_avg',chans=emissionChannels)  # make the HZ7_avg file
casa.immoments(axis='spec', imagename='HZ7_mom_includepix', moments=[0,1,2], outfile='HZ7_mom_mask1',mask='HZ7_avg>=0.0001',stretch=True,chans=emissionChannels)

###########################################################################################################################
#### Option 3: Finding the rms of the moment 0 map and masking all pixels in the data cube then create the moment maps ####
###########################################################################################################################

#create a simple moment map using the selected channels
casa.immoments(axis='spec',imagename=infile,moments=[0],outfile='HZ7_Collapsed',chans=emissionChannels,stretch=True) # Make the colapsed model 
casa.exportfits(imagename='HZ7_Collapsed',fitsimage='HZ7_Collapsed.fits',history=False)

hdu = fits.open(infile)
mask = fits.open('HZ7_Collapsed.fits')[0].data[0][0]  # the moment map we just created. 

for channel in hdu[0].data[0]:
	channel[mask<sigma*rms_val] = np.nan

hdu.writeto('HZ7_Trim_masked.fits',overwrite=True)

#create the moment maps using the masked cube
casa.immoments(axis='spec', imagename='HZ7_Trim_masked.fits', moments=[0,1,2], outfile='HZ7_mom_mask2',chans=emissionChannels)

######################
#### Organization ####
######################

# export all casa images to fits files
for name in names:
	casa.exportfits(imagename=name, fitsimage=name+'.fits', history=False)

os.system('mv HZ7_*.fits data/') # move all the fits files to the data folder 
os.system('rm -rf HZ7_*') #remove all the casa images that remain 