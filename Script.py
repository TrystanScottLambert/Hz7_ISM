import casa 
import os
import astropy.io.fits as fits
import numpy as np
import shutil

names = ['HZ7_mom_includepix.integrated','HZ7_mom_includepix.weighted_coord',
			'HZ7_mom_includepix.weighted_dispersion_coord','HZ7_mom_mask1.integrated',
			'HZ7_mom_mask1.weighted_coord','HZ7_mom_mask1.weighted_dispersion_coord',
			'HZ7_mom_mask2.integrated','HZ7_mom_mask2.weighted_coord',
			'HZ7_mom_mask2.weighted_dispersion_coord']

#for name in names:
#	shutil.rmtree(name+'/')

#rm -rf HZ7_mom*
rms_val = 6.1493023326e-05
sigma = 2.5
#infile = 'HZ7_cut.fits.fits'
infile = 'HZ7_local_restfrq.fits'
chan0 = '51~75'
# only inlucding 3.5sigma and above for the moments. 
#casa.immoments(axis='spec', imagename=infile, moments=[0,1,2], outfile='HZ7_mom_includepix',includepix=[sigma*1.45e-4,100],chans=chan0)
casa.immoments(axis='spec', imagename=infile, moments=[0,1,2], outfile='HZ7_mom_includepix',includepix=[sigma*rms_val,100],chans=chan0)

#using a mask created by collapsing the line in casaviewer. Here I manually select that everything above 0.1 mJy is real. 
#casa.immoments(axis='spec', imagename=infile, moments=[0,1,2], outfile='HZ7_mom_includepix_masked',mask='hz7_avg_narrow>=0.0001',stretch=True,chans=chan0)
casa.immoments(axis='spec', imagename='HZ7_mom_includepix', moments=[0,1,2], outfile='HZ7_mom_mask1',mask='HZ7_Collapsed.fits>=0.0001',stretch=True,chans=chan0)

print 'got here 1'
#in this step I run the MaskCube Python script
ff = fits.open(infile)
mask = fits.open('hz7_avg.fits')[0].data[0][0]

for i in range(len(ff[0].data[0])):
	#ff[0].data[0][i][mask<1e-4] = np.nan
	ff[0].data[0][i][mask<sigma*6.15e-5] = np.nan

ff.writeto('HZ7_Trim_masked.fits',overwrite=True)

#New moments, now I am using the Python masked cube. If you take a look to the python cube, I am setting everything below 1e-4 as Nan
#now I need to use the channels where I think the line is. 
casa.immoments(axis='spec', imagename='HZ7_Trim_masked.fits', moments=[0,1,2], outfile='HZ7_mom_mask2',chans=chan0)

full_stats = casa.imstat('hz7_avg',region='noise_region.crtf')
full_rms = full_stats['rms'][0]
print '!!!!!!!!!!!!!!!!!!!'
print full_rms
print '!!!!!!!!!!!!!!!!!!!'


for name in names:
	casa.exportfits(imagename=name, fitsimage=name+'.fits', history=False)
