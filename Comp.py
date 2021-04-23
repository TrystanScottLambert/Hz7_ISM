# Script comparing HZ7 u-v emission to [CII] emission
#
#

import numpy as np 
import pylab as plt
from astropy.io import fits
import astropy.units as u 
from astropy.wcs import WCS
from astropy.visualization import simple_norm

infile = 'test_mom.fits'
moment0 = fits.open(infile)
data_CII_all = moment0[0].data[0][0]
Header_CII_all = moment0[0].header
wcs_all = WCS(Header_CII_all,naxis=2)


infile = 'HZ7_mom_mask2.integrated.fits'
moment0 = fits.open(infile)
data_CII = moment0[0].data[0][0]
Header_CII = moment0[0].header
levels = np.arange(0.05,0.25,0.05)

#infile = 'hst_13641_07_wfc3_ir_f160w_sci.fits'
#infile = 'hst_13641_07_wfc3_ir_f105w_sci.fits'
#hdu = fits.open(infile)
#data_uv = hdu[1].data
#header_uv = hdu[1].header 

infile = 'gaia_corrected_rescaled_acs.fits'
infile = 'Hubble_Image.fits'
hdu = fits.open(infile)
data_uv = hdu[0].data
header_uv = hdu[0].header

wcs = WCS(header_uv,naxis=2)
wcs1 = WCS(header_uv,naxis=2)
fig = plt.figure()
ax = fig.add_subplot(projection=wcs)
norm = simple_norm(data_uv,'log')
#ax.imshow(data_uv,norm=norm,vmax = 0.7,vmin = 0.3)
ax.imshow(data_uv)
ax.contour(data_CII, transform = ax.get_transform(wcs_all), levels = levels,cmap='Greys',alpha=0.5)
plt.show()
