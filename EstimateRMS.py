''' Function to estimate the background rms of an image'''

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground


def estimate_rms(array_2D):
	sigma_clip = SigmaClip(sigma=3.)
	bkg_estimator = MedianBackground()
	bkg = Background2D(
		array_2D, 
		array_2D.shape, 
		sigma_clip = sigma_clip, 
		bkg_estimator = bkg_estimator,
		)
	return bkg.background_rms_median


if __name__ == '__main__':
	from astropy.io import fits 
	infile = '/home/trystan/Desktop/Work/Hz7_ISM/data/HZ7_Collapsed.fits'
	hdu = fits.open(infile)
	data = hdu[0].data[0][0]
	rms = estimate_rms(data)

