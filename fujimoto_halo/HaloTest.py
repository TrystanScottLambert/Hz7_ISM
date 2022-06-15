""" Script to test if a galaxy has a [CII] halo as defined by 
Fujimoto et. al., 2019
"""

import numpy as np
import pylab as plt
from astropy.io import fits
import Cosmology as cosmo
import pylab as plt


def calc_equiv_radius(semi_major_axis, semi_minor_axis):
	""" calc_equiv_radius is a function that calculates the radius of a circle
	with the same area as that of an ellipse with the given semi_major_axis
	and semi_minor_axis. These have to be single numbers and a single number
	is returned. 
	"""
	radius = np.sqrt((semi_minor_axis*semi_major_axis) / (4*np.log(2)))
	if semi_major_axis < 0 or semi_minor_axis < 0:
		raise ValueError('Cannot have negative axis')
	
	return radius


def create_circular_mask(shape, center, radius, direction='in'):
	""" code taken from:
	https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array
	shape = (height,width) 
	"""
	Y, X = np.ogrid[:shape[0], :shape[1]]
	dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

	if direction == 'in':
		mask = dist_from_center <= radius
	elif direction == 'out':
		mask = dist_from_center > radius
	else:
		raise AttributeError('direction must be either "in" or "out"')

	return mask

def cutout_annulus(array_2D, shape, center, smaller_radius, larger_radius):
	outer_mask = create_circular_mask(shape, center, larger_radius)
	lcl_data = array_2D.copy()
	lcl_data[~outer_mask] = 0 

	inner_mask = create_circular_mask(shape, center, smaller_radius, direction='out')
	annulus_data = lcl_data.copy()
	annulus_data[~inner_mask] = 0 
	plt.imshow(annulus_data)
	plt.show()
	return annulus_data

class HaloCheck():
	""" class which will determine if the galaxy has a halo accoring to 
	Fujimoto et. al. 2019 """
	def __init__(self, infile, center, redshift, H0 = 70, omega_vacuum = 0.7, omega_matter = 0.3):
		self.center = center
		self.redshift = redshift
		self.H0 = H0
		self.omega_vacuum = omega_vacuum
		self.omega_matter = omega_matter
		self.hdu_list = fits.open(infile)
		self.header = self.hdu_list[0].header 
		self.data = self.hdu_list[0].data[0][0]  # get rid of the stokes cubes
		self.shape = self.data.shape
		self.bpa = self.header['BPA']
		self.bmaj = self.header['BMAJ']
		self.bmin = self.header['BMIN']
		self.deg_per_pix = np.abs(self.header['CDELT1'])
		ellipse_area_pix = np.pi * (self.bmin/self.deg_per_pix) * (self.bmaj/self.deg_per_pix) 
		self.beam_area_pix = ellipse_area_pix / (np.log(2) * 4)

		self.arcsec_per_pix = self.deg_per_pix * 3600
		self.aperture_radius_deg = calc_equiv_radius(self.bmaj,self.bmin)
		self.aperture_radius_arcsec = self.aperture_radius_deg * 3600
		self.aperture_radius_pix = self.aperture_radius_deg / self.deg_per_pix

		self.on_sky_scale = cosmo.Cosmology(
			H0 = self.H0, 
			omega_vacuum = self.omega_vacuum, 
			omega_matter = self.omega_matter,
			).calc_on_sky_scale(self.redshift)

		self.halo_limit_arcsec = 10 / self.on_sky_scale
		self.halo_limit_pix = self.halo_limit_arcsec / self.arcsec_per_pix
		self.halo_data_array = cutout_annulus(self.data, self.shape, self.center, self.aperture_radius_pix, self.halo_limit_pix)
		self.number_of_halo_pixels = len(np.concatenate(self.halo_data_array)) - len(np.where(self.halo_data_array == 0)[0])
		self.halo_sum = np.sum(self.halo_data_array)
		self.do_test()


	def calc_rms(self):
		""" calculate an rms in a annuli which extendes further 
		than 10kpc, but needs to be less than the bounds of 
		the image"""
		rms_annulus = cutout_annulus(self.data, self.shape, self.center, self.halo_limit_pix + 5, self.halo_limit_pix + 35)
		number_of_rms_pixels = len(np.concatenate(rms_annulus)) - len(rms_annulus[np.where(rms_annulus == 0)[0]])
		rms = np.sqrt((1./number_of_rms_pixels) * np.sum(rms_annulus**2))
		return rms

	def do_test(self):
		""" calculating the S/N according to Ana """
		self.rms = self.calc_rms()

		self.flux_density_of_halo = self.halo_sum / self.beam_area_pix
		n_beam = self.number_of_halo_pixels / self.beam_area_pix
		self.halo_noise = self.rms * np.sqrt(n_beam)

		self.signal_to_noise = self.flux_density_of_halo / self.halo_noise

		if self.signal_to_noise > 4:
			self.result = True
			print('[CII] Halo Confirmed')
		else:
			self.result = False
			print('No detection of [CII] Halo')

	def why(self, show_plots = False):
		print(f'Given a user cosmology of H0 = {self.H0} kms/s/Mpc, omega_vacuum = {self.omega_vacuum}, and omega_matter = {self.omega_matter}')
		print(f'A Redshift of {self.redshift} would results in a on sky projection of {self.on_sky_scale} "/kpc.')
		print(f'This means that a radius of 10 kpc would be {self.halo_limit_arcsec} ", or {self.halo_limit_pix} pixels. \n\n')
		print(f'The beam area of the image is {self.beam_area_pix} pixels.')
		print(f'Therefore the radius of a circular aperture with the same area would be {self.aperture_radius_pix} pixels.\n\n')
		print(f'The flux density of the annulus is equal to the sum of all the pixels in the annulus divided by the number of pixels in a beam.')
		print(f'This is calculated as ({self.halo_sum}) / ({self.beam_area_pix}) = {self.flux_density_of_halo} Jy/kms \n\n')
		print(f'The rms is calculated using a larger annulus from {self.halo_limit_pix + 5} pixels to {self.halo_limit_pix +35} pixels.')
		print(f'The rms was calculated to be {self.rms} Jy/kms \n\n')
		print(f'The noise of the image is equal to the rms of the image multiplied by the square root of the number of beam areas in the annulus area.')
		print(f'I.e. noise = rms * sqrt( number of pixels in annulus / number of pixels in beam ),')
		print(f'This results in {self.rms} * sqrt({self.number_of_halo_pixels} / {self.beam_area_pix}) = {self.halo_noise} Jy/kms. \n\n')
		print(f'The signal to noise would then be the flux density divided by the noise.')
		print(f'I.e. {self.flux_density_of_halo} / {self.halo_noise} = {self.signal_to_noise}.')
		if self.result == True:
			print(f'This is greater than 4, and thus the halo test (per Fujimoto et. al., (2021)) passes.')
		else:
			print(f'This is less than 4, and thus the halo test (per Fujimoto et. al., (2021)) fails.')
		pass


if __name__ == '__main__':
	infile  = 'data/HZ7_Collapsed.fits'
	center = (156,139) 
	check = HaloCheck(infile,center,5.25)
