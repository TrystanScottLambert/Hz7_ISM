##############################################
#                                            #
# Script to center ALMA cubes on any channel #
#                                            #
##############################################

import numpy as np 
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from spectral_cube import SpectralCube


class Cube:
	"""Main cube object we use for calculations."""
	CII_REST_FREQ = 1897420620253.1646
	def __init__(self,input_file,output_file,central_channel):
		self.central_channel = central_channel
		self.input_file = input_file
		self.output_file = output_file
		self.hdu = fits.open(self.input_file)

	def centralize_cube(self):
		"""Centralizes the cube."""
		optimum_frequency = self.find_optimum_frequency()
		self.hdu[0].header['RESTFRQ'] = optimum_frequency
		self.redshift = (-optimum_frequency+self.CII_REST_FREQ)/optimum_frequency
		self.hdu.writeto(self.output_file,overwrite=True)

	def find_optimum_frequency(self):
		"""Finding the velocity at 0."""
		left, right = self.find_left_right()
		val = 10 
		tolerance = 0.000001
		while np.abs(val) > tolerance:
			midpoint = np.mean([left,right])
			val = self.calc_vel_from_freq(midpoint)
			if val > 0:
				right = midpoint
			elif val < 0:
				left = midpoint
		return midpoint

	def find_left_right(self):
		"""Finding the upper and lower vals for Newton method."""
		starting_value = self.CII_REST_FREQ
		current_vel_value = self.calc_vel_from_freq(starting_value)
		if current_vel_value < 0:
			while current_vel_value <0:
				starting_value = starting_value*10
				current_vel_value = self.calc_vel_from_freq(starting_value)
			return CII_REST_FREQ, starting_value

		elif current_vel_value > 0:
			while current_vel_value > 0:
				starting_value = starting_value/10
				current_vel_value = self.calc_vel_from_freq(starting_value)
			return starting_value, self.CII_REST_FREQ

		elif current_vel_value == 0:
			print('Already Centered')

	def calc_vel_from_freq(self, frequency):
		""""""
		cube = SpectralCube.read(self.hdu)
		cube_vel = cube.with_spectral_unit(u.km/u.s, rest_value = frequency* u.Hz, velocity_convention='optical')
		wcs = WCS(cube_vel.header)
		_,_,current_velRange = wcs.pixel_to_world_values(np.arange(len(cube_vel)),np.arange(len(cube_vel)),np.arange(len(cube_vel)))
		current_vel = current_velRange[self.central_channel]
		return current_vel


def center_cube_on_channel(infile, outfile, channel):
	return Cube(infile, outfile, channel).centralize_cube()

def main():
	''' main function'''
	#c = Cube('data/HZ7_Combined.fits','data/HZ7_Centered.fits',64)
	#c = Cube('data/Jorge_cut_HZ7/HZ7_COMB_CUBE_15_JvM_cut.fits','data/Jorge_cut_HZ7/HZ7_COMB_CUBE_15_JvM_cut_centered.fits',64)
	#c.centralize_cube()
	infile = '../data/HZ7_Combined.fits'
	center_cube_on_channel(infile, infile.replace('Combined','Centered.fits'), 64)

if __name__ == "__main__":
	main()