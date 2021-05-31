#######################################################
#
# Script to center HZ7_Combined.fits on [CII] line
#
#######################################################

import numpy as np 
from astropy.io import fits 
from astropy import units as u 
from spectral_cube import SpectralCube
from astropy.wcs import WCS

class Cube:
	ciiRestFreq = 1897420620253.1646
	def __init__(self,input_file,output_file,central_channel):
		self.central_channel = central_channel
		self.input_file = input_file
		self.output_file = output_file
		self.hdu = fits.open(self.input_file)

	def centralizeCube(self):
		optimumFrequency = self.findOptimumFrequency()
		self.hdu[0].header['RESTFRQ'] = optimumFrequency
		self.redshift = (-optimumFrequency+self.ciiRestFreq)/optimumFrequency
		self.hdu.writeto(self.output_file,overwrite=True)

	def findOptimumFrequency(self):
		left, right = self.findLeftRight()
		val = 10 
		tolerance = 0.000001
		while np.abs(val) > tolerance:
			midpoint = np.mean([left,right])
			val = self.calculateVelocityFromFrequency(midpoint)
			if val > 0:
				right = midpoint
			elif val < 0:
				left = midpoint
		return midpoint

	def findLeftRight(self):
		startingValue = self.ciiRestFreq
		currentVelValue = self.calculateVelocityFromFrequency(startingValue)
		if currentVelValue < 0:
			while currentVel <0:
				startingValue = startingValue*10
				currentVelValue = self.calculateVelocityFromFrequency(startingValue)
			return ciiRestFreq, startingValue

		elif currentVelValue > 0:
			while currentVelValue > 0:
				startingValue = startingValue/10
				currentVelValue = self.calculateVelocityFromFrequency(startingValue)
			return startingValue, self.ciiRestFreq

		elif currentVelValue == 0:
			print('Already Centered')

	def calculateVelocityFromFrequency(self, frequency):
		cube = SpectralCube.read(self.hdu)
		cube_vel = cube.with_spectral_unit(u.km/u.s, rest_value = frequency* u.Hz, velocity_convention='optical')
		wcs = WCS(cube_vel.header)
		_,_,currentVelRange = wcs.pixel_to_world_values(np.arange(len(cube_vel)),np.arange(len(cube_vel)),np.arange(len(cube_vel)))
		currentVel = currentVelRange[self.central_channel]
		return currentVel


if __name__ == "__main__":
	c = Cube('data/HZ7_Combined.fits','data/HZ7_Centered.fits',62)
	c.centralizeCube()
