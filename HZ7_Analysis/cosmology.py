''' module to deal with cosmological calculations'''
import numpy as np 
from scipy import integrate


SPEED_OF_LIGHT = 299792.458 # speed of light in km/s
MILIRAD_PER_ARCSEC = 206.264806

class Cosmology():
  ''' main cosmology class used for creating a cosmology'''

  def __init__(self,H0 = 70, omega_vacuum = 0.7, omega_matter = 0.3):
    self.H0 = H0
    self.OL = omega_vacuum
    self.OM = omega_matter
    self.hubble_distance = SPEED_OF_LIGHT/self.H0

    if (H0 < 0) or (omega_vacuum > 1) or (omega_vacuum < 0) or (omega_matter > 1) or (omega_matter < 0):
      raise ValueError('Select Valid Cosmology')

    if omega_vacuum + omega_matter != 1:
      raise ValueError('omega_vacuum and omega_matter must equate to 1')

  def _calc_hubble_parameter(self, z):
    return np.sqrt(self.OM*(1+z)**3 + (1-self.OM-self.OL)**2 + self.OL)

  def _calc_comoving_distance(self, z):
    return self.hubble_distance * integrate.quad(lambda x: 1./self._calc_hubble_parameter(x), 0, z)[0]

  def _calc_comoving_distance_array(self, z_array):
    comoving_distances = [self._calc_comoving_distance(z) for z in z_array] 
    return np.array(comoving_distances)

  def calc_luminosity_distance(self, z): 
    ''' calculates the luminoisty distance ''' 
    if isinstance(z, (list, tuple, np.ndarray)):
      self.luminosity_distance = (1 + np.array(z)) * self._calc_comoving_distance_array(z)
    else:
      self.luminosity_distance = (1 + z) * self._calc_comoving_distance(z) 
      if z < 0:
        raise ValueError('z has to be positive')

    return self.luminosity_distance

  def calc_on_sky_scale(self, z):
    ''' Works out the deg -> kpc value'''
    if isinstance(z, (list, tuple, np.ndarray)):
      self.angular_distance = self._calc_comoving_distance_array(z) / (1 + z)
    else:
      self.angular_distance = self._calc_comoving_distance(z) / (1+z)
      if z < 0:
        raise ValueError('z has to be positive')

    return self.angular_distance / MILIRAD_PER_ARCSEC

if __name__ == '__main__':
  c = Cosmology()
  x = c.calc_on_sky_scale(5.25)
  print(x)