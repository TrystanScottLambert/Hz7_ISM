
import Cosmology as cosmo
import unittest
import numpy as np


class TestCosmology(unittest.TestCase):
	def test_LD_values(self):
		"""
		Test that luminosity densities return
		correct numerical values for defult vals
		"""
		c = cosmo.Cosmology()
		result = c.calc_luminosity_distance(1)
		self.assertEqual(round(result,6),round(6607.6576117749355,6))

		result = c.calc_luminosity_distance(2)
		self.assertEqual(round(result,6),round(15539.58622322812,6))
	
	def test_LD_zerocase(self):
		"""
		Test that LD can handle zero
		"""
		c = cosmo.Cosmology()
		result = c.calc_luminosity_distance(0)
		self.assertEqual(result,0)

	def test_cosmology_limits(self):
		""" 
		Test that the limits on weird cosmology works
		"""
		with self.assertRaises(ValueError) as context:
			cosmo.Cosmology(H0=-20)
		self.assertTrue('Select Valid Cosmology'in str(context.exception))

		with self.assertRaises(ValueError) as context:
			cosmo.Cosmology(omega_matter=0.9,omega_vacuum=0.2)
		self.assertTrue('omega_vacuum and omega_matter must equate to 1' in str(context.exception))


		with self.assertRaises(ValueError) as context:
			cosmo.Cosmology(omega_matter=0.2,omega_vacuum=0.2)
		self.assertTrue('omega_vacuum and omega_matter must equate to 1' in str(context.exception))


		with self.assertRaises(ValueError) as context:
			cosmo.Cosmology(omega_matter=-0.2,omega_vacuum=2.2)
		self.assertTrue('Select Valid Cosmology' in str(context.exception))


		with self.assertRaises(ValueError) as context:
			cosmo.Cosmology(omega_vacuum=-0.2,omega_matter =2.2)
		self.assertTrue('Select Valid Cosmology' in str(context.exception))

	def test_LD_negcase(self):
		"""
		Test to check that calc luminosity distance
		can handle negative input
		"""
		c = cosmo.Cosmology()
		with self.assertRaises(ValueError) as context:
			c.calc_luminosity_distance(-2)
		self.assertTrue('z has to be positive' in str(context.exception))

	def test_array_version_LD(self):
		"""
		Test if the luminosity distance can handle
		list-like objects
		"""
		c = cosmo.Cosmology()
		array = np.linspace(0.2,2.3,20)
		result_1 = c.calc_luminosity_distance(array)
		result_2 = [c.calc_luminosity_distance(value) for value in array]
		self.assertEqual(list(result_1),result_2)

		result_1 = c.calc_luminosity_distance(2.3)
		result_2 = c.calc_luminosity_distance([2.3])
		self.assertEqual(result_1,result_2)

	def test_calc_on_sky(self):
		""" Test if on sky scale gives correct values """
		c = cosmo.Cosmology()
		self.assertEqual(round(c.calc_on_sky_scale(1),1),8.0)
		self.assertEqual(round(c.calc_on_sky_scale(5.25),1),6.1)

		""" Test if on sky scale can handle negatives """
		with self.assertRaises(ValueError) as context:
			c.calc_on_sky_scale(-2)
		self.assertTrue('z has to be positive' in str(context.exception))

		""" Test if on sky scale can handle arrays """
		array = np.linspace(0.2,2.3,20)
		result_1 = c.calc_on_sky_scale(array)
		result_2 = [c.calc_on_sky_scale(value) for value in array]
		self.assertEqual(list(result_1),result_2)

		
if __name__ == '__main__':
	unittest.main()