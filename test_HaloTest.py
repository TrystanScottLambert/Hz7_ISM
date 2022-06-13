""" unit tests for HaloTest.py """

import unittest 
import HaloTest as HT

infile = 'data/HZ7_Collapsed.fits'  # using the HZ7 test case
center = (156,139)
redshift = 5.25

class TestHaloTest(unittest.TestCase):

	def test_calc_equiv_radius(self):
		""" check that the values are correct """
		self.assertEqual(round(HT.calc_equiv_radius(2.4, 5.6), 3), 3.666)

		""" check that zero works correctly """
		self.assertEqual(HT.calc_equiv_radius(0, 4.2), 0)
		self.assertEqual(HT.calc_equiv_radius(4, 0), 0)
		self.assertEqual(HT.calc_equiv_radius(0, 0), 0)

		""" make sure no negative values are allowed """
		with self.assertRaises(ValueError) as context:
			HT.calc_equiv_radius(-2, 3)
		self.assertTrue('Cannot have negative axis' in str(context.exception))

	def test_beam_area(self):
		""" the beam area of HZ7 should alwas be 157.04
		according to casaviewer
		"""
		check = HT.HaloCheck(infile, center, redshift)
		self.assertEqual(round(check.beam_area_pix,2), 157.04)

	#def test_circularize_beam(self):



if __name__ == '__main__':
	unittest.main()