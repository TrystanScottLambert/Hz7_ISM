"""Testing that the moment_maps can be done without casa.

The tests may seem arbitary but the moments were all made using CASA with as close to
possible values as possible and compared using Ds9, CASAviewer and CARTA and found to
have reproduced the moment maps (with the tinest deviations which would make a computer
declare the tests as failed. Thus, once the moments were manually compared and determined
to be successful, those cubes were taken to be our test cubes.
"""

import unittest
import sys
import os
import numpy as np
from astropy.io import fits
sys.path.append('../')
from HZ7_Analysis import generate_moments


# test_rms = 0.021253365673497466
TEST_CUBE = 'test_cubes/HZ7_Centered.fits'

TEST_MOMENTS_MASKED = {
    'moment0': 'test_cubes/testing_Masked_integrated.fits',
    'moment1': 'test_cubes/testing_Masked_weighted_coord.fits',
    'moment2': 'test_cubes/testing_Masked_weighted_dispersion_coord.fits',

}

TEST_MOMENTS_NOT_MASKED = {
    'moment0': 'test_cubes/testing_integrated.fits',
    'moment1': 'test_cubes/testing_weighted_coord.fits',
    'moment2': 'test_cubes/testing_weighted_dispersion_coord.fits'
}

class test_NonMaskedMomentMaker(unittest.TestCase):
    """Testing that the non masked moments work."""

    def test_generate_standard_moments(self):
        """Only need to test this function seeing as it makes use of the other."""

        test_moment = generate_moments.NonMaskedMomentMaker(TEST_CUBE, 'local_testing')
        test_moment.generate_standard_moment_maps(50,80)
        for value in TEST_MOMENTS_NOT_MASKED.values():
            test_moment_map = fits.open(value)
            local_moment_map = fits.open(value.replace('test_cubes/testing', 'local_testing'))
            np.testing.assert_array_equal(test_moment_map[0].data, local_moment_map[0].data)
            os.remove(value.replace('test_cubes/testing', 'local_testing'))


class test_MaskedMomentMaker(unittest.TestCase):
    """Testing that the masked moments work."""

    def test_generate_standard_moments(self):
        """Only need to test this function seeing as it makes use of the other."""
        
        test_moment = generate_moments.MaskedMomentMaker(TEST_CUBE, 'local_testing')
        test_moment.generate_standard_moment_maps(50,80)
        for value in TEST_MOMENTS_MASKED.values():
            test_moment_map = fits.open(value)
            local_moment_map = fits.open(value.replace('test_cubes/testing', 'local_testing'))
            np.testing.assert_array_equal(test_moment_map[0].data, local_moment_map[0].data)
            os.remove(value.replace('test_cubes/testing', 'local_testing'))


if __name__ == '__main__':
    unittest.main()