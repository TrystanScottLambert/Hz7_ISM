""" Module which makes moment 0, 1, and 2 cubes"""

from spectral_cube import SpectralCube
import astropy.units as u
from abc import ABC, abstractmethod
from astropy.io import fits
from calc_channels import calc_rms
import os

MOMENT_NAMES = {
    0: 'integrated',
    1: 'weighted_coord',
    2: 'weighted_dispersion_coord'
}

class MomentMaker(ABC):
    def __init__(self, infile: str, outfile: str) -> None:
        self.infile = infile
        self.outfile = outfile
    
    @abstractmethod
    def generate_moment_map(self, start_channel: int, end_channel: int, order: int):
        """ Makes a moment map for the give order."""
        pass
    
    def generate_standard_moment_maps(self, start_channel: int, end_channel: int, order: int) -> None:
        for i in range(3):
            self.generate_moment_map(start_channel, end_channel, order = i)
    
class NonMaskedMomentMaker(MomentMaker):
    def generate_moment_map(self, start_channel: int, end_channel: int, order: int):
        cube = SpectralCube.read(self.infile)
        cube_kms  = cube.with_spectral_unit(u.km / u.s, velocity_convention = 'optical')
        sub_cube_kms = cube_kms.spectral_slab(
            cube_kms.spectral_axis[start_channel],
            cube_kms.spectral_axis[end_channel])

        moment_map = sub_cube_kms.moment(order=order)
        moment_map.write(f'{self.outfile}_{MOMENT_NAMES[order]}.fits', overwrite = True)
    
class MaskedMomentMaker(MomentMaker):
    def generate_moment_map(self, start_channel: int, end_channel: int, order: int):
        #create a test moment-0 map to work out the rms
        temp_maker = NonMaskedMomentMaker(self.infile, 'temp')
        temp_maker.generate_moment_map(start_channel, end_channel, order = 0)

        temp_string = f'temp_{MOMENT_NAMES[0]}.fits'
        moment0 = fits.open(temp_string)
        rms = calc_rms(moment0[0].data, 20)
        mask = moment0[0].data > 3 * rms

        cube = SpectralCube.read(self.infile)
        cube_kms  = cube.with_spectral_unit(u.km / u.s, velocity_convention = 'optical')
        sub_cube_kms = cube_kms.spectral_slab(
            cube_kms.spectral_axis[start_channel],
            cube_kms.spectral_axis[end_channel])

        masked_cube = sub_cube_kms.with_mask(mask)
        moment_map = masked_cube.moment(order=order)
        moment_map.write(f'{self.outfile}_Masked_{MOMENT_NAMES[order]}.fits', overwrite = True)

        delete_file(temp_string)


def delete_file(file_name: str) -> None:
    if os.path.exists(file_name):
        os.remove(file_name)

if __name__ == '__main__':
    infile = 'data/HZ7_Centered.fits'
    test = NonMaskedMomentMaker(infile, 'delete')
    test.generate_moment_map(50, 80, order=0)
    another_test = MaskedMomentMaker(infile, 'delete')
    another_test.generate_moment_map(50, 80, order=1)
