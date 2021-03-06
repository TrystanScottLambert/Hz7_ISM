""" Module which makes moment 0, 1, and 2 cubes"""

import os
from abc import ABC, abstractmethod
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
from calc_channels import calc_rms


MOMENT_NAMES = {
    0: 'integrated',
    1: 'weighted_coord',
    2: 'weighted_dispersion_coord'
}

class MomentMaker(ABC):
    """General class defning what a moment map should have."""
    def __init__(self, infile: str, outfile: str) -> None:
        self.infile = infile
        self.outfile = outfile

    @abstractmethod
    def generate_moment_map(self, start_channel: int, end_channel: int, order: int) -> None:
        """ Makes a moment map for the give order."""

    def generate_standard_moment_maps(self, start_channel: int, end_channel: int) -> None:
        """Creates the moment 0, 1, and 2 maps."""
        for i in range(3):
            self.generate_moment_map(start_channel, end_channel, order = i)

class NonMaskedMomentMaker(MomentMaker):
    """Manages non masked moment maps."""

    def generate_moment_map(self, start_channel: int, end_channel: int, order: int):
        cube = SpectralCube.read(self.infile)
        cube_kms  = cube.with_spectral_unit(u.km / u.s, velocity_convention = 'optical')
        sub_cube_kms = cube_kms.spectral_slab(
            cube_kms.spectral_axis[start_channel],
            cube_kms.spectral_axis[end_channel])

        if order == 2:
            moment_map = sub_cube_kms.linewidth_sigma()
        else:
            moment_map = sub_cube_kms.moment(order=order)
        moment_map.write(f'{self.outfile}_{MOMENT_NAMES[order]}.fits', overwrite = True)

class MaskedMomentMaker(MomentMaker):
    """Manages masked moment maps."""

    def generate_moment_map(self, start_channel: int, end_channel: int, order: int):
        mid_channel = (start_channel + end_channel) // 2

        cube = SpectralCube.read(self.infile)
        cube_kms  = cube.with_spectral_unit(u.km / u.s, velocity_convention = 'optical')
        rms = calc_rms(cube_kms[mid_channel], 20)
        sub_cube_kms = cube_kms.spectral_slab(
            cube_kms.spectral_axis[start_channel],
            cube_kms.spectral_axis[end_channel])

        masked_cube = sub_cube_kms.with_mask(sub_cube_kms > 3*rms)
        if order == 2:
            moment_map = masked_cube.linewidth_sigma()
        else:
            moment_map = masked_cube.moment(order=order)
        moment_map.write(f'{self.outfile}_Masked_{MOMENT_NAMES[order]}.fits', overwrite = True)


def main():
    """Main function."""
    FILE_NAME = 'data/HZ7_Centered.fits'
    masked_test = MaskedMomentMaker(FILE_NAME, 'testing')
    masked_test.generate_standard_moment_maps(50, 80)

    not_masked_test = NonMaskedMomentMaker(FILE_NAME, 'testing')
    not_masked_test.generate_standard_moment_maps(50, 80)


if __name__ == '__main__':
    main()
