""" Module which makes moment 0, 1, and 2 cubes"""

from spectral_cube import SpectralCube
import astropy.units as u

MOMENT_NAMES = {
    0: 'integrated',
    1: 'weighted_coord',
    2: 'weighted_dispersion_coord'
}

def generate_moment_map(fitsfile: str, outfile: str,
 starting_channel: int, ending_channel :int, order: int) -> None:
    """ makes moments maps for given order"""

    cube = SpectralCube.read(fitsfile)
    cube_kms  = cube.with_spectral_unit(u.km / u.s, velocity_convention = 'optical')
    moment_map = cube_kms.moment(order=order)
    moment_map.write(f'{outfile}_{MOMENT_NAMES[order]}.fits')

def generate_stanadrd_moments(fitsfile: str, outfile: str,
 starting_channel: int, ending_channel: int) -> None:
    """ This will generate the standard moment 0, moment 1, and moment 2 maps."""

    for i in range(3):
        generate_moment_map(fitsfile, outfile,
         starting_channel, ending_channel, order = i)

    