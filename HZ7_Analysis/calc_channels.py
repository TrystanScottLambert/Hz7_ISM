"""Script for calculating the channels used for generating moment maps.

Following the method used in Schouws et. al., (2022) and
Endsely et. al., (2022).
https://arxiv.org/pdf/2202.04080.pdf
https://arxiv.org/pdf/2206.00018.pdf 
"""

from matplotlib.pyplot import plot
import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube
import line_fitting as lf
import pylab as plt

def calc_rms(array_2d, padding):
    """Function to quickly estimate rms of a 2d_array."""
    center = locate_center(array_2d)
    max_radius = center[0] - 2
    min_radius = max_radius - padding
    rms_data = cutout_annulus(array_2d, center, min_radius, max_radius)
    n = len(np.where(rms_data != 0)[0])
    return np.sqrt(np.sum(rms_data**2)/n)


def find_keyword(hdulist, keyword):
    """Find the index in a hdulist with the keyword in it."""
    index = -1
    for i, hdu in enumerate(hdulist):
        try:
            list(hdu.header[keyword])
            index = i
        except KeyError:
            pass

    return index


def calc_equiv_radius(semi_major_axis, semi_minor_axis):
    """Radius of a circle with the same area as that of an ellipse with the.

    given semi_major_axis and semi_minor_axis. These have to be single
    numbers and a single number is returned in native units.
    """
    radius = np.sqrt((semi_minor_axis*semi_major_axis) / (4*np.log(2)))
    if semi_major_axis < 0 or semi_minor_axis < 0:
        raise ValueError('Cannot have negative axis')

    return radius


def split_beam_values_and_units(beam_string):
    """ take the string in pattern '0.001 arcseconds' and
    return the 0.001 as a float value with astropy units
    """
    vals = beam_string.split(' ')
    value = float(vals[0])
    unit = u.Unit(vals[1])
    return value * unit


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


def cutout_data(array_2d, center, larger_radius):
    """ creates an array of the same size but with everything outside of the larger_radius made 0"""
    outer_mask = create_circular_mask(array_2d.shape, center, larger_radius)
    lcl_data = array_2d.copy()
    lcl_data[~outer_mask] = 0
    return lcl_data


def cutout_annulus(array_2d, center, smaller_radius, larger_radius):
    """ creates an array of the same size but only with an annulus of data available between smaller radius and larger radius """
    outer_mask = create_circular_mask(array_2d.shape, center, larger_radius)
    lcl_data = array_2d.copy()
    lcl_data[~outer_mask] = 0

    inner_mask = create_circular_mask(
        array_2d.shape, center, smaller_radius, direction='out')
    annulus_data = lcl_data.copy()
    annulus_data[~inner_mask] = 0
    return annulus_data


def locate_center(array_2d):
    """ returns the central "pixel" of a 2d array"""
    number_of_rows, number_of_columns = array_2d.shape
    return (int(round(number_of_rows/2)), int(round(number_of_columns/2)))


class RadioCube():
    """ main class for calculations"""
    def __init__(self, infile):
        readin_spectral_cube = SpectralCube.read(infile)
        self.spectral_cube = readin_spectral_cube.with_spectral_unit(
            u.km/u.s, velocity_convention='optical')
        self.initial_moment0 = self.spectral_cube.moment0()
        self.bmaj, self.bmin, self.BPA = self.get_beam_info()
        self.equiv_radius = calc_equiv_radius(self.bmaj, self.bmin)
        self.wcs_3d = self.spectral_cube.wcs
        self.wcs_2d = self.initial_moment0.wcs
        self.deg_per_pix = self.spectral_cube.header['CDELT2'] * u.deg
        self.equiv_radius_pix = self.equiv_radius.to(u.deg) / self.deg_per_pix
        self.beam_size = (np.pi * (self.bmaj.value /self.deg_per_pix) *
                          (self.bmin.value / self.deg_per_pix))/(4 * np.log(2))
        self.number_of_deviations = 3

    def get_beam_info(self):
        """ Strip the bmaj, bmin, and Bpa, from the spectral cube moment0 map
        with the units included
        """
        beam_info = self.initial_moment0.header['Beam']
        bmaj = split_beam_values_and_units(
            beam_info.split('BMAJ=')[-1].split(' BMIN')[0])
        bmin = split_beam_values_and_units(
            beam_info.split('BMIN=')[-1].split(' BPA')[0])
        bpa = split_beam_values_and_units(
            beam_info.split('BPA=')[-1])
        return bmaj, bmin, bpa

    def generate_profile_fit(self, mask):
        """ works out the integrated profile and perfoms a gaussian fit"""
        sums = []
        for frame in self.spectral_cube:
            lcl_data = frame.value.copy()
            lcl_data[~mask] = 0
            sums.append(np.sum(lcl_data))
        sums = np.array(sums) 

        fit = lf.GuassianFit(self.spectral_cube.spectral_axis.value, sums,
                              self.spectral_cube.spectral_axis.unit,
                              self.spectral_cube.unit * u.beam)
        return fit

    def calculate_channels(self, center):
        """ finds the channels that we need to collapse the moment map"""
        current_mask = create_circular_mask(self.initial_moment0.shape,
                                            center,
                                            self.equiv_radius_pix.value)
        profile = self.generate_profile_fit(current_mask)

        current_fwfm, current_mean = profile.fwfm[0], profile.mean[0]
        start_vel, end_vel = current_mean - current_fwfm/2, current_mean + current_fwfm/2
        new_start_vel, new_end_vel = np.nan, np.nan
        iterations = 0

        while (start_vel != new_start_vel) and (end_vel != new_end_vel):
            new_start_vel, new_end_vel = start_vel, end_vel
            current_subcube = self.spectral_cube.spectral_slab(
                start_vel * u.km / u.s, end_vel * u.km / u.s)
            current_moment0 = current_subcube.moment0()

            current_rms = calc_rms(current_moment0.value, 20)
            current_mask = current_moment0.value > self.number_of_deviations * current_rms
            profile = self.generate_profile_fit(current_mask)
            plottting_data = current_moment0.value.copy()
            plottting_data[~current_mask] = np.nan

            current_fwfm, current_mean = profile.fwfm[0], profile.mean[0]
            start_vel, end_vel = current_mean - current_fwfm/2, current_mean + current_fwfm/2
            iterations += 1

        start_channel = self.spectral_cube.closest_spectral_channel(
            start_vel * u.km / u.s)
        end_channel = self.spectral_cube.closest_spectral_channel(
            end_vel * u.km / u.s)
        
        self.integrated_emission_profile = profile
        self.moment0_used_for_profile = plottting_data
        return start_channel, end_channel

if __name__ == '__main__':
    infile = '../data/HZ7_Centered.fits'
    center = (156, 140)
    test = RadioCube(infile)

    start_channel, end_channel = test.calculate_channels(center)
    print(start_channel, end_channel)
