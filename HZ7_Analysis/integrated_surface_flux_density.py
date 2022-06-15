"""Module to calculate the integrated surface flux density.

This is done by taking a moment-0 map and placing an aperture.
One then sums the flux within this aperture and corrects for beam (hence flux density).
The uncertainty will only ever be dependent on the size of the aperture. 
"""

from email.mime import application
import numpy as np
from numpy import ndarray 
import pylab as plt
from astropy.io import fits
import integrated_line_emission.calc_channels as calc_channels
import astropy.units as u 

def flux_density_uncertainty(rms_cube, number_pix_aperture, number_pix_beam):
    return rms_cube * np.sqrt(number_pix_aperture / number_pix_beam)

class Moment0:
    def __init__(self, infile):
        self.hdulist = fits.open(infile)
        self.header = self.hdulist[0].header
        self.data = self.hdulist[0].data[0][0]
        self.unit = self.header['BUNIT']
        if self.unit != u.Jy / u.beam / u.km / u.s:
            print('OHHHHH NOOO')
        self.bmaj = self.header['BMAJ']
        self.bmin = self.header['BMIN']
        self.deg_per_pix = self.header['CDELT2']
        self.square_deg_per_pix = self.deg_per_pix**2
        self.beam_size = (np.pi * (self.bmaj / self.deg_per_pix) * (self.bmin / self.deg_per_pix))/(4 * np.log(2))
        self.rms_thickness = 25
        self.rms = self._calc_rms()

    def read_beam_info(self):
        bmaj = self.header['BMAJ']
        bmin = self.header['BMIN']

    def work_out_flux_density(self, center, radius):
        """Work out flux density for circular aperture with center and radius."""
        cut_out_data = calc_channels.cutout_data(self.data, center, radius)
        #plt.imshow(cut_out_data)
        #plt.savefig(f'delete_when_done_{radius}.png')
        n = len(np.where(cut_out_data != 0)[0])
        print(n)
        flux = np.sum(cut_out_data)
        flux_density = flux / self.beam_size
        uncertainty = flux_density_uncertainty(self.rms, n, self.beam_size)
        return flux_density, uncertainty
    
    def apply_multiple_apertures(self, center, start_radius, end_radius):
        flux_densities = []
        areas = []
        for radius in np.linspace(start_radius, end_radius,100):
            flux_densities.append(self.work_out_flux_density(center, radius)[0])
            areas.append(np.pi * (radius**2))
        return np.array(flux_densities), np.array(areas)

    def _calc_rms(self):
        """Gets RMS of the current moment0."""
        center = calc_channels.locate_center(self.data)
        max_radius = center[0] - 2
        min_radius = max_radius - self.rms_thickness
        rms_data = calc_channels.cutout_annulus(self.data, center, min_radius, max_radius)
        n = len(np.where(rms_data != 0)[0])
        return np.sqrt(np.sum(rms_data**2)/n)



if __name__ == '__main__':
    #infile = 'data/HZ7_Collapsed.fits'
    #center = (154, 138) 

    infile = '../data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits'
    center = (303, 286)

    cube = Moment0(infile)
    x = cube.work_out_flux_density(center, 60)
    fluxes, areas = cube.apply_multiple_apertures(center, 5, 200)
    radiis = np.sqrt(areas / np.pi)
    uncertainties = cube.rms * np.sqrt((areas / cube.beam_size))
    optimum = (fluxes/np.max(fluxes)) - (uncertainties/np.max(uncertainties))
    max_optimum = np.max(optimum)
    ideal_x = np.where(optimum == max_optimum)[0]

    ideal_radii = radiis[ideal_x]

    plt.errorbar(radiis, fluxes, yerr=uncertainties, fmt = 'ro', capsize=3, ms=2,label ='Flux Density Measurements')
    plt.axvline(ideal_radii)
    plt.axhline(0.71, ls = '--', color='k', lw=2, label='Capak Value')
    plt.axhspan(0.71 - 0.07, 0.71 + 0.07, color='k', alpha=0.5) #capak value
    plt.xlabel('Radius', fontsize=20)
    plt.ylabel('Surface Flux Density', fontsize=20)
    plt.legend(fontsize=20)
    plt.show()

    plt.plot(radiis, uncertainties)
    #plt.axvline(ideal_radii)
    plt.show()

    plt.plot(radiis, fluxes/np.max(fluxes) - uncertainties/np.max(uncertainties))
    #plt.axvline(ideal_radii)
    plt.show()