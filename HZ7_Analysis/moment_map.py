"""Moment Map Class."""

from typing import Tuple
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from fits_images import FitsImage
from calc_channels import cutout_annulus, cutout_data, calc_rms
from sersic_fitter import SersicFitter
import pylab as plt

ARCSECONDS_IN_DEGREE = 3600.

def flux_density_uncertainty(rms_cube, number_pix_aperture, number_pix_beam):
    """Working out the uncertainty of an aperture size."""
    return rms_cube * np.sqrt(number_pix_aperture / number_pix_beam)


class MomentMap(FitsImage):
    """Main moment0 map class."""
    def __init__(self, infile):
        super().__init__(infile)
        self.hdul = fits.open(infile)
        self.header = self.hdul[0].header
        self.data = self.hdul[0].data

    @property
    def wcs(self) -> WCS:
        """World Coordinate System of the image."""
        return WCS(self.header)

    @property
    def arcsec_per_pix(self) -> float:
        """Number of arcseconds in a pixel."""
        return self.header['CDELT2'] * ARCSECONDS_IN_DEGREE

    @property
    def beam_size(self) -> float:
        """Beam size in pixels. (number of pixels in a beam)."""
        bmin = self.header['BMIN']
        bmaj = self.header['BMAJ']
        deg_per_pix = self.header['CDELT2']
        return (np.pi * (bmaj/deg_per_pix) * (bmin/deg_per_pix))/(4 * np.log(2))

    def generate_radial_profile(self, center: Tuple[int, int], radii_pixs: np.ndarray):
        """Creates the radial profile based on an array of radii."""
        fluxes, uncertainties, plotting_radii, areas = [], [], [], []
        for i in range(len(radii_pixs) - 1):
            if radii_pixs[i] == 0:
                plotting_radii.append(0)
                val = self.calc_aperture_flux_density(center, radii_pixs[i+1])
                pass
            else:
                plotting_radii.append((radii_pixs[i] + radii_pixs[i+1]) / 2)
                val = self.calc_annulus_flux_density(center, radii_pixs[i], radii_pixs[i+1])

            #plotting_radii.append((radii_pixs[i] + radii_pixs[i+1]) / 2)
            fluxes.append(val[0])
            uncertainties.append(val[1])
            areas.append(val[2])
        return np.array(plotting_radii), np.array(fluxes), np.array(uncertainties), np.array(areas)

    def _calc_sum(self, masked_data: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Works out the flux density and uncertainty for a masked region."""

        number_of_pixels = len(np.where(masked_data != 0)[0])
        flux = np.sum(masked_data)
        flux_density = flux / self.beam_size
        uncertainty = flux_density_uncertainty(
            calc_rms(self.data, 20), number_of_pixels, self.beam_size)
        region_area = number_of_pixels * self.arcsec_area_per_pix
        return flux_density, uncertainty, region_area

    def calc_annulus_flux_density(self, center: Tuple[int, int],
                min_radius: float, max_radius: float) -> Tuple[float, float]:
        """Calculates the flux density of a annulus with errors"""

        cut_out_data = cutout_annulus(self.data, center, min_radius, max_radius)
        return self._calc_sum(cut_out_data)

    def get_surface_profile_params(self, center: Tuple[int, int], radii: np.ndarray) -> None:
        """Fits sersic profiles to surface profile data for a given fits image."""
        x_vals, y_vals, y_uncertainties, areas = self.generate_radial_profile(center, radii)
        y_plotting = y_vals / areas
        y_plotting_uncertainties = y_uncertainties / areas
        fit_free_params, fit_n1_params = SersicFitter(x_vals, y_plotting, y_plotting_uncertainties).fit_sersics()
        return x_vals, y_plotting, y_plotting_uncertainties, fit_free_params, fit_n1_params

    def calc_aperture_flux_density(self, center: Tuple[int, int],
                                    radius: float) -> Tuple[float, float]:
        """Calculates the flux density of a annulus with errors"""

        cut_out_data = cutout_data(self.data, center, radius)
        return self._calc_sum(cut_out_data)

    def _plot_settings(self) -> None:
        plt.ylabel('Flux Density / Area [Jy.km/s /arcsecond'+r'$^{2}]$')
        plt.minorticks_on()
        plt.axhline(0)
        plt.tick_params(which='both', width=2,direction='in') 
        plt.tick_params(which='major', length=4, direction='in')


if __name__ == '__main__':
    INFILE = 'data/HZ7_integrated.fits'
    moment0 = MomentMap(INFILE)
    CENTER = (154, 138)
    RADII = np.arange(1, 100, 10)
    moment0.do_radial_profile(CENTER, RADII)
    
