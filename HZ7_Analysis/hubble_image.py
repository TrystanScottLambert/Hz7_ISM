"""Module for handling Hubble Image data."""

from typing import Tuple
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from calc_channels import cutout_annulus, cutout_data
from fits_images import FitsImage
from sersic_fitter import SersicFitter
import pylab as plt

ARCSECONDS_IN_DEGREE = 3600.

def _calc_local_rms(hubble_image: np.ndarray, center: Tuple[int, int]) -> float:
    """Calculates the rms using a 20 pixel anulus, starting at 20pix radius around the source."""
    data = cutout_annulus(hubble_image, center, 20, 40)
    number_of_pixels = len(np.where(data != 0)[0])
    return np.sqrt(np.sum(data**2)/number_of_pixels)


class HubbleImage(FitsImage):
    """Main class for reading in a hubble image."""
    def __init__(self, infile) -> None:
        super().__init__(infile)
        self.hdul = fits.open(infile)
        self.header = self.hdul[1].header
        self.data = self.hdul[1].data

    @property
    def wcs(self) -> WCS:
        """The world coordinate system."""
        return WCS(self.hdul[1].header)

    @property
    def arcsec_per_pix(self) -> float:
        """Number of arcsecond in a single pixel."""
        return self.hdul[1].header['CD2_2'] * ARCSECONDS_IN_DEGREE

    @property
    def photflam(self) -> float:
        """Inverse sensitivity. Should be in units of ergs/cm2/Ang/electron"""
        return self.hdul[0].header['PHOTFLAM']

    @property
    def photplam(self) -> float:
        """Pivot wavelength in Angstroms"""
        return self.hdul[0].header['PHOTPLAM']

    @property
    def photfnu(self) -> float:
        """Flux conversion in units of Jy * sec / electron"""
        return self.hdul[0].header['PHOTFNU']

    @property
    def photzpt(self) -> float:
        """ST magnitude zero point."""
        return self.hdul[0].header['PHOTZPT']

    def calculate_stmags(self, fluxes: np.ndarray) -> np.ndarray:
        """Works out the ST magnitudes for given flux measurements.

        See https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry
        for more details, specifically section 9.1.1
        """

        return -2.5 * np.log10(fluxes * self.photflam) + self.photzpt

    def calculate_abmags(self, fluxes: np.ndarray) -> np.ndarray:
        """Converts the STmags into AB mags."""
        return self.calculate_stmags(fluxes) - 5 * np.log10(self.photplam) + 18.692

    def _calc_sum(self, masked_data: np.ndarray, center: Tuple[int, int]):
        """Works out the flux density and uncertainty for a masked region."""

        number_of_pixels = len(np.where(masked_data != 0)[0])
        flux_sum = np.sum(masked_data)
        region_area = number_of_pixels * self.arcsec_area_per_pix
        rms = _calc_local_rms(self.data, center)
        flux_uncertainty = np.sqrt(number_of_pixels*rms**2 + flux_sum)
        return flux_sum, flux_uncertainty, region_area

    def calc_sum_anulus(self, center: Tuple[int, int],
                start_radius: float, end_radius: float) -> Tuple[float, float]:
        """Calculate the sum with uncertainty using an anulus."""

        cut_out_data = cutout_annulus(self.data, center, start_radius, end_radius)
        return self._calc_sum(cut_out_data, center)

    def calc_sum_aperture(self, center: Tuple[int, int], radius: float):
        """Calculate the sum with uncertainty using a circular aperture."""
        cut_out_data = cutout_data(self.data, center, radius)
        return self._calc_sum(cut_out_data, center)

    def generate_radial_profile(self, center: Tuple[int, int], radii: np.ndarray):
        """Creates a surface profile"""
        fluxes, uncertainties, plotting_radii, areas = [], [], [], []
        for i in range(len(radii) - 1):
            if radii[i] == 0:
                plotting_radii.append(0)
                val = self.calc_sum_aperture(center, radii[i+1])
            else:
                plotting_radii.append((radii[i] + radii[i+1]) / 2)
                val = self.calc_sum_anulus(center, radii[i], radii[i+1])

            fluxes.append(val[0])
            uncertainties.append(val[1])
            areas.append(val[2])
        return np.array(plotting_radii), np.array(fluxes), np.array(uncertainties), np.array(areas)
    
    def get_surface_profile_params(self, center: Tuple[int, int], radii: np.ndarray):
        """Fits sersic profiles to surface profile data for a given fits image."""
        x_vals, y_vals, y_uncertainties, areas = self.generate_radial_profile(center, radii)
        y_plotting = y_vals / areas

        surface_brightness = self.calculate_stmags(y_plotting)
        surface_brightness_uncertainty = 2.5 * (y_uncertainties / (y_vals  * np.log(10)))
        fit_free_params, fit_n1_params = SersicFitter(x_vals, surface_brightness, surface_brightness_uncertainty).fit_sersics()

        return x_vals, surface_brightness, surface_brightness_uncertainty, fit_free_params, fit_n1_params
    
    def _plot_settings(self):
        plt.gca().invert_yaxis()
        plt.ylabel('Surface brightness [mag /arcsecond'+r'$^{2}]$')
        plt.minorticks_on()
        plt.tick_params(which='both', width=2,direction='in') 
        plt.tick_params(which='major', length=4, direction='in')


if __name__ == '__main__':
    INFILE = 'data/ASCImages/convolved/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved.fits'
    F105W = HubbleImage(INFILE)
    CENTER = (1150,1046)
    RADII = np.arange(0, 22, 1)
    params = F105W.get_surface_profile_params(CENTER, RADII)  
    F105W.plot_surface_profile(*params)
