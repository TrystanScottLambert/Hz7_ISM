"""Base class for fits images."""

from abc import ABC, abstractmethod
from typing import List, Tuple
import numpy as np
from astropy.wcs import WCS
from sersic_fitter import sersic, sersic_n1
from pylab import plt
from astropy.cosmology import FlatLambdaCDM

class FitsImage(ABC):
    """Main class for fits images."""
    def __init__(self, infile: str) -> None:
        self.infile = infile

    @property
    @abstractmethod
    def wcs(self) -> WCS:
        """World Coordinate System of the image."""

    @property
    @abstractmethod
    def arcsec_per_pix(self) -> float:
        """Number of arcseconds in a pixel."""

    @abstractmethod
    def generate_radial_profile(self, center: Tuple[int, int],
                radii_pixs: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Creates the radial profile based on an array of radii."""

    @property
    def arcsec_area_per_pix(self) -> float:
        """Area of a pixel in square arcseonds"""
        return self.arcsec_per_pix**2

    @abstractmethod
    def get_surface_profile_params(self, center: Tuple[int, int], radii: np.ndarray):
        """Fits sersic profiles to surface profile data for a given fits image."""
    
    @abstractmethod
    def _plot_settings(self) -> None:
        """Sets all the plot details that we need."""

    def plot_surface_profile(self, x_plotting, y_plotting,
     y_plotting_uncertainties, fit_free_params, fit_n1_params) -> None:
        x_for_fit_plots = np.linspace(x_plotting[0], x_plotting[-1], 1000)
        plt.plot(x_for_fit_plots * self.arcsec_per_pix, sersic(x_for_fit_plots, *fit_free_params[0]), label = f'n = {fit_free_params[0][-2]}')
        plt.plot(x_for_fit_plots * self.arcsec_per_pix, sersic_n1(x_for_fit_plots, *fit_n1_params[0]), label = 'n = 1')
        plt.errorbar(x_plotting * self.arcsec_per_pix, y_plotting, yerr=y_plotting_uncertainties, fmt='rs', ms=5, capsize=3)
        plt.legend()
        plt.xlabel('Radius [arcseconds]')
        self._plot_settings()
        plt.show()

    def do_radial_profile(self, center: Tuple[int, int], radii: np.ndarray) -> None:
        """Plots a radial profile for the image """
        params = self.get_surface_profile_params(center, radii)
        self.plot_surface_profile(*params)
    
    def convert_pixels_to_arcsec(self, array_of_pixels: np.ndarray) -> np.ndarray:
        """Takes and arry of pixels and converts it to arcseconds."""
        return array_of_pixels * self.arcsec_per_pix
    
    def convert_arcsecs_to_kpc(self, array_of_arcsecs: np.ndarray, redshift: float) -> np.ndarray:
        """Converts an array of arcsecs into kpc"""
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(redshift)
        return array_of_arcsecs * kpc_per_arcsec

