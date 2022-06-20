"""Sersic Fitter Routine for Surface profiles and integrated surface flux density."""

from typing import Callable, Tuple
import numpy as np
from scipy.optimize import curve_fit

def fit_function(function: Callable, x_data: np.ndarray,
            y_data: np.ndarray, y_uncertainty: np.ndarray, nstd: float = 3.):
    """Fits a function to x and y data and y_uncertainty."""

    popt, pcov = curve_fit(function, x_data, y_data, sigma = y_uncertainty)
    perr = np.sqrt(np.diag(pcov))
    popt_up = popt + nstd * perr
    popt_down = popt - nstd * perr
    return popt, pcov , perr, popt_up, popt_down

def sersic(radius: float, peak_emission: float,
            effective_radius: float, sersic_index: float, offset: float) -> np.ndarray:
    """Sersic profile with a free n index."""
    const = 2.0 * sersic_index - 0.324
    return peak_emission * np.exp(
        -const * ((radius/effective_radius)**(1.0/sersic_index) - 1.0)) + offset

def sersic_n1(radius: float, peak_emission: float,
              effective_radius: float, offset: float) -> np.ndarray:
    """Sersic prfile with a fixed n=1 index."""
    return sersic(radius, peak_emission, effective_radius, sersic_index = 1, offset = offset)

class SersicFitter:
    """Sersic fitter class."""
    def __init__(self, x_values: np.ndarray, y_values: np.ndarray,
                 y_uncertainties: np.ndarray) -> None:

        self.x_values = x_values
        self.y_values = y_values
        self.y_uncertainties = y_uncertainties

    def fit_sersic(self, sersic_function: Callable) -> Tuple:
        """Function to fit a seris profile."""
        fit = fit_function(sersic_function, self.x_values, self.y_values, self.y_uncertainties)
        return fit

    def fit_sersics(self) -> Tuple:
        """Fitting both a free index and a n=1 sersic fits."""
        fit_sersic_free = self.fit_sersic(sersic)
        fit_sersic_n1 = self.fit_sersic(sersic_n1)
        return fit_sersic_free, fit_sersic_n1
