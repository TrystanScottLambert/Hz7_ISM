"""Sersic Fitter Routine for Surface profiles and integrated surface flux density."""

from typing import Callable, Tuple
import numpy as np
from scipy.optimize import curve_fit

def fit_function(function: Callable, x_data: np.ndarray,
            y_data: np.ndarray, y_uncertainty: np.ndarray, nstd: float = 3., initial_guess: np.ndarray = None):
    """Fits a function to x and y data and y_uncertainty."""

    popt, pcov = curve_fit(function, x_data, y_data, sigma = y_uncertainty, p0 = initial_guess)
    perr = np.sqrt(np.diag(pcov))
    popt_up = popt + nstd * perr
    popt_down = popt - nstd * perr
    return popt, pcov , perr, popt_up, popt_down

def bn(n: float) -> float:
    """Calculates the bn term in a serisc profile for n > 0.36"""
    bn = 2*n -1./3 + 4./(405*n) + 46./(25515*n**2) + 131./(1148175*n**3) - 2194697./(30690717750*n**4)
    return bn

def sersic(radius: float, peak_emission: float,
            effective_radius: float, sersic_index: float, offset: float) -> np.ndarray:
    """Sersic profile with a free n index."""
    return peak_emission * np.exp(
        -(2*sersic_index -1./3) * ((radius/effective_radius)**(1.0/sersic_index) - 1.0)) + offset

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
        guess_amp = 0.5 * np.max(self.y_values)
        guess_reff = 0.5 * np.max(self.x_values)
        guess_sersic_index = 0.75
        guess_offset = np.min(self.y_values)
        if sersic_function == sersic_n1:
            initial_guess = np.array([guess_amp, guess_reff, guess_offset])
        else:
            initial_guess = np.array([guess_amp, guess_reff, guess_sersic_index, guess_offset])
        fit = fit_function(sersic_function, self.x_values, self.y_values, self.y_uncertainties, initial_guess=initial_guess)
        return fit

    def fit_sersics(self) -> Tuple:
        """Fitting both a free index and a n=1 sersic fits."""
        fit_sersic_free = self.fit_sersic(sersic)
        fit_sersic_n1 = self.fit_sersic(sersic_n1)
        return fit_sersic_free, fit_sersic_n1
