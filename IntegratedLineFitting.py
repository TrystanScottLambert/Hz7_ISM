"""
Class to fit guassian lines to [CII] emission 
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import astropy.units as u 
import pylab as plt
import Cosmology as cosmo

def gaussian_function(x, a, x0, sigma, c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def fit_function(function,xData,yData,yUncertainty,nstd=3):
	popt,pcov = curve_fit(function,xData,yData,sigma=yUncertainty)
	perr = np.sqrt(np.diag(pcov))
	poptUp = popt + nstd * perr 
	poptDown = popt - nstd * perr
	return popt,pcov,perr,poptUp,poptDown


def calculate_luminosity(integrated_flux, redshift):
	'''
	Convert integrated flux into a physical L_{/odot} luminosity using:
	Solomon and Vanden Bout
	https://ui.adsabs.harvard.edu/abs/2005ARA%26A..43..677S/abstract
	'''
	rest_freq = 1900.54 #Ghz
	const = 1.04e-3
	lumonosity_distance = cosmo.Cosmology().calc_luminosity_distance(redshift)
	return const * rest_freq * integrated_flux * ((1 + redshift)**-1) * (lumonosity_distance**2)
		

def calculate_luminosity_uncertainty(luminosity, integrated_flux, redshift, u_integrated_flux, u_redshift):
	terms = (u_integrated_flux / integrated_flux)**2 + (u_redshift / redshift)**2
	return luminosity * np.sqrt(terms)


class GuassianFit():
	def __init__(self, x, y, x_unit, y_unit, nstd = 1, y_uncertainty = None):
		self.x_unit = x_unit
		self.y_unit = y_unit
		self.x_conversion_factor = self._get_conversion_to_kms(x_unit)
		self.y_conversion_factor = self._get_conversion_to_mJy(y_unit)
		self.x = x * self.x_conversion_factor 
		self.y = y * self.y_conversion_factor
		if y_uncertainty != None:
			self.y_uncertainty = y_uncertainty * self.y_conversion_factor
		else:
			self.y_uncertainty = y_uncertainty
		self.nstd = 1
		self.fit_params = self.fit_params()


	def _get_conversion_to_kms(self, x_unit):
		return x_unit.to(u.km / u.s)

	def _get_conversion_to_mJy(self, y_unit):
		return y_unit.to(u.mJy)

	def fit_params(self):
		fit_parameters = fit_function(
			gaussian_function,
			self.x,
			self.y, 
			self.y_uncertainty,
			self.nstd,
			)

		self.background = fit_parameters[0][-1]

		'''fit_parameters = fit_function(
			gaussian_function,
			self.x,
			self.y - self.background, 
			self.y_uncertainty,
			self.nstd,
			)'''
		return fit_parameters

	@property
	def a(self):
		return (self.fit_params[0][0], self.fit_params[2][0])
	
	@property
	def mean(self):
		return (self.fit_params[0][1], self.fit_params[2][1])

	@property
	def sigma(self):
		return (self.fit_params[0][2], self.fit_params[2][1])
	
	@property
	def FWHM(self):
		const = 2 * np.sqrt(2 * np.log(2))
		fwhm = const * self.sigma[0]
		fwhm_error = const * self.sigma[1]
		return (fwhm, fwhm_error)

	@property
	def integral(self):
		const = np.sqrt(2 * np.pi)
		integral = const * self.a[0] * self.sigma[0]
		error_term = np.sqrt(
			(self.a[1]/self.a[0])**2 + (self.sigma[1]/self.sigma[0])**2
			)
		error = const * error_term
		return (integral / 1000, error_term / 1000)

	@property
	def fwfm(self):
		"""The maximum range theoretically containing all of the emission."""
		return (4 * self.sigma[0], 4 * self.sigma[1])

	def plot_line(self, outfile_name):
		fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
		plt.step(self.x, self.y - self.background, where = 'mid', color = 'k', lw = 1)
		plt.axhline(0,color='k',lw=1,ls=':')
		plt.axvline(self.mean[0], color = 'k', lw=1, ls=':')
		model_x = np.linspace(np.min(self.x), np.max(self.x), 1000)
		plt.plot(model_x,gaussian_function(model_x, *self.fit_params[0]),color='k',lw=2)

		binwidth = np.abs(self.x[1]-self.x[0])
		y_min, y_max = np.min(self.y) - 0.5, np.max(self.y) + 0.5
		plt.ylim(y_min,y_max)
		
		y_range = y_max-y_min
		zero_point = np.abs(float(y_min)/y_range)

		for i in range(len(self.x)):
			if self.y[i] > 0:
				plt.axvspan(self.x[i]-binwidth/2,self.x[i]+binwidth/2,ymax=np.abs(float(self.y[i]-y_min - self.background))/y_range,ymin =zero_point,color='yellow')
			elif self.y[i] < 0:	
				plt.axvspan(self.x[i]-binwidth/2,self.x[i]+binwidth/2,ymin=np.abs(float(self.y[i]-y_min  - self.background))/y_range,ymax =zero_point,color='yellow')



		plt.xlabel('mJy', fontsize = 12)
		plt.ylabel(r'km s$^{-1}$', fontsize = 12)
		plt.tick_params(axis = 'both',labelsize = 10)
		plt.minorticks_on()
		plt.tick_params(which = 'both', width = 2) 
		plt.tick_params(which = 'major', length = 8, direction = 'in') 
		plt.tick_params(which = 'minor', length = 4, direction = 'in')
		plt.legend(fontsize = 8, frameon = False)
		plt.savefig(outfile_name, bbox_inches = 'tight', transparent = False)



if __name__ == '__main__':
	infile = 'data/JVM_Corrected_Integrated_Spectrum.txt'
	infile = '/home/trystan/Desktop/DELETEME.txt'
	x,y = np.loadtxt(infile,unpack=True)
	fit = GuassianFit(x, y, u.km / u.s, u.Jy)
	fit.plot_line('DELETEME.png')
