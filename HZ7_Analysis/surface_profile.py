"""Module used to calculate the surface density profiles """

from abc import abstractmethod
import numpy as np 
import pylab as plt
from astropy.io import fits 
from scipy.optimize import curve_fit
from typing import Tuple
from convolution import get_hdu_data_index
from calc_channels import calc_rms
import astropy

ARCSECONDS_IN_DEGREE = 3600

def log_uncertainty(value,uncertainty):
	cut = np.where(value==0)[0]
	local_value = value.copy()
	local_value[cut] = 0.000001
	return 0.434*(uncertainty/value)

def sersic(R, Ie, Re, m):
	bm = 2.0*m - 0.324
	return Ie * np.exp(-bm * ((R/Re)**(1.0/m) - 1.0))

def sersic1(R, Ie, Re):
	return sersic(R,Ie,Re,1)

def gaussianFunction(x, a, x0, sigma, c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def fitFunction(function,xData,yData,yUncertainty,nstd=3):
	popt,pcov = curve_fit(function,xData,yData,sigma=yUncertainty)
	perr = np.sqrt(np.diag(pcov))
	poptUp = popt + nstd * perr 
	poptDown = popt - nstd * perr
	return popt,pcov,perr,poptUp,poptDown

def sum_annulus(center, data, radius_min, radius_max):
	xs = []
	ys = []
	for i in range(len(data)):
		for j in range(len(data[0])):
			xs.append(i)
			ys.append(j)
	xs = np.array(xs)
	ys = np.array(ys)
	rs = np.sqrt((center[0] - xs)**2 + (center[1] - ys)**2)
	cut = np.where((rs <= radius_max) & (rs >= radius_min))[0]
	xs, ys, rs = xs[cut], ys[cut], rs[cut]

	val = []
	for i in range(len(xs)):
		val.append(data[ys[i]][xs[i]])

	sum_val = np.nansum(val)
	number_of_pixels = len(xs)
	return sum_val, number_of_pixels

def get_arcsec_per_pix(header_object: astropy.io.fits.header.Header) -> float:
	"""Function to get the CDELT or CD2_2 value from a header"""
	try:
		deg_per_pix = header_object['CD2_2']
	except KeyError:
		deg_per_pix = header_object['CDELT2']
	return deg_per_pix * ARCSECONDS_IN_DEGREE



class FitsImage():
	def __init__(self, fits_file):
		self.hdulist = fits.open(fits_file)
		hdu_index = get_hdu_data_index(self.hdulist)
		self.data = self.hdulist[hdu_index].data
		self.header = self.hdulist[hdu_index].header
		self.arcsec_per_pixel = get_arcsec_per_pix(self.header)

	@abstractmethod
	def sum_region(self, center: Tuple[int, int], radius_min: float,
	 radius_max: float, sigma: float) -> Tuple[float, float, float]:
		"""Method to calculate the sum of a region."""
		pass

	@abstractmethod
	def _calculate_anulus_error(self):
		"""Method which will calculate the error of the annulus."""
		pass

	def sum_regions(self, center, radiis, sigma, offset=False):
		areas = []
		sums = []
		uncertainties = []
		new_radiis = []
		for i in range(len(radiis)-1):
			current_sum, current_area, current_uncertainty = self.sum_region(center, radiis[i], radiis[i+1], sigma,)
			sums.append(current_sum)
			areas.append(current_area)
			uncertainties.append(current_uncertainty)
			if i == 0:
				new_radiis.append(0)
			else:
				new_radiis.append((radiis[i+1]+radiis[i])/2)
		areas = np.array(areas)
		sums = np.array(sums)
		uncertainties = np.array(uncertainties)
		new_radiis = np.array(new_radiis)

		new_radiis_arc = new_radiis * self.arcsec_per_pixel

		vals = sums/areas
		if offset is not False:
			offset = np.min(vals)
			vals = vals - offset
		return vals, uncertainties, new_radiis_arc

	def plot_annuli(self, center, radiis):
		fig = plt.figure(figsize=(20,16)) 
		ax = fig.add_subplot()
		for radii in radiis:	
			circle = plt.Circle(center,radii,facecolor='None',color='r',fill=False)
			ax.add_artist(circle)
		ax.imshow(self.data)
		plt.show()


class Moment0(FitsImage):
	def __init__(self, moment0_file_name):
		FitsImage.__init__(self, moment0_file_name)
		self.bmaj = self.header['BMAJ'] * ARCSECONDS_IN_DEGREE
		self.bmin = self.header['BMIN'] * ARCSECONDS_IN_DEGREE	
		self.beam_area_arcsec = 2*(np.pi*self.bmaj*self.bmin/(8.0*np.log(2.0)))
		self.beam_area_pixels = self.beam_area_arcsec/self.arcsec_per_pixel**2
		self.data = self.data / self.beam_area_pixels 

	def sum_region(self, center, radius_min, radius_max, sigma):
		sum_val, number_of_pixels = sum_annulus(center, self.data, radius_min, radius_max)
		annuli_area = number_of_pixels * (self.arcsec_per_pixel**2)
		error = self._calculate_annulus_error(sigma, annuli_area, self.beam_area_arcsec)
		return sum_val, annuli_area, error

	def _calculate_annulus_error(self, sigma, annulus_area, beam_area_arcsec):
		term = np.sqrt(annulus_area / beam_area_arcsec)
		error = sigma * term / annulus_area
		return error


class HubbleImage(FitsImage):
	def __init__(self, hubble_file_name):
		FitsImage.__init__(self, hubble_file_name)
		self.exptime = 2611.737182
		self.data = (self.exptime * self.data) / 2.5 

	def _calculate_annulus_error(self, sigma, number_pixels, annulus_sum, annulus_area):
		term1 = (sigma*np.sqrt(number_pixels))**2
		term2 = (annulus_sum / self.exptime)
		error = np.sqrt(term1 + term2)/annulus_area 
		return error 

	def sum_region(self, center, radius_min, radius_max, sigma):
		sum_val, number_of_pixels = sum_annulus(center, self.data, radius_min, radius_max)
		annuli_area = number_of_pixels * self.arcsec_per_pixel**2
		error = self._calculate_annulus_error(sigma, number_of_pixels, sum_val, annuli_area)
		return sum_val, annuli_area, error


def calculate_pixel_radii(binwidth_arcsec, max_arcsec, pixel_conversion):
	binwidth_pixel = round(binwidth_arcsec / pixel_conversion)
	max_pix = round(max_arcsec / pixel_conversion)
	radiis = np.arange(0, max_pix + binwidth_pixel, binwidth_pixel)
	return radiis


def set_up_plot():
	fig = plt.figure(figsize = (3.54*2, 3.54), dpi = 600)
	ax_main = fig.add_subplot(111)
	plt.tick_params(axis='both',labelsize=10)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4, direction='in')
	return fig, ax_main

def finish_plot(file_name, axis):
	axis.set_xlabel('Radius [kpc]', fontsize = 12)
	axis.set_ylabel(r'Normalized flux/area', fontsize = 12)
	plt.tick_params(axis='both',labelsize=10)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4, direction='in')
	axis.legend(fontsize=8,frameon=False)
	plt.savefig(file_name, bbox_inches = 'tight', transparent = False)
	plt.show()


def plot_CII(figure, axis, arcseconds, y_vals, y_uncertainties, sersic_params, sersic_params_upper, sersic_params_lower,normalize = False):
	x = np.linspace(0, 2, 1000)

	if normalize == True:
		norm_factor = np.max(sersic(x, *sersic_params))
	elif normalize == False:
		norm_factor = 1 
	else:
		assert ValueError

	axis.errorbar(arcseconds * arcsecondFactor, y_vals / norm_factor, yerr = y_uncertainties / norm_factor, fmt = 'ro', label = '[CII] Data', capsize = 3, ms = 2)
	axis.plot(x * arcsecondFactor, sersic(x, *sersic_params) / norm_factor, label = '[CII] Sersic Fit', color = 'k')
	axis.fill_between(x * arcsecondFactor, sersic(x, *sersic_params_lower) / norm_factor,sersic(x, *sersic_params_upper) / norm_factor, alpha = 0.25)

	ax_arcsecond = axis.twiny()
	ax_arcsecond.plot(x, sersic(x, *sersic_params) / norm_factor, alpha=0.)
	ax_arcsecond.set_xlabel('Radius ["]', fontsize = 12)


def plot_UV(figure,axis, arcseconds, y_vals, y_uncertainties,sersic_params,label = None,color = 'k',normalize = False):
	x = np.linspace(0, 2, 1000)
	if normalize == True:
		norm_factor = np.max(sersic(x, *sersic_params))
	elif normalize == False:
		norm_factor = 1 
	else:
		assert ValueError


	axis.errorbar(arcseconds * arcsecondFactor, y_vals / norm_factor, yerr = y_uncertainties / norm_factor, fmt = 'o', capsize = 3, ms = 2,label = label,color = color)
	axis.plot(x * arcsecondFactor,sersic(x, *sersic_params) / norm_factor, color = color)


if __name__ == '__main__':
	RC = Moment0('data/HZ7_integrated.fits')
	sigma = calc_rms(RC.data, 20)
	binwidth_arcsec = 0.2
	max_arcsec = 1.8
	center_pix = (155,139)
	arcsecondFactor = 6.117 # kpc/arcsecond from Ned Wright Calculator
	x = np.linspace(0, 2, 1000)
	
	radiis = calculate_pixel_radii(binwidth_arcsec, max_arcsec, RC.arcsec_per_pixel)
	cii_sums, cii_uncertainties, plotting_radiis_arc = RC.sum_regions(center_pix, radiis, sigma)
	nstd = 1
	CII_params_free = fitFunction(sersic, plotting_radiis_arc, cii_sums, cii_uncertainties,nstd)
	CII_params_one = fitFunction(sersic1, plotting_radiis_arc, cii_sums, cii_uncertainties,nstd)

	print(f'CII: Re = {CII_params_free[0][1]} +- {CII_params_free[2][1]}')	
	print(f'CII: Re = {CII_params_free[0][1] * arcsecondFactor} +- {CII_params_free[2][1] * arcsecondFactor}')

	fig, ax = set_up_plot()

	plot_CII(fig,ax,plotting_radiis_arc, cii_sums, cii_uncertainties, CII_params_free[0], CII_params_free[3], CII_params_free[4],normalize = True)

	
	prefix = 'data/ASCImages/convolved/'
	filters =[
	'hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved.fits',
	'hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected_convolved.fits',
	'hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected_convolved.fits',
	]

	colors = ['g', 'b', 'm']
	labels = [r'$Y_{\rm 105W}$',r'$J_{\rm 125W}$',r'$H_{\rm 160W}$']

	filters = [prefix+filt for filt in filters]
	center_pix = (1150,1046)
	hubble_sigma = 0.39
	effective_radiis = []
	for i, _filter in enumerate(filters):
		HI = HubbleImage(_filter)
		radiis = calculate_pixel_radii(0.12, 1.5, HI.arcsec_per_pixel)

		HI.sum_regions(center_pix, radiis, sigma)
		vals, uncertainties, radiis = HI.sum_regions(center_pix, radiis, hubble_sigma, offset = True)
		print('#############')
		print('filter', _filter)
		print(vals)
		print(radiis)
		print('############')

		popt, pcov = curve_fit(sersic, radiis, vals, sigma = uncertainties)
		effective_radiis.append(popt[1])
		perr = np.sqrt(np.diag(pcov))
		norm_factor = np.max(sersic(x, *popt))

		print(f'{_filter.split("_f")[-1].split("w_")[0]}: Re = {popt[1]} +- {perr[1]}')
		print(f'{_filter.split("_f")[-1].split("w_")[0]} [kpc]: Re = {popt[1] * arcsecondFactor} +- {perr[1] * arcsecondFactor}')
		
		plot_UV(fig, ax, radiis, vals, uncertainties, popt,label = labels[i],normalize = True,color = colors[i],)


	ax.axvline(CII_params_free[0][1]*arcsecondFactor,ls='--',color='k',lw=1.5,alpha=0.7)
	ax.axvline(np.mean(effective_radiis)*arcsecondFactor,color='k',lw=1,alpha=0.5,ls='-.')
	finish_plot('SurfaceDensities.png', ax)
	