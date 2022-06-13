import numpy as np 
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy.stats import norm 
from astropy import modeling
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit

import numpy as np 
import pylab as plt
from astropy.io import fits
import astropy.units as u 
from astropy.wcs import WCS
from astropy import modeling
from tqdm import tqdm
from astropy.convolution import convolve, Box1DKernel
from astropy.modeling.models import Sersic1D
from RegscorePy import * 
from scipy.optimize import curve_fit
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 


def gaussianFunction(x, a, x0, sigma, c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c


def fitFunction(function,xData,yData,yUncertainty,nstd=3):
	popt,pcov = curve_fit(function,xData,yData,sigma=yUncertainty)
	perr = np.sqrt(np.diag(pcov))
	poptUp = popt + nstd * perr 
	poptDown = popt - nstd * perr
	return popt,pcov,perr,poptUp,poptDown

def gaussian(x, mu, sig, c):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) + c

### Reading in Beam information so that we can correct to Jy. 
fits_file = fits.open('data/HZ7_Combined.fits')
beam_data = fits_file[1].data
corrections = np.array([(beam[0]*(1./0.02599)*beam[1]*(1./0.02599)*np.pi)/(4*np.log(2)) for beam in beam_data])

infile = 'data/integratedSpectrum.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000

s3 = convolve(y,Box1DKernel(3))
s5 = convolve(y,Box1DKernel(5))

mean, std = 80,20
noise_region_1 = [20,40]
noise_region_2 = [110,130]

def get_closest(array,value):
	vals = np.abs(array - value)
	closest_idx=np.where(vals==np.min(vals))[0]
	return closest_idx[0]

fitter = modeling.fitting.LevMarLSQFitter()
model = modeling.models.Gaussian1D(amplitude=20,mean=0,stddev=200)



def print_set(x,y,show_sigma = False,save_name=None):
	fig = plt.figure(figsize=(3.54,3.54),dpi=600)
	fitted_model = fitter(model,x,y)
	gy = norm.pdf(x,mean,std)
	
	## Capak Plotting #############################################
	plt.step(x,y,where='mid',color='k',lw=1)
	binwidth = np.abs(x[1]-x[0])
	plt.axhline(0,color='k',lw=1,ls=':')
	y_min,y_max = -1,3
	y_min, y_max = -0.5, 0.5
	plt.ylim(y_min,y_max)
	
	y_range = y_max-y_min
	zero_point = np.abs(float(y_min)/y_range)

	for i in range(len(x)):
		if y[i] > 0:
			plt.axvspan(x[i]-binwidth/2,x[i]+binwidth/2,ymax=np.abs(float(y[i]-y_min))/y_range,ymin =zero_point,color='yellow')
		elif y[i] < 0:	
			plt.axvspan(x[i]-binwidth/2,x[i]+binwidth/2,ymin=np.abs(float(y[i]-y_min))/y_range,ymax =zero_point,color='yellow')

	##########################################################
	X = np.linspace(min(x),max(x),10000)
	plt.plot(X,fitted_model(X),color='k',lw=2)
	#plt.axvline(fitted_model.mean.value,color='k',ls='dashed',lw=1)
	plt.xlabel('Offset Velocity [km/s]',fontsize=10)
	plt.ylabel('Flux Density [mJy]',fontsize=10)

	plt.tick_params(axis='both',labelsize=10)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=4, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=2, direction='in')
	#print out the FWHM range and noise things
	half_FWHM = fitted_model.fwhm/2
	print('FWHM: ', fitted_model.fwhm)
	print('Mean (Gauss): ', fitted_model.mean.value )
	print('FWHM Range: ', fitted_model.mean.value-half_FWHM, fitted_model.mean.value + half_FWHM)
	#plt.axvspan(fitted_model.mean.value-half_FWHM, fitted_model.mean.value + half_FWHM,color='g',alpha=0.3)

	noise1 = y[get_closest(x,noise_region_1[0]):get_closest(x,noise_region_1[1])]
	noise2 = y[get_closest(x,noise_region_2[0]):get_closest(x,noise_region_2[1])]
	
	noise = np.append(noise1,noise2)
	print('Noise mean: ', np.mean(noise))
	print('Noise std:', np.std(noise))
	print('1sigma level:', np.mean(noise)+np.std(noise))
	print('3sigma level:', np.mean(noise)+3*np.std(noise))
	print('5sigma level:', np.mean(noise)+5*np.std(noise))
	print('\n\n')

	if show_sigma == True:
		plt.axhline(np.std(noise),color='k',ls=':',alpha=0.3)
		plt.axhline(3*np.std(noise),color='k',ls=':',alpha=0.3)
		plt.axhline(5*np.std(noise),color='k',ls=':',alpha=0.3)
		plt.axvspan(noise_region_1[0],noise_region_1[1],color='r',alpha=0.2)
		plt.axvspan(noise_region_2[0],noise_region_2[1],color='r',alpha=0.2)

	if save_name != None:
		plt.savefig(save_name,bbox_inches='tight',frameon=False)
	

infile = 'data/JVM_Corrected_Integrated_Spectrum.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000#/corrections  #convert into mJy / pixel 

s3 = convolve(y,Box1DKernel(3))
s5 = convolve(y,Box1DKernel(5))

print_set(x,y,save_name='CII_Integrated_Spectrum.png')



infile = '/home/trystan/Desktop/5.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000#/corrections  #convert into mJy / pixel 

s3 = convolve(y,Box1DKernel(3))
s5 = convolve(y,Box1DKernel(5))

print_set(x,y,save_name='five.png')
