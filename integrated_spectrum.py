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

def sersic(R, Ie, Re, m):
	bm = 2.0*m - 0.324
	return Ie * np.exp(-bm * ((R/Re)**(1.0/m) - 1.0))


def sersic1(R, Ie, Re):
	return sersic(R,Ie,Re,1)


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
model = modeling.models.Gaussian1D(amplitude=20,mean=200,stddev=200)
def print_set(x,y,show_sigma = False,save_name=None):
	fig = plt.figure(figsize=(16,9))
	#plt.ylim(-0.0002,0.0007)
	fitted_model = fitter(model,x,y)
	gy = norm.pdf(x,mean,std)
	#plt.plot(x,y,marker='o',color='k',lw=3)
	
	## Capak Plotting #############################################
	plt.step(x,y,where='mid',color='k',lw=1.5)
	#plt.scatter(x,y)
	binwidth = np.abs(x[1]-x[0])
	plt.axhline(0,color='k',lw=1,ls=':')
	y_min,y_max = -10,10
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
	plt.plot(X,fitted_model(X),color='k',lw=3)
	plt.axvline(fitted_model.mean.value,color='k',ls='dashed',lw=2)
	plt.xlabel('Radio Velocity [km/s]',fontsize=30)
	plt.ylabel('mJy',fontsize=30)
	#plt.xlabel('Channels',fontsize=30)
	#plt.ylabel('mean Jy/Beam',fontsize=30)
	
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4, direction='in')
	#print out the FWHM range and noise things
	half_FWHM = fitted_model.fwhm/2
	print('Mean (Gauss): ', fitted_model.mean.value )
	print('FWHM Range: ', fitted_model.mean.value-half_FWHM, fitted_model.mean.value + half_FWHM)
	plt.axvspan(fitted_model.mean.value-half_FWHM, fitted_model.mean.value + half_FWHM,color='g',alpha=0.3)

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
		plt.savefig(save_name)
	
	#plt.show()

#print_set(x,y,show_sigma=True,save_name='Channel_IS.png')
#print_set(x,s3,show_sigma=True,save_name='Channel_IS_3.png')
#print_set(x,s5,show_sigma=True,save_name='Channel_IS_5.png')
############ Plotting a long list of profiles using multiple differnet anuli. #######################

'''for i in range(1,11):
	infile = 'is_'+str(i)+'_sum.txt'
	print infile
	x,y = np.loadtxt(infile,unpack=True)
	arg = np.argsort(x)
	x,y = x[arg]/1000,1000*(y[arg]/corrections[arg])
	
	fitter = modeling.fitting.LevMarLSQFitter()
	model = modeling.models.Gaussian1D(amplitude=20,mean=200,stddev=200)
	print_set(x,y,save_name='is_'+str(i)+'.png')'''

############## Plotting one very nice spectrum with corrected beams to be compared to Capak. #################
infile = 'data/integratedSpectrum.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000#/corrections  #convert into mJy / pixel 

s3 = convolve(y,Box1DKernel(3))
s5 = convolve(y,Box1DKernel(5))

print_set(x,y,save_name='CII_Integrated_Spectrum.png')
#print_set(x,s3,save_name='CII_Integrated_Spectrum_3.png')
#print_set(x,s5,save_name='CII_Integrated_Spectrum_5.png')

corrections = np.array([(beam[0]*(1./0.02599)*beam[1]*(1./0.02599)*np.pi)/(4*np.log(2)) for beam in beam_data])
corrections = corrections[21:]
x,y = np.loadtxt('data/integratedSpectrum.txt',unpack=True)
arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000 #convert into mJy / pixel

#print_set(x,y,save_name='CII_flux_density.png') 


##################### PLotting old Integrated Spectrum vs the New Integrated spectrum ####################################

infile = 'data/integratedSpectrum.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000  #convert into mJy / pixel 

print_set(x,y)


infile = 'data/integratedSpectrum.txt'
x,y = np.loadtxt(infile,unpack=True)
arg = np.argsort(x)
arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000  #convert into mJy / pixel 

print_set(x,y)
plt.show()


infile = 'data/integratedSpectrum.txt'
x,y = np.loadtxt(infile,unpack=True)

arg = np.argsort(x)
x = x[arg]
y = y[arg]*1000  #convert into mJy / pixel 
model = modeling.models.Gaussian1D(amplitude=10,mean=0,stddev=20)
print_set(x,y)
plt.show()

fitter = modeling.fitting.LevMarLSQFitter()
model = modeling.models.Gaussian1D(amplitude=8,mean=60,stddev=20)  #note using the mean and standed deviation of the redshifts. 
fitted_model = fitter(model,x,y)
gauss_1d_aic = aic.aic(y,fitted_model(x),3)   # work out the aic score for the fitted model

gg_init = modeling.models.Gaussian1D(amplitude=8, mean=60, stddev=20) + modeling.models.Gaussian1D(amplitude=8, mean=60, stddev=20)
fitter = modeling.fitting.SLSQPLSQFitter()
gg_fit = fitter(gg_init, x, y)
gauss_2d_aic = aic.aic(y,gg_fit(x),6)

part_a_gauss = modeling.models.Gaussian1D(amplitude=gg_fit.amplitude_0.value,mean=gg_fit.mean_0.value, stddev=gg_fit.stddev_0.value)
part_b_gauss = modeling.models.Gaussian1D(amplitude=gg_fit.amplitude_1.value,mean=gg_fit.mean_1.value, stddev=gg_fit.stddev_1.value)

x_lin = np.linspace(0,130,1000)  #in arcseconds

plt.scatter(x,y,color='r')
plt.plot(x,y,color='k')
plt.plot(x_lin,fitted_model(x_lin),label='Gaussian Fit',lw=5,color='k')
plt.plot(x_lin,gg_fit(x_lin),label='2 Gaussian Fit',lw=3,color='m')
plt.plot(x_lin,part_a_gauss(x_lin),lw=2,ls='--',color='r',label='part_b_gauss')
plt.plot(x_lin,part_b_gauss(x_lin),lw=2,ls=':',color='b',label='part_a_gauss')	

plt.legend()
plt.tick_params(axis='both',labelsize=20)
plt.minorticks_on()
plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
plt.tick_params(which='minor', length=4, direction='in')
plt.show()

