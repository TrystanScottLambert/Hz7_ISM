#############################
#
# Surface Density Profile 
# Trystan Scott Lambert 
# 
############################

# To do: 
# Make sure the fitting routine is acceptable
# Uncertainties of the actual parameters. 
# (possible mcmc for uncertainties)

# import packages
import numpy as np 
import pylab as plt
from astropy.io import fits
import astropy.units as u 
from astropy.wcs import WCS
from astropy import modeling
from astropy.convolution import convolve, Box1DKernel
from astropy.modeling.models import Sersic1D
from RegscorePy import * 
from scipy.optimize import curve_fit
from tqdm import tqdm
from astropy.visualization import simple_norm
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 

arcsecondFactor = 6.117 # kpc/arcsecond from Ned Wright Calculator

def log_uncertainty(value,uncertainty):
	cut = np.where(value==0)[0]
	local_value = value.copy()
	local_value[cut] = 0.000001
	return 0.434*(uncertainty/value)

# from the rms region in CASA
sigma_val = 6.1493023326e-05
sigma_val = 2.3e-2


# sersic profile 
def sersic(R, Ie, Re, m):
	bm = 2.0*m - 0.324
	return Ie * np.exp(-bm * ((R/Re)**(1.0/m) - 1.0))

# sersic profile with an index of 1
def sersic1(R, Ie, Re):
	return sersic(R,Ie,Re,1)

def gaussianFunction(x, a, x0, sigma, c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def twoGaussianFunctions(x, a1, x01, sigma1, a2, x02, sigma2):
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2)) + a2*np.exp(-(x-x02)**2/(2*sigma2**2))


def sum_region(center,radius_min,radius_max,data_array,error=False):
	xs = []
	ys = []
	for i in range(len(data_array)):
		for j in range(len(data_array[0])):
			xs.append(i)
			ys.append(j)
	xs = np.array(xs)
	ys = np.array(ys)
	rs = np.sqrt((center[0]-xs)**2+(center[1]-ys)**2)
	cut = np.where((rs<=radius_max) & (rs>=radius_min))[0]
	xs = xs[cut]
	ys = ys[cut]
	rs = rs[cut]

	val = []
	for i in range(len(xs)):
		val.append(data_array[ys[i]][xs[i]])

	sum_val = np.nansum(val)
	annuli_area = len(xs)*arcsec_area
	#print(annuli_area)

	if error==False:
		return sum_val, len(xs)*arcsec_area
	elif error == True:
		return sum_val,annuli_area,sigma_val*(np.sqrt(annuli_area/AreaBeam))/annuli_area


value=0
infile = 'data/Jorge_cut_HZ7/data/HZ7_Collapsed.fits' # 'data/HZ7_Collapsed.fits'
infile = 'data/HZ7_Collapsed.fits'
moment0 = fits.open(infile)
data = moment0[0].data[0][0]   # open file and strip data into the file
header = moment0[0].header
PixSize = header['CDELT2']*3600
BMAJ = header['BMAJ']*3600
BMIN = header['BMIN']*3600
AreaBeam = 2*(np.pi*BMAJ*BMIN/(8.0*np.log(2.0))) #arcsec^2
AreaBeamPixels = AreaBeam / (PixSize**2)
data = data/AreaBeamPixels


degree_per_pixel = np.abs(float(moment0[0].header['CDELT1']))   # figure out the degrees per pixel from the headers 
arcsec_per_pixel = degree_per_pixel*3600  					# convert to arcseconds per pixel 
arcsec_area = arcsec_per_pixel**2  # area of a single pixel     # work out the area per pixel in square arcseconds

circles=[]
binwidth_pix_arc = 0.259999999
max_arc = 2.079999999
binwidth_pix_arc = 0.2
max_arc = 1.8

binwidth_pix = round(binwidth_pix_arc/arcsec_per_pixel)#10#10
max_pix = round(max_arc/arcsec_per_pixel)



radiis = np.arange(0,max_pix+binwidth_pix,binwidth_pix)     									# list of radiis  (in pixels)
plotting_radiis = np.append(np.array([0]),np.arange(1.5*binwidth_pix,max_pix,binwidth_pix)) 	
radii_arcsec = radiis*arcsec_per_pixel  														# converting radiis into arcseconds
plotting_radiis_arcsec = plotting_radiis*arcsec_per_pixel

center_pix = (155,139)
#center_pix = (303,289) # FOR JORGE IMAGE

wcs = WCS(moment0[0].header,naxis=2)  					# getting the world co-ordinate system  (for plotting mostly)
center_pix_world = wcs.wcs_pix2world([center_pix],0)    # getting the wcs of the center pixels

fig = plt.figure(figsize=(20,16))   			# plotting the circluar apetures on to the images
ax = fig.add_subplot(121)#,projection=wcs)
for radii in radiis:	
	circle = plt.Circle(center_pix,radii,facecolor='None',color='r',fill=False)
	ax.add_artist(circle)
	ax.set_xlabel('Dec [dd:mm:ss]',fontsize=12)
	ax.set_ylabel('RA [hh:mm:ss]',fontsize=12)
ax.imshow(data)
#plt.show()

areas = []
sums = []
new_radiis = []
uncertainties =[]
for i in range(len(radiis)-1):
	current_sum,current_area,current_uncertainty = sum_region(center_pix,radiis[i],radiis[i+1],data,error=True)
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
new_radiis_arc = new_radiis*arcsec_per_pixel

vals = sums/areas
vals_uncertainties = uncertainties#/areas

full_gauss_radii = np.append(-1*new_radiis_arc[::-1][:-1],new_radiis_arc)
full_gauss_radii_red_circles = full_gauss_radii
full_gauss_vals = np.append(vals[::-1][:-1],vals)
full_guass_vals_uncertainties = np.append(vals_uncertainties[::-1][:-1],vals_uncertainties)
full_guass_vals_uncertainties_red_circles = full_guass_vals_uncertainties

offset = np.mean(full_gauss_vals[:2])
full_gauss_vals_red_circles = full_gauss_vals
full_gauss_vals = full_gauss_vals - offset



ax1 = fig.add_subplot(122)
ax1.errorbar(full_gauss_radii,full_gauss_vals,yerr=full_guass_vals_uncertainties,fmt='ro',lw=1,capsize=3,ms=1)
ax1.set_xlabel('Radius ["]',fontsize=12)
ax1.set_ylabel(r'$\Sigma_{C_{II}}$/${\rm area}$',fontsize=12)
ax1.tick_params(axis='both',labelsize=10)
ax1.minorticks_on()
ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
ax1.tick_params(which='minor', length=4, direction='in')
#ax1.set_xlim(left=0)
ax1.grid()


def fitFunction(function,xData,yData,yUncertainty,nstd=3):
	popt,pcov = curve_fit(function,xData,yData,sigma=yUncertainty)
	perr = np.sqrt(np.diag(pcov))
	poptUp = popt + nstd * perr 
	poptDown = popt - nstd * perr
	return popt,pcov,perr,poptUp,poptDown


x = np.linspace(-2,2,1000)  

popt,pcov,perr,popt_up,popt_dw = fitFunction(gaussianFunction,full_gauss_radii,full_gauss_vals,full_guass_vals_uncertainties,3)
fit = gaussianFunction(x, *popt)
fit_up = gaussianFunction(x, *popt_up)
fit_dw = gaussianFunction(x, *popt_dw)

#popt,pcov,perr,popt_up,popt_dw = fitFunction(twoGaussianFunctions,full_gauss_radii,full_gauss_vals,full_guass_vals_uncertainties,3)
#fit = twoGaussianFunctions(x, *popt)
#fit_up = twoGaussianFunctions(x, *popt_up)
#fit_dw = twoGaussianFunctions(x, *popt_dw)


plt.plot(x,fit,lw=2,color='k')
plt.fill_between(x,fit_up,fit_dw,alpha=2.5)
plt.show()


x = np.linspace(0,2,1000)
nstd = 1
poptCII,pcovCII,perrCII,popt_upCII,popt_dwCII = fitFunction(sersic,new_radiis_arc,vals,vals_uncertainties,nstd)
print('FIND ME')
print('##########################')
print(new_radiis_arc)
print(vals)
print('##########################')
popt1,pcov1,perr1,popt_up1,popt_dw1 = fitFunction(sersic1,new_radiis_arc,vals,vals_uncertainties,nstd)

fig = plt.figure(figsize=(3.54*2,3.54),dpi=600)
plt.errorbar(new_radiis_arc*arcsecondFactor,vals,yerr=vals_uncertainties,fmt='ro',label='Data',capsize=3,ms=2)
plt.plot(x*arcsecondFactor,sersic(x,*poptCII),label=r'$n = $'+f'{round(poptCII[-1],2)} ' + r'$R_{e} = $'f'{round(poptCII[1]*arcsecondFactor,2)} kpc')
plt.plot(x*arcsecondFactor,sersic1(x,*popt1),lw=2,color='k',label=r'$n=1$, $R_{e}$ = ' + f'{round(popt1[1]*arcsecondFactor,2)} kpc')
plt.fill_between(x*arcsecondFactor,sersic(x,*popt_upCII),sersic(x,*popt_dwCII),alpha=0.25,label=f'{nstd}' + r'$\sigma$')
plt.xlim(left=0)
plt.xlabel('Radius [kpc]',fontsize=12)
plt.ylabel(r'flux/area [mJy]',fontsize=12)
plt.tick_params(axis='both',labelsize=10)
plt.minorticks_on()
plt.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
plt.tick_params(which='minor', length=4, direction='in')
plt.legend(fontsize=8,frameon=False)
plt.savefig('SurfaceDensity.png',bbox_inches='tight',transparent=False)
plt.show()


################################################################################################################
# from the rms region in CASA
sigma_val = 0.002128163
EXPTIME = 4056. # exposure time in seconds. 

def sum_region(center,radius_min,radius_max,data_array,error=False):
	xs = []
	ys = []
	for i in range(len(data_array)):
		for j in range(len(data_array[0])):
			xs.append(i)
			ys.append(j)
	xs = np.array(xs)
	ys = np.array(ys)
	rs = np.sqrt((center[0]-xs)**2+(center[1]-ys)**2)
	cut = np.where((rs<=radius_max) & (rs>=radius_min))[0]
	xs = xs[cut]
	ys = ys[cut]
	rs = rs[cut]

	val = []
	for i in range(len(xs)):
		val.append(data_array[ys[i]][xs[i]])

	sum_val = np.nansum(val)
	annuli_area = len(xs)*arcsec_area
	#print(annuli_area)

	if error==False:
		return sum_val, len(xs)*arcsec_area
	elif error == True:
		return sum_val,annuli_area,np.sqrt(np.power(sigma_val*np.sqrt(len(val)),2) + (sum_val/EXPTIME))/annuli_area #
		#return sum_val,annuli_area,sigma_val*np.sqrt(len(val))/annuli_area




 #################################### Preparing the different filters ##################################################

prefix = '/home/trystan/Desktop/Work/Hz7_ISM/data/ASCImages/deresolved/'
filters =['hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected_convolved_deresolved.fits',
			'hst_13641_07_wfc3_ir_f125w_sci_gaia_corrected_convolved_deresolved.fits',
			'hst_13641_07_wfc3_ir_f160w_sci_gaia_corrected_convolved_deresolved.fits']

filters = [prefix+filt for filt in filters]

colors = ['g','b','m','c','r']
labels = [r'$Y_{\rm 105W}$',r'$J_{\rm 125W}$',r'$H_{\rm 160W}$','total']
value=0

fig = plt.figure(figsize=(3.54*1.5,3.54),dpi=600)
ax1 = fig.add_subplot(111)	
#ax1.plot(sers_x_c,sersic(sers_x_c,*popt_c)/np.max(sersic(sers_x_c,*popt_c)),lw=3,color='b',label=f'[CII], n = {round(popt_c[-1],3)} \n reff = {round(popt_c[1],3)}')
#ax1.errorbar(CII_radii,CII_vals/np.max(CII_vals),yerr=CII_undertainties/np.max(CII_vals),fmt='bo',lw=3)

#f = open('CII_SurfaceProfile.txt','w')
#f.write('# CII Surface Profile data \n')
#f.write('# Radii, normalized flux, flux uncertainty \n')
#for j in range(len(CII_radii)):
#	f.write(f'{CII_radii[j]} {CII_vals[j]/np.max(CII_vals)} {CII_undertainties[j]/np.max(CII_vals)} \n')
#f.close()

### Looping over all the filters ###
averageRadii = []
for filt in filters:
	print(filt)
	infile = filt
	moment0 = fits.open(infile)
	data = moment0[0].data   # open file and strip data into the file
	if filt != 'gaia_corrected_rescaled_acs.fits':
		data = data/2.5    # divide by the gain to get counts per second. 
	header = moment0[0].header
	PixSize = header['CDELT2']*3600
	degree_per_pixel = np.abs(float(moment0[0].header['CDELT2']))   # figure out the degrees per pixel from the headers 
	wcs = WCS(moment0[0].header,naxis=2)  
	arcsec_per_pixel = degree_per_pixel*3600  						# convert to arcseconds per pixel 
	arcsec_area = arcsec_per_pixel**2  # area of a single pixel     # work out the area per pixel in square arcseconds

	circles=[]
	binwidth_pix_arc = 0.12
	max_arc = 1.5

	binwidth_pix = round(binwidth_pix_arc/arcsec_per_pixel)#10#10
	max_pix = round(max_arc/arcsec_per_pixel)



	radiis = np.arange(0,max_pix+binwidth_pix,binwidth_pix)     									# list of radiis  (in pixels)
	plotting_radiis = np.append(np.array([0]),np.arange(1.5*binwidth_pix,max_pix,binwidth_pix)) 	
	radii_arcsec = radiis*arcsec_per_pixel  														# converting radiis into arcseconds
	plotting_radiis_arcsec = plotting_radiis*arcsec_per_pixel

	center_pix = (152,141)

						# getting the world co-ordinate system  (for plotting mostly)
	center_pix_world = wcs.wcs_pix2world([center_pix],0)    # getting the wcs of the center pixels


	areas = []
	sums = []
	new_radiis = []
	uncertainties =[]
	for i in range(len(radiis)-1):
		current_sum,current_area,current_uncertainty = sum_region(center_pix,radiis[i],radiis[i+1],data,error=True)
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
	new_radiis_arc = new_radiis*arcsec_per_pixel


	vals = sums/areas
	vals_uncertainties = uncertainties#/areas
	if filt != 'gaia_corrected_rescaled_acs.fits':
		vals = vals - np.min(vals)
		vals_uncertainties = vals_uncertainties - np.min(vals)

	full_gauss_radii = np.append(-1*new_radiis_arc[::-1][:-1],new_radiis_arc)
	full_gauss_vals = np.append(vals[::-1][:-1],vals)
	full_guass_vals_uncertainties = np.append(vals_uncertainties[::-1][:-1],vals_uncertainties)

	fitter = modeling.fitting.LevMarLSQFitter()
	model = modeling.models.Gaussian1D(amplitude=np.max(vals),mean=0.0,stddev=0.2)  #note using the mean and standed deviation of the redshifts. 
	fitted_model = fitter(model,full_gauss_radii,full_gauss_vals)
	gauss_1d_aic = aic.aic(full_gauss_vals,fitted_model(full_gauss_radii),3)   # work out the aic score for the fitted model

	gg_init = modeling.models.Gaussian1D(amplitude=np.max(vals), mean=0., stddev=0.2) + modeling.models.Gaussian1D(amplitude=np.max(vals), mean=0., stddev=0.2)
	fitter = modeling.fitting.SLSQPLSQFitter()
	gg_fit = fitter(gg_init, full_gauss_radii, full_gauss_vals)
	gauss_2d_aic = aic.aic(full_gauss_vals,gg_fit(full_gauss_radii),6)

	part_a_gauss = modeling.models.Gaussian1D(amplitude=gg_fit.amplitude_0.value,mean=gg_fit.mean_0.value, stddev=gg_fit.stddev_0.value)
	part_b_gauss = modeling.models.Gaussian1D(amplitude=gg_fit.amplitude_1.value,mean=gg_fit.mean_1.value, stddev=gg_fit.stddev_1.value)

	x = np.linspace(-2,2,1000)  #in arcseconds


	sers_x = np.linspace(0,2,1000)
	popt,pcov = curve_fit(sersic,plotting_radiis_arcsec,vals)
	print('#################')
	print('SERSIC VALS: ', vals)
	print('radii: ', plotting_radiis_arcsec)
	print('################')
	#ax1.plot(sers_x,sersic(sers_x,*popt),lw=3,color='g',label=f'n = {round(popt[-1],3)}')
	sersic_aic = aic.aic(vals,sersic(plotting_radiis_arcsec,*popt),3)

	popt1,pcov1 = curve_fit(sersic1,plotting_radiis_arcsec,vals)
	#ax1.plot(sers_x,sersic1(sers_x,*popt1),lw=2,color='k',label=f'n = 1')
	sersic_aic = aic.aic(vals,sersic1(plotting_radiis_arcsec,*popt1),3)


############ Reading in the UV image #####################


	print('AIC Scores: ')
	print(f'1 Gauss: {gauss_1d_aic}')
	print(f'2 Gauss: {gauss_2d_aic}')
	print(f'Sersic: {sersic_aic}')
	print('\n\n')
	print('Sersic Best Fit: ')
	print(popt)
	print('\n\n\n\n')
	print(f'{filt.split("_f")[-1].split("w_")[0]}: Re = {popt[1]} +- {perr[1]}')
	print(f'{filt.split("_f")[-1].split("w_")[0]}: n = {popt[-1]} +- {perr[-1]}')
	print('\n\n\n\n')


	sers_x = np.linspace(0,2,1000)
	print('HELLO PLEASE FIND ME CARALHO')
	print(f'radiis: {plotting_radiis_arcsec}, vals: {vals}, uncertainties: {vals_uncertainties}')
	popt,pcov = curve_fit(sersic,plotting_radiis_arcsec,vals,sigma=vals_uncertainties)
	perr = np.sqrt(np.diag(pcov))
	print(f'{filt.split("_f")[-1].split("w_")[0]}: {popt[0]} +- {perr[0]}')
	
	sersic_aic = aic.aic(vals,sersic(plotting_radiis_arcsec,*popt),3)

	popt1,pcov1 = curve_fit(sersic1,plotting_radiis_arcsec,vals)
	# Chose to plot n best fit of n = 1. 
	#ax1.plot(sers_x,sersic1(sers_x,*popt1)/np.max(sersic1(sers_x,*popt1)),lw=3,color=colors[value],label=f'{labels[value]}, n = {1} \n reff = {round(popt1[1],3)}')
	norm_factor = np.max(sersic(sers_x,*popt))
	ax1.plot(sers_x*arcsecondFactor,sersic(sers_x,*popt)/norm_factor,color=colors[value],lw=2,label=f'{labels[value]}')#, n = {round(popt[-1],2)} \n ' + r'$R_{e}$ = '+f'{round(popt[1],2)}')
	ax1.errorbar(full_gauss_radii*arcsecondFactor,full_gauss_vals/norm_factor,color=colors[value],yerr=full_guass_vals_uncertainties/norm_factor,fmt='o',lw=3,capsize=3,ms=2)

	sersic_aic = aic.aic(vals,sersic1(plotting_radiis_arcsec,*popt1),3)
	if filt == 'resampled_gaia_corrected_hst_13641_07_wfc3_ir_f160w_sci (1).fits' :
		f = open('F160W_SurfaceProfile.txt','w')
		f.write('# F160W Surface Profile Data \n')
		f.write('# Radii, normalzied flux, flux uncertainty \n')
		for j in range(len(full_gauss_radii)):
			print('WIRITING')
			f.write(f'{full_gauss_radii[j]} {full_gauss_vals[j]/np.max(full_gauss_vals)}  {full_guass_vals_uncertainties[j]/np.max(full_gauss_vals)} \n')
		f.close()

	value+=1
	averageRadii.append(popt[1])

averageOpticalRadius = np.mean(averageRadii)

norm_factor = np.max(sersic(sers_x,*poptCII))
ax1.plot(sers_x*arcsecondFactor,sersic(sers_x,*poptCII)/norm_factor,lw=2.5,color='k',label='[CII] Sersic Fit')#+r'$n = $'+f'{round(poptCII[-1],2)} ' + r'$R_{e} = $'f'{round(poptCII[1],2)}')
ax1.errorbar(full_gauss_radii_red_circles*arcsecondFactor,full_gauss_vals_red_circles/norm_factor,color=colors[value],yerr=full_guass_vals_uncertainties_red_circles/norm_factor,fmt='o',lw=1,capsize=3,ms=1,label='[CII] Data')
plt.fill_between(x*arcsecondFactor,sersic(x,*popt_upCII)/norm_factor,sersic(x,*popt_dwCII)/norm_factor,alpha=0.25)
ax1.set_xlim(left=0)
ax1.axvline(poptCII[1]*arcsecondFactor,ls='--',color='k',lw=1.5,alpha=0.7)
ax1.axvline(averageOpticalRadius*arcsecondFactor,color='k',lw=1,alpha=0.5,ls='-.')
ax1.set_xlabel('Radius [kpc]',fontsize=12)
ax1.set_ylabel(r'Normalized flux/${\rm area}$',fontsize=12)
ax1.tick_params(axis='both',labelsize=10)
ax1.minorticks_on()
ax1.legend(fontsize=8,frameon=False)
ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
ax1.tick_params(which='major', length=4, direction='in') #If you want tick marks on the opposite side also, add right=True
ax1.tick_params(which='minor', length=2, direction='in')
plt.savefig('SurfaceDensities.png',bbox_inches='tight',transparent=False)
plt.show()


