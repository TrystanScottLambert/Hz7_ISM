#############################
#
# Surface Density Profile 
# Trystan Scott Lambert 
# 
############################

# import packages
import numpy as np 
import pylab as plt
from astropy.io import fits
import astropy.units as u 
from astropy.wcs import WCS
from astropy import modeling
from astropy.convolution import convolve, Box1DKernel
from astropy.modeling.models import Sersic1D
import TEST_Surface_Profile as TEST
from RegscorePy import * 
from scipy.optimize import curve_fit
from tqdm import tqdm
from astropy.visualization import simple_norm

def log_uncertainty(value,uncertainty):
	cut = np.where(value==0)[0]
	local_value = value.copy()
	local_value[cut] = 0.000001
	return 0.434*(uncertainty/value)

# from the rms region in CASA
sigma_val = 6.1493023326e-05


# sersic profile 
def sersic(R, Ie, Re, m):
	bm = 2.0*m - 0.324
	return Ie * np.exp(-bm * ((R/Re)**(1.0/m) - 1.0))

# sersic profile with an index of 1
def sersic1(R, Ie, Re):
	return sersic(R,Ie,Re,1)

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
	print(annuli_area)

	if error==False:
		return sum_val, len(xs)*arcsec_area
	elif error == True:
		return sum_val,annuli_area,sigma_val*(np.sqrt(annuli_area/AreaBeam))/annuli_area

'''def sum_region(center,radius_min,radius_max,data_array,error=False):
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
	print(annuli_area)

	if error==False:
		return sum_val, len(xs)*arcsec_area
	elif error == True:
		print('FUCK ANNULU',annuli_area)
		print('FUCK SIGMA_VAL: ', val)
		return sum_val,annuli_area,np.std(val)/annuli_area#*sigma_val*(np.sqrt(annuli_area/AreaBeam))/annuli_area'''

def calc_rms(center,radius_min,radius_max,data_array):
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

	rms = np.nanstd(val)
	return rms 


def get_xyz(image_data):
	x,y,z = [],[],[]
	for i in range(len(image_data)):
		for j in range(len(image_data[0])):
			z.append(image_data[i][j])
			x.append(j)
			y.append(i)
	return np.array(x),np.array(y),np.array(z)

def make_hist_xy(image_data,method='sum'):
	x,y,z = get_xyz(image_data)
	unique_x,unique_y = np.unique(x),np.unique(y)
	x_vals,y_vals = [],[]
	 
	for x_val in unique_x:
		print('HERE')
		cut = np.where(x==x_val)[0]
		x_vals.append(z[cut])

	for y_val in unique_y:
		cut = np.where(y==y_val)[0]
		y_vals.append(z[cut])

	if method=='sum':
		final_x = np.array([np.nansum(x_val) for x_val in x_vals])
		final_y = np.array([np.nansum(y_val) for y_val in y_vals])
	
	elif method=='mean':
		final_x = np.array([np.nanmean(x_val) for x_val in x_vals])
		final_y = np.array([np.nanmean(y_val) for y_val in y_vals])

	elif method=='max':
		final_x = np.array([np.max(x_val) for x_val in x_vals])
		final_y = np.array([np.max(y_val) for y_val in y_vals])
	
	else:
		print('Method choice is incorrect. Use "sum", "mean", or "max"')

	return unique_x,final_x,unique_y,final_y


#infile = 'hz7_full_mom0.fits'
#for value in tqdm(range(1,17)):

value=0
infile = 'HZ7_mom_mask2.integrated.fits'
infile = 'HZ7_mom_includepix.integrated.fits'
infile = 'test_mom.fits'
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
arcsec_per_pixel = degree_per_pixel*3600  						# convert to arcseconds per pixel 
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

#center_pix = (155,139)
#center_pix = (148,135)
center_pix = TEST.Locate_Center(data,plot=True)
#center_pix = (49,40)
#center_pix = (len(data)/2,len(data[0])/2)    # define where we center the expansion

wcs = WCS(moment0[0].header,naxis=2)  					# getting the world co-ordinate system  (for plotting mostly)
center_pix_world = wcs.wcs_pix2world([center_pix],0)    # getting the wcs of the center pixels

fig = plt.figure(figsize=(20,16))   			# plotting the circluar apetures on to the images
ax = fig.add_subplot(121,projection=wcs)
for radii in radiis:	
	circle = plt.Circle(center_pix,radii,facecolor='None',color='r',fill=False)
	ax.add_artist(circle)
	ax.set_xlabel('Dec [dd:mm:ss]',fontsize=30)
	ax.set_ylabel('RA [hh:mm:ss]',fontsize=30)
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
	new_radiis.append((radiis[i+1]+radiis[i])/2)
areas = np.array(areas)
sums = np.array(sums)
uncertainties = np.array(uncertainties)
new_radiis = np.array(new_radiis)
new_radiis_arc = new_radiis*arcsec_per_pixel

vals = sums/areas
vals_uncertainties = uncertainties/areas

full_gauss_radii = np.append(-1*new_radiis_arc[::-1],new_radiis_arc)
full_gauss_vals = np.append(vals[::-1],vals)
full_guass_vals_uncertainties = np.append(vals_uncertainties[::-1],vals_uncertainties)


ax1 = fig.add_subplot(122)
ax1.errorbar(full_gauss_radii,full_gauss_vals,yerr=full_guass_vals_uncertainties,fmt='ro',lw=3)
ax1.set_xlabel('Radius ["]',fontsize=30)
ax1.set_ylabel(r'$\Sigma_{C_{II}}$/${\rm area}$',fontsize=30)
ax1.tick_params(axis='both',labelsize=20)
ax1.minorticks_on()
ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
ax1.tick_params(which='minor', length=4, direction='in')
#ax1.set_xlim(left=0)
ax1.grid()


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

ax1.plot(x,fitted_model(x),label='Gaussian Fit',lw=3,color='k')
ax1.plot(x,gg_fit(x),label='2 Gaussian Fit',lw=3,color='m')
ax1.plot(x,part_a_gauss(x),lw=2,ls='--',color='r')
ax1.plot(x,part_b_gauss(x),lw=2,ls='--',color='b')	

plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)

sers_x = np.linspace(0,2,1000)
popt,pcov = curve_fit(sersic,plotting_radiis_arcsec,vals)
ax1.plot(sers_x,sersic(sers_x,*popt),lw=3,color='g',label=f'n = {round(popt[-1],3)}')
sersic_aic = aic.aic(vals,sersic(plotting_radiis_arcsec,*popt),3)

popt1,pcov1 = curve_fit(sersic1,plotting_radiis_arcsec,vals)
ax1.plot(sers_x,sersic1(sers_x,*popt1),lw=2,color='k',label=f'n = 1')
sersic_aic = aic.aic(vals,sersic1(plotting_radiis_arcsec,*popt1),3)

ax1.set_xlim(left=0)
ax1.errorbar(full_gauss_radii,full_gauss_vals,yerr=full_guass_vals_uncertainties,fmt='ro',lw=3)
ax1.set_xlabel('Radius ["]',fontsize=30)
ax1.set_ylabel(r'$\Sigma_{C_{II}}$/${\rm area}$',fontsize=30)
ax1.tick_params(axis='both',labelsize=20)
ax1.minorticks_on()
ax1.legend(fontsize=20)
ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
ax1.tick_params(which='minor', length=4, direction='in')
plt.show()

sers_x_c, popt_c, popt1_c = sers_x, popt, popt1
CII_radii, CII_vals, CII_undertainties = full_gauss_radii, full_gauss_vals, full_guass_vals_uncertainties

############ Reading in the UV image #####################


print('AIC Scores: ')
print(f'1 Gauss: {gauss_1d_aic}')
print(f'2 Gauss: {gauss_2d_aic}')
print(f'Sersic: {sersic_aic}')
print('\n\n')
print('Sersic Best Fit: ')
print(popt)



########################################################################################################################################################

# from the rms region in CASA
sigma_val = 0.002128163

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
	print(annuli_area)

	if error==False:
		return sum_val, len(xs)*arcsec_area
	elif error == True:
		return sum_val,annuli_area,np.sqrt(np.power(sigma_val*np.sqrt(len(val)),2) + (sum_val/4056.))/annuli_area #
		#return sum_val,annuli_area,sigma_val*np.sqrt(len(val))/annuli_area




 #################################### Preparing the different filters ##################################################


filters =['gaia_corrected_rescaled_acs.fits',
  'resampled_gaia_corrected_hst_13641_07_wfc3_ir_f105w_sci (1).fits',
  'resampled_gaia_corrected_hst_13641_07_wfc3_ir_f160w_sci (1).fits',
  'resampled_gaia_corrected_hst_13641_07_wfc3_ir_f125w_sci.fits',
  'resampled_gaia_corrected_hst_13641_07_wfc3_ir_total_sci.fits']

colors = ['g','r','m','c','y']
labels = ['ACS','105w','160w','125w','total']
value=0

fig = plt.figure()
ax1 = fig.add_subplot(111)	
ax1.plot(sers_x_c,sersic(sers_x_c,*popt_c)/np.max(sersic(sers_x_c,*popt_c)),lw=3,color='b',label=f'[CII], n = {round(popt_c[-1],3)} \n reff = {round(popt_c[1],3)}')
ax1.errorbar(CII_radii,CII_vals/np.max(CII_vals),yerr=CII_undertainties/np.max(CII_vals),fmt='bo',lw=3)

f = open('CII_SurfaceProfile.txt','w')
f.write('# CII Surface Profile data \n')
f.write('# Radii, normalized flux, flux uncertainty \n')
for j in range(len(CII_radii)):
	f.write(f'{CII_radii[j]} {CII_vals[j]/np.max(CII_vals)} {CII_undertainties[j]/np.max(CII_vals)} \n')
f.close()

### Looping over all the filters ###
for filt in filters:
	infile = filt
	moment0 = fits.open(infile)
	data = moment0[0].data   # open file and strip data into the file
	if filt != 'gaia_corrected_rescaled_acs.fits':
		data = data/2.5    # divide by the gain to get counts per second. 
	header = moment0[0].header
	PixSize = header['CDELT2']*3600
	degree_per_pixel = np.abs(float(moment0[0].header['CDELT1']))   # figure out the degrees per pixel from the headers 
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

	'''fig = plt.figure(figsize=(20,16))   			# plotting the circluar apetures on to the images
	ax = fig.add_subplot(121,projection=wcs)
	for radii in radiis:	
		circle = plt.Circle(center_pix,radii,facecolor='None',color='r',fill=False)
		ax.add_artist(circle)
		ax.set_xlabel('Dec [dd:mm:ss]',fontsize=30)
		ax.set_ylabel('RA [hh:mm:ss]',fontsize=30)
	norm = simple_norm(data,'log')
	if filt == 'gaia_corrected_rescaled_acs.fits':
		ax.imshow(data,norm=norm,vmin=0,vmax=0.03)
	else:
		ax.imshow(data)'''


	areas = []
	sums = []
	new_radiis = []
	uncertainties =[]
	for i in range(len(radiis)-1):
		current_sum,current_area,current_uncertainty = sum_region(center_pix,radiis[i],radiis[i+1],data,error=True)
		sums.append(current_sum)
		areas.append(current_area)
		uncertainties.append(current_uncertainty)
		new_radiis.append((radiis[i+1]+radiis[i])/2)
	areas = np.array(areas)
	sums = np.array(sums)
	uncertainties = np.array(uncertainties)
	new_radiis = np.array(new_radiis)
	new_radiis_arc = new_radiis*arcsec_per_pixel

	vals = sums/areas
	vals_uncertainties = uncertainties/areas
	if filt != 'gaia_corrected_rescaled_acs.fits':
		vals = vals - np.min(vals)
		vals_uncertainties = vals_uncertainties - np.min(vals)

	full_gauss_radii = np.append(-1*new_radiis_arc[::-1],new_radiis_arc)
	full_gauss_vals = np.append(vals[::-1],vals)
	full_guass_vals_uncertainties = np.append(vals_uncertainties[::-1],vals_uncertainties)

	'''ax1 = fig.add_subplot(122)
	ax1.errorbar(full_gauss_radii,full_gauss_vals,yerr=full_guass_vals_uncertainties,fmt='ro',lw=3)
	ax1.set_xlabel('Radius ["]',fontsize=30)
	ax1.set_ylabel(r'$\Sigma_{C_{II}}$/${\rm area}$',fontsize=30)
	ax1.tick_params(axis='both',labelsize=20)
	ax1.minorticks_on()
	ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	ax1.tick_params(which='minor', length=4, direction='in')
	#ax1.set_xlim(left=0)
	ax1.grid()'''


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

	'''ax1.plot(x,fitted_model(x),label='Gaussian Fit',lw=3,color='k')
	ax1.plot(x,gg_fit(x),label='2 Gaussian Fit',lw=3,color='m')
	ax1.plot(x,part_a_gauss(x),lw=2,ls='--',color='r')
	ax1.plot(x,part_b_gauss(x),lw=2,ls='--',color='b')	

	plt.show()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)'''

	sers_x = np.linspace(0,2,1000)
	popt,pcov = curve_fit(sersic,plotting_radiis_arcsec,vals)
	#ax1.plot(sers_x,sersic(sers_x,*popt),lw=3,color='g',label=f'n = {round(popt[-1],3)}')
	sersic_aic = aic.aic(vals,sersic(plotting_radiis_arcsec,*popt),3)

	popt1,pcov1 = curve_fit(sersic1,plotting_radiis_arcsec,vals)
	#ax1.plot(sers_x,sersic1(sers_x,*popt1),lw=2,color='k',label=f'n = 1')
	sersic_aic = aic.aic(vals,sersic1(plotting_radiis_arcsec,*popt1),3)

	'''ax1.set_xlim(left=0)
	ax1.errorbar(full_gauss_radii,full_gauss_vals,yerr=full_guass_vals_uncertainties,fmt='ro',lw=3)
	ax1.set_xlabel('Radius ["]',fontsize=30)
	ax1.set_ylabel(r'$\Sigma_{UV}$/${\rm area}$',fontsize=30)
	ax1.tick_params(axis='both',labelsize=20)
	ax1.minorticks_on()
	ax1.legend(fontsize=20)
	ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	ax1.tick_params(which='minor', length=4, direction='in')
	plt.show()'''

############ Reading in the UV image #####################


	print('AIC Scores: ')
	print(f'1 Gauss: {gauss_1d_aic}')
	print(f'2 Gauss: {gauss_2d_aic}')
	print(f'Sersic: {sersic_aic}')
	print('\n\n')
	print('Sersic Best Fit: ')
	print(popt)


	

	sers_x = np.linspace(0,2,1000)
	popt,pcov = curve_fit(sersic,plotting_radiis_arcsec,vals)
	
	sersic_aic = aic.aic(vals,sersic(plotting_radiis_arcsec,*popt),3)

	popt1,pcov1 = curve_fit(sersic1,plotting_radiis_arcsec,vals)
	# Chose to plot n best fit of n = 1. 
	#ax1.plot(sers_x,sersic1(sers_x,*popt1)/np.max(sersic1(sers_x,*popt1)),lw=3,color=colors[value],label=f'{labels[value]}, n = {1} \n reff = {round(popt1[1],3)}')
	ax1.plot(sers_x,sersic(sers_x,*popt)/np.max(sersic(sers_x,*popt)),color=colors[value],lw=3,label=f'{labels[value]}, n = {round(popt[-1],3)} \n reff = {round(popt[1],3)}')
	
	ax1.errorbar(full_gauss_radii,full_gauss_vals/np.max(full_gauss_vals),color=colors[value],yerr=full_guass_vals_uncertainties/np.max(full_gauss_vals),fmt='o',lw=3)
	sersic_aic = aic.aic(vals,sersic1(plotting_radiis_arcsec,*popt1),3)
	if filt == 'resampled_gaia_corrected_hst_13641_07_wfc3_ir_f160w_sci (1).fits' :
		f = open('F160W_SurfaceProfile.txt','w')
		f.write('# F160W Surface Profile Data \n')
		f.write('# Radii, normalzied flux, flux uncertainty \n')
		for j in range(len(full_gauss_radii)):
			print('WIRITING')
			f.write(f'{full_gauss_radii[j]} {full_gauss_vals[j]/np.max(full_gauss_vals)}  {full_guass_vals_uncertainties[j]/np.max(full_gauss_vals)} \n')
		f.close()
	#ax1.semilogy()
	ax1.set_xlim(left=0)
	ax1.set_xlabel('Radius ["]',fontsize=30)
	ax1.set_ylabel(r'Normalized $\Sigma$ flux/${\rm area}$',fontsize=30)
	ax1.tick_params(axis='both',labelsize=20)
	ax1.minorticks_on()
	ax1.legend(fontsize=20)
	ax1.tick_params(which='both', width=2) #this has to be a separate line because width can't be set when using which='minor'
	ax1.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	ax1.tick_params(which='minor', length=4, direction='in')
	value+=1
plt.show()