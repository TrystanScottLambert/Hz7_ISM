'''Script to plot the famous sersic profile plot. '''

from bdb import effective
import numpy as np
import pylab as plt
from mcmc_sersisc_fit import sersic_fit
from mcmc_sersisc1_fit import sersic1_fit
from sersic_fitter import sersic, sersic_n1
from hubble_image import HubbleImage
from moment_map import MomentMap
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 

REDSHIFT = 5.2532
hubble_images = [
    'data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits',
	'data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits',
	'data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits',
]

labels = [
    r'Y$_{105W}$',
    r'J$_{125W}$', 
    r'H$_{160W}$',
]

colors = [
    'g', 
    'b',
    'm',
]

shapes = [
    's',
    '^',
    'd',
    'o',
]

def stringarize_params(param_index, param_array, uncertainty_array):
    term_1 = str(round(param_array[param_index],3))
    term_2 = '+' + str(round(uncertainty_array[param_index][0],3))
    term_3 = '-' + str(round(uncertainty_array[param_index][1],3))
    return '$' + term_1 + '^{' + term_2 +'}_{'+ term_3 +'}$'

def stringarize_kpc(param_array, uncertainty_array, fits_object):
    er_params = np.array([param_array[1], uncertainty_array[1][0], uncertainty_array[1][1]])
    er_params_kpc = fits_object.convert_arcsecs_to_kpc(er_params, REDSHIFT)
    terms = [str(round(er_params_kpc[i].value,3)) for i in range(3)]
    return '$' + terms[0] + '^{+' + terms[1] +'}_{-'+ terms[2] +'}$'


def sersic_fits(sersic_function, sersic_fitter):
    fig = plt.figure(figsize=(3.54*2,3.54),dpi=600)
    ax = fig.add_subplot()

    bins_arcsecs = np.arange(0.05, 2, 0.1)

    hubble_radii = []
    latex_parm_strings = []
    for i, hubble_image in enumerate(hubble_images):
        filter = HubbleImage(hubble_image)
        HUBBLE_CENTER = (1150,1046)
        #RADII = np.arange(1, 20, 2)
        RADII = bins_arcsecs / filter.arcsec_per_pix

        hubble_x, y_vals, y_uncertainties, areas = filter.generate_radial_profile(HUBBLE_CENTER, RADII)
        y_plotting = (y_vals / areas) 

        #surface_brightness = filter.calculate_stmags(y_plotting)
        #surface_brightness_uncertainty = 2.5 * (y_uncertainties / (y_vals  * np.log(10)))

        hubble_x_arcsecs = filter.convert_pixels_to_arcsec(hubble_x)
        hubble_x_kpc = filter.convert_arcsecs_to_kpc(hubble_x_arcsecs, REDSHIFT)

        params, uncertainties = sersic_fitter(hubble_x_arcsecs, y_plotting, y_uncertainties)
        effective_radius_arcsec = stringarize_params(1, params, uncertainties)
        effective_radius_kpc = stringarize_kpc(params, uncertainties, filter)
        sersic_index = stringarize_params(2, params, uncertainties)

        if sersic_fitter == sersic1_fit:
            sersic_index = '1'

        latex_parm_strings.append(f'{labels[i]} & {effective_radius_arcsec} & {effective_radius_kpc} & {sersic_index}')
        hubble_radii.append(params[1])

        x_plotting_arcsec = np.linspace(0, hubble_x_arcsecs[-1], 10000)
        x_plotting_kpc = filter.convert_arcsecs_to_kpc(x_plotting_arcsec, REDSHIFT)

        offset = np.mean(y_plotting[-3:])
        norm_factor = sersic_function(0, *params) - offset
        ax.plot(x_plotting_arcsec, (sersic_function(x_plotting_arcsec, *params) - offset) / norm_factor, color = colors[i], lw=1.5)
        ax.errorbar(hubble_x_arcsecs, (y_plotting - offset) / norm_factor, yerr = y_uncertainties / norm_factor, fmt=shapes[i], capsize=3, ms = 4, color = colors[i],  label = labels[i])

    INFILE = 'data/HZ7_integrated.fits'
    moment0 = MomentMap(INFILE)
    CENTER = (154, 138)
    RADII = bins_arcsecs / moment0.arcsec_per_pix


    x_vals, y_vals, y_uncertainties, areas = moment0.generate_radial_profile(CENTER, RADII)
    y_plotting = y_vals / areas
    y_plotting_uncertainties = y_uncertainties / areas

    x_vals_arcseconds = moment0.convert_pixels_to_arcsec(x_vals)
    x_vals_kpc = moment0.convert_arcsecs_to_kpc(x_vals_arcseconds, REDSHIFT)
    radio_params, radio_uncertainties = sersic_fitter(x_vals_arcseconds, y_plotting, y_plotting_uncertainties)

    effective_radius_arcsec = stringarize_params(1, radio_params, radio_uncertainties)
    effective_radius_kpc = stringarize_kpc(radio_params, radio_uncertainties, moment0)
    sersic_index = stringarize_params(2, radio_params, radio_uncertainties)

    if sersic_fitter == sersic1_fit:
        sersic_index = '1'

    radio_string = [f'[CII] & {effective_radius_arcsec} & {effective_radius_kpc} & {sersic_index}']
    hubble_radii.append(params[1])

    x_plotting = np.linspace(0, x_vals_arcseconds[-1], 10000)
    x_plotting_kpc = moment0.convert_arcsecs_to_kpc(x_plotting, REDSHIFT)

    norm_factor = sersic_function(0, *radio_params)

    ax.plot(x_plotting, sersic_function(x_plotting, *radio_params) / norm_factor, color = 'k', lw=2.5, label = 'Sersic best fit')
    ax.errorbar(x_vals_arcseconds, y_plotting / norm_factor, yerr = y_plotting_uncertainties / norm_factor, fmt = 'ro', capsize = 5, ms = 4.5, label = '[CII] data')

    mean_hubble_radii = np.mean(hubble_radii)
    sigma_radii = np.std(hubble_radii)
    ax.axvline(np.mean(hubble_radii), ls='--', color='k', alpha=0.3)
    ax.axvspan(mean_hubble_radii - sigma_radii, mean_hubble_radii + sigma_radii, color='k', alpha=0.2)
    ax.axvline(radio_params[1], ls = '-.', lw=2, color='k')
    ax.axvspan(radio_params[1] - radio_uncertainties[1][0], radio_params[1] + radio_uncertainties[1][1], color='k', alpha=0.2)
    ax.set_ylabel('Normalalized Flux / Area [arbitary units]')
    ax.set_xlabel('Radius ["]')

    ax2 = ax.twiny()
    ax2.plot(x_plotting_kpc, sersic_function(x_plotting, *radio_params), alpha=0, color='k')
    ax2.set_xlabel('Radius [kpc]')

    plt.minorticks_on()
    plt.ylim(0,1.2)
    ax.tick_params(which='both', width=2,direction='in') 
    ax.tick_params(which='major', length=4, direction='in')
    ax2.tick_params(which='both', width=2,direction='in') 
    ax2.tick_params(which='major', length=4, direction='in')
    ax.legend(fontsize=8,frameon=False)
    plt.savefig('surface_densities.png',bbox_inches='tight',transparent=False)
    latex_strings = radio_string + latex_parm_strings
    for string in latex_strings:
        print(string + 2*'\\')
    print('\hline')
if __name__ == '__main__':
    #sersic_fits(sersic, sersic_fit)
    sersic_fits(sersic_n1, sersic1_fit)