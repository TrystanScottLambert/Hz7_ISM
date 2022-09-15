"""Main pipeline for analysis of HZ7"""

import center_cube_on_CII
import calc_channels
import generate_moments
from plot_moment_maps import MomentMap
from align import GaiaAlignment
from convolution import generateConvolvedImages

import pylab as plt
import warnings

warnings.filterwarnings("ignore")  # mainly coming from spectral_cube

########################## INPUTS ################################

'''combined_cube = 'data/HZ7_Combined.fits'
centered_cube = 'data/HZ7_Centered.fits'
central_channel = 64
centeral_pixel_position = (156, 140)'''

combined_cube = 'data/Jorge_cut_HZ7/HZ7_COMB_CUBE_15_JvM_cut.fits'
centered_cube = 'data/Jorge_cut_HZ7/HZ7_COMB_CUBE_15_JvM_cut_centered.fits'
central_channel = 63
centeral_pixel_position = (156,140)

hubble_images = [
    'data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits',
	'data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits',
	'data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits',
]
#################################################################

# center the combined cube on the emission from [CII]
print('\n\n')
print(f'Centering the cube emission on channel {central_channel}...')
center_cube_on_CII.center_cube_on_channel(combined_cube, centered_cube, central_channel)

# calculate the channels we will use to generate the moment maps.
print()
print('Calculating the channels to collapse the centered cube...')
radio_cube = calc_channels.RadioCube(centered_cube)
start_channel, end_channel = radio_cube.calculate_channels(centeral_pixel_position)
print(f'\t Channels determined to be: {start_channel} {end_channel}')

# plot the integrated line emission profile and fit 
print()
print(f'Generating the integrated line emission profile from {centered_cube}')
line_emission = radio_cube.integrated_emission_profile
line_emission.plot_line('plots/integrated_emission_line.png')
plt.imshow(radio_cube.moment0_used_for_profile)
plt.savefig('plots/moment0_map_used_for_mask.png')
print(f'\t -> File saved in plots folder as "integrated_emission_line.png')
print(f'\t -> Diagnostic plot "moment0_map_used_for_mask is also saved"')

# generate the moment maps
print('Generating the moment maps...')
print('\t Making moment maps (0,1,2)...')
not_masked_moments = generate_moments.NonMaskedMomentMaker(centered_cube,'data/HZ7')
not_masked_moments.generate_standard_moment_maps(start_channel, end_channel)
print('\t Making masked moment maps (0,1,2)...')
masked_moments = generate_moments.MaskedMomentMaker(centered_cube, 'data/HZ7')
masked_moments.generate_standard_moment_maps(start_channel, end_channel)
print('Moment maps created and are in the data folder')

# Create the plots for the moment maps
print()
print('Plotting the moment maps...')
moment_0 = MomentMap('data/HZ7_integrated.fits', order = 0)
moment_1 = MomentMap('data/HZ7_Masked_weighted_coord.fits', order = 1)
moment_2 = MomentMap('data/HZ7_Masked_weighted_dispersion_coord.fits', order = 2)

print('\t plotting 0')
moment_0.plot_moment('plots/moment0.png', 'data/HZ7_integrated.fits', cmap = 'inferno')
print('\t plotting 1')
moment_1.plot_moment('plots/moment1.png', 'data/HZ7_integrated.fits', vmin = -100, vmax = 100, cmap = 'coolwarm')
print('\t plotting 2')
moment_2.plot_moment('plots/moment2.png', 'data/HZ7_integrated.fits', cmap = 'inferno')
print('Moment maps saved in plotting as moment0.png, moment1.png, and moment2.png')

# Create cubes which we can compare to hubble

"""First thing we want to do is align the cubes using GAIA"""
print('Aligning hubble images with GAIA (user interaction required)')
gaia_corrected_file_names = []
for hubble_image in hubble_images:
    gaia_corrected_file_name = hubble_image.replace(
        'raw', 'gaiacorrected').split('.fits')[0] + '_gaia_corrected.fits'

    GaiaAlignment(hubble_image).applySimpleCorrection(gaia_corrected_file_name)
    gaia_corrected_file_names.append(gaia_corrected_file_name)
print('Images aligned are saved in the data/ASCImages/gaiacorrected folder')

"""Next thing we need to do is convolve the images to the radio beam"""
print('Convolving hubble images to Radio Cube beam size (this might take a minute)')
convolved_file_names = []
for gaia_corrected_file_name in gaia_corrected_file_names:
    convolved_file_name = gaia_corrected_file_name.replace(
        'gaiacorrected', 'convolved').split('.fits')[0] + '_convolved.fits'
    convolved_file_names.append(convolved_file_name)
   
    generateConvolvedImages(
        gaia_corrected_file_name, 'data/HZ7_integrated.fits', convolved_file_name, 'TEST_RADIO_CONVOLVED.fits'
    )
print('Convolution done. Files saved in data/ASCImages/convolved')

"""Make Radio Overlay Image on top of hubble image with continuum"""

