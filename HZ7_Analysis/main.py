"""Main pipeline for analysis of HZ7"""

import center_cube_on_CII
import calc_channels
import warnings
import pylab as plt

warnings.filterwarnings("ignore")  # mainly coming from spectral_cube


####### INPUTS #######
combined_cube = 'data/HZ7_Combined.fits'
centered_cube = 'data/HZ7_Centered.fits'
central_channel = 64
centeral_pixel_position = (156, 140)
######################

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

