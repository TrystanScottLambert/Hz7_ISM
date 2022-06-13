####################################################
#
# Simple script to remove first couple channels
#
####################################################
import casa 
import os

# Change these as required
original_file = 'data/Jorge_raw_HZ7/HZ7_COMB_CUBE_15_JvM.fits'
start_channel = 21  # use casa viewer to determine these numbers 
end_channel = 143
output_file = 'data/Jorge_cut_HZ7/HZ7_COMB_CUBE_15_JvM_cut.fits'
channels = str(start_channel)+'~'+str(end_channel)

#importing fits image
casa.importfits(fitsimage=original_file, imagename='T11.im', overwrite=True)

#cut cube
casa.imsubimage(imagename="T11.im", outfile="T11_cut.im", chans=channels)

casa.imhead('T11_cut.im', mode='summary')

#exporting fits file and cleaning up
casa.exportfits(imagename='T11_cut.im', fitsimage=output_file, overwrite=True)

os.system('rm -rf T11_cut.im')
os.system('rm -rf T11.im')