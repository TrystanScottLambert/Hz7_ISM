from mpl_toolkits.axes_grid1 import ImageGrid
from qubefit.qfutils import qubebeam
import numpy as np 
import pylab as plt
import warnings
warnings.filterwarnings("ignore")

def PlotPlot(Cube):
  # model, data and residual
  Mom0Data = Cube.calculate_moment(moment=0)
  Mom0RMS = Mom0Data.calculate_sigma()
  Mom0Model = Cube.calculate_moment(moment=0, use_model=True)
  Mom0Res = Mom0Data.data - Mom0Model.data

  Mask = Mom0Data.mask_region(value=3*Mom0RMS, applymask=False)
  CubeClip = Cube.mask_region(value=2*np.sqrt(Cube.variance[:, 0, 0]))
  Mom1Data = CubeClip.calculate_moment(moment=1)
  Mom1Data = Mom1Data.mask_region(mask=Mask)
  Mom1Model = Cube.calculate_moment(moment=1, use_model=True)
  Mom1Model = Mom1Model.mask_region(mask=Mask)
  Mom1Res = Mom1Data.data - Mom1Model.data

  # ranges to plot
  vrange = np.array([-3.2 * Mom0RMS, 11 * Mom0RMS])
  levels = np.insert(3 * np.power(np.sqrt(2), np.arange(0, 5)), 0, [-4.2, -3]) * Mom0RMS
  vrange2 = np.array([-220, 220])

  # figure specifics
  #fig = plt.figure(figsize=(7.2, 4.67))
  fig = plt.figure(figsize=(14, 9))
  grid = ImageGrid(fig, (0.1, 0.53, 0.80, 0.45), nrows_ncols=(1, 3), axes_pad=0.15,
                   cbar_mode='single', cbar_location='right', share_all=True)
  labels = ['data', 'model', 'residual']
  for ax, label in zip(grid, labels):
      ax.set_xticks([]); ax.set_yticks([])
      ax.text(0.5, 0.87, label, transform=ax.transAxes, fontsize=20,
              color='k', bbox={'facecolor': 'w', 'alpha': 0.8, 'pad': 2},
              ha='center')

  # the moment-zero images
  ax = grid[0]
  im = ax.imshow(Mom0Data.data, cmap='RdYlBu_r', origin='lower', vmin=vrange[0],
                 vmax=vrange[1])
  ax.contour(Mom0Data.data, levels=levels, linewidths=1, colors='black')
  ax.add_artist(qubebeam(Mom0Data, ax, loc=3, pad=0.3, fill=None, hatch='/////',
                         edgecolor='black'))
  ax = grid[1]
  ax.imshow(Mom0Model.data, cmap='RdYlBu_r', origin='lower', vmin=vrange[0],
            vmax=vrange[1])
  ax.contour(Mom0Model.data, levels=levels, linewidths=1, colors='black')
  ax = grid[2]
  ax.imshow(Mom0Res, cmap='RdYlBu_r', origin='lower',vmin=vrange[0], vmax=vrange[1])
  ax.contour(Mom0Res, levels=levels, linewidths=1, colors='black')
  plt.colorbar(im, cax=grid.cbar_axes[0], ticks=np.arange(-10, 10) * 0.2)

  # the moment-one images
  grid = ImageGrid(fig, (0.1, 0.06, 0.80, 0.45), nrows_ncols=(1, 3),
                   axes_pad=0.15, cbar_mode='single', cbar_location='right',
                   share_all=True)
  labels = ['data', 'model', 'residual']
  for ax, label in zip(grid, labels):
      ax.set_xticks([]); ax.set_yticks([])
      ax.text(0.5, 0.87, label, transform=ax.transAxes, fontsize=20,
              color='k', bbox={'facecolor': 'w', 'alpha': 0.8, 'pad': 2},
              ha='center')

  ax = grid[0]
  im = ax.imshow(Mom1Data.data, cmap='Spectral_r', origin='lower', vmin=vrange2[0],
                 vmax=vrange2[1])
  ax = grid[1]
  ax.imshow(Mom1Model.data, cmap='Spectral_r', origin='lower', vmin=vrange2[0],
            vmax=vrange2[1])
  ax = grid[2]
  ax.imshow(Mom1Res, cmap='Spectral_r', origin='lower',vmin=vrange2[0], vmax=vrange2[1])
  plt.colorbar(im, cax=grid.cbar_axes[0], ticks=np.arange(-10, 10) * 100)

  # some labels
  fig.text(0.5, 0.49, 'Mean velocity (km s$^{-1}$)',
           fontsize=14, ha='center')
  fig.text(0.5, 0.96,
           'Velocity-integrated flux density (Jy km s$^{-1}$ beam$^{-1}$)',
           fontsize=14, ha='center')
  plt.show()