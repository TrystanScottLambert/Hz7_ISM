import numpy as np 
import pylab as plt
from mpl_toolkits.axes_grid1 import ImageGrid

def PlotChannelMaps(Cube):
    fig = plt.figure(figsize=(14, 11))
    grid = ImageGrid(fig, (0.06, 0.09, 0.40, 0.94), nrows_ncols=(4, 3),
                     axes_pad=0.0, cbar_mode='single', cbar_location='top',
                     share_all=True)
    CubeRMS = np.sqrt(Cube.variance[:, 0, 0]) * 1E3
    MedRMS = np.median(CubeRMS)
    vrange = [-3 * MedRMS, 11 * MedRMS]
    vel = Cube.get_velocity()
    for idx, channel in enumerate(np.arange(12)):
        ax = grid[idx]
        ax.set_xticks([]); ax.set_yticks([])
        im = ax.imshow(Cube.data[channel, :, :] * 1E3, cmap='RdYlBu_r',
                       origin='lower', vmin=vrange[0], vmax=vrange[1])
        levels = (np.insert(2 * np.power(np.sqrt(2), np.arange(0, 5)), 0, -3) *
                  CubeRMS[channel])
        con = ax.contour(Cube.model[channel, :, :] * 1E3, levels=levels,
                         linewidths=1, colors='black')
        velocity = str(round(vel[channel])) + ' km s$^{-1}$'
        ax.text(0.5, 0.85, velocity, transform=ax.transAxes, fontsize=12,
                color='black', bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 2},
                ha='center')
    cb = plt.colorbar(im, cax=grid.cbar_axes[0], orientation='horizontal', pad=0.1)
    cb.ax.tick_params(axis='x',direction='in',labeltop=True, labelbottom=False)

    grid = ImageGrid(fig, (0.50, 0.09, 0.40, 0.94), nrows_ncols=(4, 3),
                     axes_pad=0.0, cbar_mode='single', cbar_location='top',
                     share_all=True)
    for idx, channel in enumerate(np.arange(12)):
        ax = grid[idx]
        ax.set_xticks([]); ax.set_yticks([])
        im = ax.imshow(Cube.data[channel, :, :] * 1E3, cmap='RdYlBu_r',
                       origin='lower', vmin=vrange[0], vmax=vrange[1])
        levels = (np.insert(2 * np.power(np.sqrt(2), np.arange(0, 5)), 0, -2) *
                  CubeRMS[channel])
        con = ax.contour((Cube.data[channel, :, :] - Cube.model[channel, :, :]) * 1E3,
                         levels=levels, linewidths=1, colors='black')
        velocity = str(round(vel[channel])) + ' km s$^{-1}$'
        ax.text(0.5, 0.85, velocity, transform=ax.transAxes, fontsize=12,
                color='black', bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 2},
                ha='center')
    cb = plt.colorbar(im, cax=grid.cbar_axes[0], orientation='horizontal', pad=0.1)
    cb.ax.tick_params(axis='x',direction='in',labeltop=True, labelbottom=False)

    fig.text(0.5, 0.96, 'Flux density (mJy beam$^{-1}$)', fontsize=14, ha='center')
    plt.show()