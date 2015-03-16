import matplotlib.pyplot as plt
import pywcsgrid2
import pywcs

import mpl_toolkits.axes_grid1.axes_grid as axes_grid
#from mpl_toolkits.axes_grid.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import pyfits
#

def setup_axes(fig, header):

    gh = pywcsgrid2.GridHelper(wcs=header)
    gh.locator_params(nbins=3)

    g = axes_grid.ImageGrid(fig, 111,
                            nrows_ncols=(4, 4),
                            ngrids=None,
                            direction='row',
                            axes_pad=0.02, add_all=True,
                            share_all=True, aspect=True,
                            label_mode='L', cbar_mode=None,
                            cbar_location='right',
                            axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh))) ##

    # make colorbar
    ax = g[-1]
    cax = inset_axes(ax,
                     width="8%", # width = 10% of parent_bbox width
                     height="100%", # height : 50%
                     loc=3,
                     bbox_to_anchor=(1.01, 0, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0.
                     )
    return g, cax
#
fits_cube = pyfits.open("kat7_my_16_List.fits")	#ic443_co.ambient.fits
header = fits_cube[0].header
#
fig = plt.figure(1, figsize=(12, 12), dpi=70)
g, cax = setup_axes(fig, header)
#
# draw images
#
cmap = plt.get_cmap('jet')	#plt.cm.gray_r
#import matplotlib.colors as mcolors
#norm = mcolors.Normalize()
images = []
#start_channel = i*nxy+dxy
for i, ax in enumerate(g):
    #channel_number = start_channel + i
    channel = fits_cube[0].data[i]	#channel_number
    im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos") 	#norm=norm
    #cb = plt.colorbar(im, cax=ax)
    images.append(im)

# make colorbar
cb = plt.colorbar(im, cax=cax)
cb.set_label("Jy / Beam")
#cb.set_ticks([0, 1, 2, 3])

# adjust norm
'''norm.vmin = -1
norm.vmax = 1
for im in images:
    im.changed()'''

plt.show()    
