#!/usr/bin/env python
"""
compare-oskar-vla-beams.py
====================

This program compares mueller and oskar beams.

python setup.py install --user
"""
#
#from __future__ import division
import glob
import pyfits 
import pywcsgrid2
import pywcs
#from kapteyn import maputils, tabarray
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
from matplotlib import pylab as plt
import numpy as np
from scipy.stats import expon
import scipy.interpolate as interpolate
import warnings
warnings.filterwarnings("ignore")
import os, sys
import shutil
#from ConfigParser import SafeConfigParser
import ConfigParser
from scipy import signal
#
#
def openFitsFile(filename):
	    """ Opens a FITS Filen using PyFITS
	    Parameters
	    ----------
	    filename:str
		path to file
	    """
	    try:
		hdu = pyfits.open(filename)
		return hdu
	    except:
		print "Error: cannot open %s"%filename
		raise
#
def Mueller_image(fig, header, nrows, ncols, i , Mueller_matrix, figname):
	    #	
	    g2, cax2 = setup_axes(fig, header, nrows, ncols, i)
	    #
	    # Plot M Images
	    #
	    cmap = plt.get_cmap('jet')	
	    images = []
	    for i, ax in enumerate(g2):
		channel = Mueller_matrix[i]			
		im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos")
		if i in [0, 5, 10, 15]:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel),  np.max(channel)], orientation = "horizontal") 
		    # ticks = [0,  np.max(channel)]
		   cax2[i].xaxis.set_ticks_position("bottom")  # ------ just added
		else:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel), np.max(channel)], orientation = "horizontal")
		   cax2[i].xaxis.set_ticks_position("bottom")  # ------ just added
		images.append(im)
	    plt.savefig(figname)
	    return None 
	    
#
def headers(Fits_filename):
	    hd = pyfits.getheader(Fits_filename)
	    return hd
#    
def figur(n1, n2, n3, n4):
    fx = plt.figure(n1, figsize=(n2, n3), dpi= n4)
    return fx
#
def setup_axes(fig, header, nrows, ncols, i):
    gh = pywcsgrid2.GridHelper(wcs=header)
    gh.locator_params(nbins=3)
    g = axes_grid.ImageGrid(fig, 111,
                            nrows_ncols=(nrows, ncols),
                            ngrids=None,
                            direction='row',
                            axes_pad=0.02, add_all=True,
                            share_all=True, aspect=True,
                            label_mode='L', cbar_mode=None,
                            cbar_location='right',
                            axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh))) ##
    # make colorbar
    cax = []
    for iter in range(i):
        #ax = g[iter]
        ax = inset_axes(g[iter],
                     width="45%", # width = 10% of parent_bbox width ---- "5%"
                     height="5%", # height : 50%  ------- "70%"
                     loc=9,
                     #bbox_to_anchor=(1.05, 0., 1, 1), # (1.01, 0, 1, 1), --------(0.025, 0, 1.05, 0.95)
                     #bbox_transform=g[iter].transAxes,
                     #borderpad=0.
                     )
        cax.append(ax)
    return g, cax
    #
# 
def run_compare_beam():
	#
#	In [25]: a = 1008000000

#	In [26]: d = 1000000.

#	In [27]: tn = a +(1024-1)*d  ==> n = (tn - a)/d + 1

	antnum1 = 5
	antnum2 = 6
	#osk_beam = 'oskarr_beam_XY__mueller_4_x_4_beam.%s' %('fits')
	for chan in [0, 100, 230, 340, 450]:    #[0, 9, 19, 29, 39, 49]:	[0, 49, 99, 149, 199, 249]
		vla_beam1 = 'holo_a5/vla_ant%d_chan%s_mueller_im.%s' %(antnum1, str(chan).zfill(4), 'fits')	
		vla_beam2 = 'holo_a6/vla_ant%d_chan%s_mueller_im.%s' %(antnum2, str(chan).zfill(4), 'fits')	
		hdu_1 = openFitsFile('%s' %(vla_beam1)) 
		hdu_2 = openFitsFile('%s' %(vla_beam2))
		beam_difference  = -hdu_2[0].data + hdu_1[0].data
		#	header
		#
		header = headers('%s' %(vla_beam1))
		#print header 
		#header.update('CTYPE1', 'RA---SIN')
		#header.update('CTYPE2', 'DEC--SIN')
		#
		
		fv3 = figur(1, 12, 12, 80)
	    	Mueller_image(fv3, header, 4, 4, 16 , hdu_1[0].data, 'vla-ant%d_chan%s_4_X_4_Images.png' %(antnum1, str(chan).zfill(4)))	    	
	    	fv3 = figur(2, 12, 12, 80)
	    	Mueller_image(fv3, header, 4, 4, 16 , hdu_2[0].data, 'vla-ant%d_chan%s_4_X_4_Images.png' %(antnum2, str(chan).zfill(4)))
	    		    		
	    	fv3 = figur(3, 12, 12, 80)
	    	
	    	Mueller_image(fv3, header, 4, 4, 16 , beam_difference, 'vla_ant%d_ant%d_chan%s_diff_4_X_4_Images.png'
	    		     %(antnum1,antnum2,str(chan).zfill(4)))
      #
	result_directory = 'vla_plots'
	if os.path.isdir("%s" %(result_directory)):
		os.popen('rm -r %s' %result_directory)
	os.makedirs('%s' %(result_directory))
	os.system('mv -f *.%s %s' %('png', result_directory))
	
	print '\n >> Done !!!!!!!!!!'
    	
    	#
if __name__ == '__main__':
	run_compare_beam()
