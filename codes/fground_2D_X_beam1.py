from ap_array_modelclass1 import *
#from __future__ import division
import glob
import pyfits 
import pywcsgrid2
import pywcs
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
from matplotlib import pylab as plt
import numpy as np
from scipy.stats import expon
import scipy.interpolate as interpolate
import os, sys
import shutil
#from ConfigParser import SafeConfigParser
import ConfigParser
import copy
import time
#
def foregroundGrid_x_beam(healpix_grid_fits, primary_beam_fits):
    hdu_grid = openFitsFile(healpix_grid_fits)
    stokes_I = hdu_grid[0].data[0]
    stokes_Q = hdu_grid[0].data[1]
    stokes_U = hdu_grid[0].data[2]
    stokes_V = hdu_grid[0].data[3]
    hdu_beam = openFitsFile(primary_beam_fits)
    primary_beam = hdu_beam[0].data
    x_new, y_new = primary_beam[0].shape
    measured_map= np.zeros(shape = (x_new, y_new), dtype=np.float64)
    header = headers(primary_beam_fits)
    fgrid_x_beam = []
    temp = []
    for iter in xrange(16):
    	if iter in range(16)[0::4]:
    	   for k1 in xrange(x_new):
    	       for k2 in xrange(y_new):
    	           measured_map[k1, k2] = np.sum(stokes_I[k1:x_new + k1, k2:y_new + k2]*primary_beam[iter])    	   	
	elif iter in range(16)[1::4]:
	     for k1 in xrange(x_new):
	         for k2 in xrange(y_new):
	    	     measured_map[k1, k2] = np.sum(stokes_Q[k1:x_new + k1, k2:y_new + k2]*primary_beam[iter])	     
	elif iter in range(16)[2::4]:
	     for k1 in xrange(x_new):
	    	 for k2 in xrange(y_new):
	    	     measured_map[k1, k2] = np.sum(stokes_U[k1:x_new + k1, k2:y_new + k2]*primary_beam[iter])	     
	else:
	     for k1 in xrange(x_new):
	         for k2 in xrange(y_new):
	    	     measured_map[k1, k2] = np.sum(stokes_V[k1:x_new + k1, k2:y_new + k2]*primary_beam[iter])	     
	fgrid_x_beam.append(measured_map)
	temp.append(copy.deepcopy(measured_map))
    return temp, header
#
def run_test():
    start = time.time()
    startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print "Start at %s" % startime
    '''if os.path.isdir("../fgridxbeamfolder"):
	   os.popen('rm -rf ../fgridxbeamfolder')
    os.makedirs('../fgridxbeamfolder')
    os.chdir('../fgridxbeamfolder')'''
    # 
    # 	Multiply Foreground Grid by True Primary Beam
    #
    tru_image_data, h = foregroundGrid_x_beam('smap00.fits', 'mueller_4_x_4_beam.fits')
    #pyfits.writeto('for_grid_x_trubeam.fits', tru_image_data, h, clobber=True) +++++++++++++
    fu1 = figur(1, 12, 12, 80)
    Mueller_image(fu1, h[0:38], 4, 4, 16 , tru_image_data ,'f_x_trubeam00.png')
    # 
    # 	Multiply Foreground Grid by Distorted Primary Beam
    #
    dis_image_data, h = foregroundGrid_x_beam('smap00.fits', 'dis_mueller_4_x_4_beam.fits')
    #pyfits.writeto('for_grid_x_dis_beam.fits', dis_image_data, h, clobber=True) +++++++++++
    fu2 = figur(2, 12, 12, 80)
    Mueller_image(fu2, h[0:38], 4, 4, 16 ,dis_image_data ,'f_x_disbeam00.png')
    #
    # 	Computes the Difference 
    #
    map_difference = []
    for iter1 in xrange(16):
        difference = dis_image_data[iter1] - tru_image_data[iter1]
        map_difference.append(difference)
    map_difference = np.array(map_difference)
    #print map_difference
    #pyfits.writeto('multiplicative_diff.fits', map_difference , h, clobber=True) +++++++++++++++++
    fu3 = figur(3, 12, 12, 80)
    Mueller_image(fu3, h, 4, 4, 16 , map_difference ,'f_x_beam00diff.png')
    # 
    # 	Computes the Measured map & the Difference 
    #
    mapp = [] 
    '''hdu_grd = openFitsFile('foreground_mapp12.fits')
    m_I = hdu_grd[0].data[0]
    m_Q = hdu_grd[0].data[1]
    m_U = hdu_grd[0].data[2]
    m_V = hdu_grd[0].data[3]
    mapp.append(m_I)
    mapp.append(m_Q)
    mapp.append(m_U)
    mapp.append(m_V)'''
    
    #
    # initializing
    tru_I = 0
    tru_Q = 0
    tru_U = 0
    tru_V = 0
    #measured_map= np.zeros(shape = (x_new, y_new), dtype=np.float64)
    for iter2 in range(4):
        tru_I += tru_image_data[iter2]
    #mapp.append(tru_I)
    for iter2 in range(4, 8):
        tru_Q += tru_image_data[iter2]
    #mapp.append(tru_Q)
    for iter2 in range(8, 12):
        tru_U += tru_image_data[iter2]
    #mapp.append(tru_U)
    for iter2 in range(12, 16):
        tru_V += tru_image_data[iter2]
    mapp.append(tru_I)
    mapp.append(tru_Q)
    mapp.append(tru_U)
    mapp.append(tru_V)
     
    #
    # initializing
    dis_I = 0
    dis_Q = 0
    dis_U = 0
    dis_V = 0
    for iter2 in range(4):
        dis_I += dis_image_data[iter2]
    #mapp.append(dis_I)
    for iter2 in range(4, 8):
        dis_Q += dis_image_data[iter2]
    #mapp.append(dis_Q)
    for iter2 in range(8, 12):
        dis_U += dis_image_data[iter2]
    #mapp.append(dis_U)
    for iter2 in range(12, 16):
        dis_V += dis_image_data[iter2]
    mapp.append(dis_I)
    mapp.append(dis_Q)
    mapp.append(dis_U)
    mapp.append(dis_V)  
    #
    dif_I =  dis_I - tru_I
    dif_Q =  dis_Q - tru_Q
    dif_U =  dis_U - tru_U
    dif_V =  dis_V - tru_V
    mapp.append(dif_I)
    mapp.append(dif_Q)
    mapp.append(dif_U)
    mapp.append(dif_V)
    mapp = np.array(mapp)
    print 'mapp', mapp.shape
    #
    #	pyfits.writeto('map_diff.fits', np.array(mapp) , h, clobber=True)
    #
    fu4 = figur(4, 12, 12, 80)
    Mueller_image(fu4, h, 3, 4, 12 , mapp ,'f_x_beam00msrdmapp.png')
    #pyfits.writeto('map_diff.fits', mapp , h, clobber=True) +++++++++++++++++++++
    ###############################################################################################
    stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print "Stop at %s" % stoptime 
    end = time.time()
    elasped_time = (end - start)/3600.0
    print "Total run time: %7.2f hours" % elasped_time
    print 'Done !!!!!!!!!!!'
if __name__ == '__main__':
    run_test()
