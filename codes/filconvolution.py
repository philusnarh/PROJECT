import pyfits as pf
import numpy as np
#from numpy import linalg as la
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import healpy as hp
#
def openFitsFile(filename):
    """ Opens a FITS Filen using PyFITS
    Parameters
    ----------
    filename:str
        path to file
    """
    try:
        hdu = pf.open(filename)
        return hdu
    except:
        print "Error: cannot open %s"%filename
        raise
#
def gridding(mapp, theta, phi, nest, sigma):
    npix = mapp.size
    nside = np.sqrt(npix/12.0)	# the nside paramter corresponding to npix
    if nside != np.floor(nside):
       raise RuntimeError('check, pixel number must be : 12*nside**2')
    nside = int(np.floor(nside))
    x = np.arange(theta, phi, 0.1)
    y = np.arange(theta, phi, 0.1)
    ra, dec = np.meshgrid(x, y)
    r = hp.pixelfunc.get_all_neighbours(nside, ra, dec)
    '''if nest:
       r = pxl._get_interpol_nest(nside, theta, phi)
    else:
       r = pxl._get_interpol_ring(nside, theta, phi)		#searches for nea	
    pl = map2[np.array(r[0:4])]	#nearest pixel
    wgt = np.array(r[4:8])	#weight values
    del r'''
    p = np.array(r[:])
    
    map2 = mapp.ravel()
    outpt = np.zeros(shape = (ra.shape[0], dec.shape[0]), dtype=np.float64)
    for i in xrange(ra.shape[0]):
    	for j in xrange(dec.shape[0]):
    		v = np.exp(-((theta - np.ones(8))**2 + (phi - np.array(p[:,i,j], dtype= np.float)/p[:,i,j].sum())**2)/(2.*sigma*2))
		#v1.append(v)
		v = v/v.sum()
	        outpt[i,j] = np.sum(v*map2[p[:,i,j]],0)
	               # return outpt
    	
    return outpt				#np.sum(map2[p]*w,0)
#

#
def run_test():
	option = OptionParser()
  	option.set_description('Adds Kernel size for convolution & FITS filename')
 	option.set_usage('%prog [options] arg1 arg2')
  	option.add_option('-r','--theta',dest='theta',default= 0.0,type=float,help='Enter Right Ascension')
  	option.add_option('-d','--phi',dest='phi',default= 1.0,type=float,help='Enter Declination')
  	option.add_option('-s','--sigma',dest='sigma',default=1.0,type=float,help='Enter the standard deviation for \
  	the Gaussian kernel, Optional')
  	option.add_option('-n','--nest',dest='nest',default= False,help='Enter nest type')
  	option.add_option('-f','--frequency',dest='frequency',default=0,type=int,help='Enter Frequency Value')
  	option.add_option('-i','--stokes',dest='stokes',default=0,type=int,help=' For Stokes Q, U & V  Enter {1, 2, 3} respectively')
  	options,args = option.parse_args(sys.argv[1:])
  	filename = args[0]
	print "Opening FITS File for gridding ..."
	print "WILL TAKE SOME TIME ..."
	#load HEALPix File
	hdu = openFitsFile(filename)	
	if (options.frequency < hdu[0].data.shape[0] and options.stokes < hdu[0].data.shape[1]):
		mapp = hdu[0].data[options.frequency][options.stokes]
	else:
		raise RuntimeError('Wrong, frequency or stokes valaue exceeds data range')
	H = gridding(mapp, options.theta,options.phi, options.nest, options.sigma)
		#('/home/narh/theo/comp_sc/emma_n_theo/plott/mkelvintheodata.fits')
	#hdu = openFitsFile('/home/narh/theo/comp_sc/emma_n_theo/plott/no_50_200MHz_256pix_K2Jansky.fits')
	#hdu = openFitsFile('/home/narh/theo/comp_sc/emma_n_theo/plott/noiseCubes/mKel2D.fits')
	#inputdata = hdu[0].data[0][1]
	#H = gridding(inputdata, 0, 2, nest = False, sigma = 1.0)
	print H
	plt.imshow(H, origin = 'lower')
	plt.grid()
	plt.show()
	
if __name__ == '__main__':
	run_test()
	
	
