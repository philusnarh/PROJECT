#!/usr/bin/env python
"""
healpix_gridding.py
===================
"""
#
import pywcsgrid2 #as pg
import pyfits as pf
import numpy as np
#from numpy import linalg as la
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import healpy as hp
from scipy import interpolate
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
def gridding(mapp, theta, phi,radius, sigma):
    npix = mapp.size 		# checks for number of pixels
    nside = np.sqrt(npix/12.0)	# the nside paramter corresponding to npix
    if nside != np.floor(nside):
       raise RuntimeError('check, pixel number must be : 12*nside**2')
    nside = int(np.floor(nside))
    map2 = mapp.ravel()
    outpt = np.zeros(shape = (theta.shape[0], theta.shape[1]), dtype=np.float64)
    for i in xrange(theta.shape[0]):
    	for j in xrange(phi.shape[0]):    		
    		vec = np.array([np.sin(theta[i,j])*np.cos(phi[i, j]), np.sin(theta[i,j])*np.sin(phi[i,j]), np.cos(theta[i,j])])
    		r  = hp.query_disc(nside, vec, 2*np.pi*radius/360./60.) 
    		v = np.exp(-((theta[i,j] - np.ones(r.shape[0]))**2 + (phi[i,j] - np.array(r, dtype = np.float)/r.sum())**2)\
    		/(2.*sigma**2))
		v = v/v.sum()
	        outpt[i,j] = np.sum(v*map2[r],0)	        	
    return outpt				
#
def run_test():
	option = OptionParser()
  	option.set_description('Adds Kernel size for convolution & FITS filename')
 	option.set_usage('%prog [options] arg1 arg2')
  	option.add_option('-s','--sigma',dest='sigma',default=1.0,type = float,help='Enter the standard deviation for \
  	the Gaussian kernel, Optional')
  	option.add_option('-k','--rad0',dest='radd',type = float, help='Enter r0 in arc_minutes')
  	option.add_option('-r','--rad',dest='rad',type = float, help='Enter ra in radians')
  	option.add_option('-d','--dec',dest='dec',type = float, help='Enter dec in radians')
  	option.add_option('-a','--nra',dest='nra',type = int, help='Enter number of pixels')
  	option.add_option('-m','--deltaRA',dest='deltaRA',default = 0.3, type = float, help='Enter Change in RA in arc_minute')
  	option.add_option('-n','--deltaDEC',dest='deltaDEC',default = 0.3, type = float,help='Enter Change in DEC in arc_minute')
  	option.add_option('-f','--frequency',dest='frequency',default=0,type=int,help='Enter Frequency Value')
  	option.add_option('-i','--stokes',dest='stokes',default=0,type=int,help=' For Stokes Q, U & V  Enter {1, 2, 3} respectively')
  	options,args = option.parse_args(sys.argv[1:])
  	filename = args[0]
	print "Opening FITS File for gridding ..."
	print "WILL TAKE SOME TIME ..."
	#load HEALPix File
	hdu = openFitsFile(filename)	
	if (options.frequency < hdu[0].data.shape[0] and options.stokes < hdu[0].data.shape[1]):
		mapp = 1e-03*hdu[0].data[options.frequency][options.stokes] 
	else:
		raise RuntimeError('Wrong, frequency or stokes valaue exceeds data range')
	x = np.arange(options.rad, options.rad + 2*np.pi/360./60.*options.deltaRA*options.nra, 2*np.pi/360./60.*options.deltaRA)
        y = np.arange(options.dec, options.dec + 2*np.pi/360./60.*options.deltaDEC*options.nra, 2*np.pi/360./60.*options.deltaDEC)
        #x = np.linspace(options.rad, options.dec, num=options.nra, endpoint=True, retstep=False, dtype=None)
        #y = np.linspace(options.rad, options.dec, num=options.nra, endpoint=True, retstep=False, dtype=None)
        ra1, dec1 = np.meshgrid(x, y)
        H = gridding(mapp, ra1, dec1,options.radd, options.sigma)
	print H.shape
	np.save('smap', H)
	#np.save('gridmapQ61', H)
	#print H
	#plt.imshow(H, origin = 'lower')
	#sc = plt.imshow(H, origin = 'lower')
	#plt.colorbar()
	#cb = plt.colorbar(sc)
	#cb.set_label('Jy/pixel')
	#h = np.load('headr.txt.npy')
	ty = pf.open('/home/narh/telescopemodel/kat7__Mueller_00.fits')
	h = ty[0].header
	plt.figure(1)  #
	ax = pywcsgrid2.subplot(111, header=h)
	im = plt.imshow(H, origin="lower", interpolation='lanczos')
	ax.locator_params(axis="y", nbins=5)
	ax.set_xlabel('$Right Ascension\, (J\, 2000)$') 
	ax.set_ylabel('$ Declination\, (J\, 2000)$')
	cbar = plt.colorbar(im,fraction=0.048, pad=0.012)
	cbar.set_label('$Jy\, /\, pixel$', rotation=270)
	ax.add_size_bar(4, 'Stokes_I, 300 MHz', loc=1)
	#pf.writeto('Foreground_Gridding__300MHz_Q.fits',H, h[0:38], clobber=True)
	#pf.writeto('F_Gridding__300MHz_128I.fits',H, h[0:38], clobber=True)
	#plt.xlabel('X')
	#plt.ylabel('Y')
	#plt.colorbar()
	plt.figure(2)  #
	outKSize = 1024
	xx = np.linspace(x.min(),x.max(),outKSize)
        yy = np.linspace(y.min(),y.max(),outKSize)
	kernelFxn = interpolate.RectBivariateSpline(np.sort(x),np.sort(y),H, kx=3,ky=3)
	kernelFxn1 = kernelFxn(xx,yy)
	im = plt.imshow(kernelFxn1,origin = 'lower', interpolation='lanczos')
	#im = plt.imshow(H, origin="lower", interpolation='lanczos')
	ax.locator_params(axis="y", nbins=5)
	ax.set_xlabel('$Right Ascension\, (J\, 2000)$') 
	ax.set_ylabel('$ Declination\, (J\, 2000)$')
	cbar = plt.colorbar(im,fraction=0.048, pad=0.012)
	cbar.set_label('$Jy\, /\, pixel$', rotation=270)
	ax.add_size_bar(4, 'Stokes_I, 300 MHz', loc=1)
	plt.grid()
	plt.show()
	#python regridconvolutionnew.py -s 10 -f 300 -i 2 -k 120 -r 0  -a 512 -d -0.61 -m -8.26 -n -14.4 560_4_786432_K2JanskypPixel.fits
	#python regridconvolutionnew.py -s 10 -f 300 -i 0 -k 120 -r 0  -a 512 -d -0.61 -m 7.26 -n 5.8  560_4_786432_K2JanskypPixel.fits 
	#python regridconvolutionnew.py -s 10 -f 300 -i 1 -k 120 -r 0  -a 512 -d -0.61 -m 9.42 -n -9.42 560_4_786432_K2JanskypPixel.fits 
if __name__ == '__main__':
	run_test()
	
	
