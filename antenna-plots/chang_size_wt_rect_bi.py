import pyfits as pf
#from pyfits import getdata, getheader
import numpy as np
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import time
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
#def rescale_beam_size(H, outKSize):
#	outKSize = 256
#	xx = np.linspace(x.min(),x.max(),outKSize)
#	yy = np.linspace(y.min(),y.max(),outKSize)
#	kernelFxn = interpolate.RectBivariateSpline(np.sort(x),np.sort(y),H, kx=3,ky=3)
#	kernelFxn1 = kernelFxn(xx,yy)
#	np.save('UU256', kernelFxn1)
#	im = plt.imshow(kernelFxn1,origin = 'lower', interpolation='lanczos')
#
def run_test():
	start = time.time()
        startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Start at %s" % startime
        #
        beam_file = 'jones_amplitude.fits' #'jones_amplitude_ideal.fits'
        hdu = openFitsFile('%s' %(beam_file))
        beam_data = hdu[0].data[0]
        beam_header = hdu[0].header
        x1 = 5./256.*(np.linspace(-256./2., 256./2., 256))
        print len(x1)
        y1 = 5./256.*(np.linspace(-256./2., 256./2., 256))
        print len(y1)
        print beam_data.shape
        outKSize = 128
        xx =  beam_header['CDELT2']*(np.linspace(- outKSize/2.,  outKSize/2.,  outKSize))
	yy =  beam_header['CDELT2']*(np.linspace(- outKSize/2.,  outKSize/2.,  outKSize))
        kernelFxn = interpolate.RectBivariateSpline(np.sort(xx), np.sort(yy),beam_data, kx=3,ky=3)
        print 'Done !!!!!!!!'
        kernelFxn1 = kernelFxn(x1,y1)
        plt.imshow(kernelFxn1)
        plt.colorbar()
        plt.figure(2)
        plt.plot(x1, kernelFxn1.diagonal())
        plt.show()
        #
       
        stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0		
        print "Total run time: %7.2f hours" % elasped_time











if __name__ == '__main__':
	run_test()
