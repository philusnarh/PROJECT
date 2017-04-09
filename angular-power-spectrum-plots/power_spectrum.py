#!/usr/bin/env python
#
import pyfits as pf
import numpy as np
from pyfits import getdata, getheader
#from numpy import linalg as la
from multiprocessing import Pool
import multiprocessing
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import healpy as hp
from scipy import interpolate
import time
import copy
import cPickle
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
def run_test():
    hdu = openFitsFile('560_4_786432_K2JanskypPixel.fits')
    mapp = hdu[0].data
    freq, stokes, npixels = mapp.shape
    #
    fig1, axes1 = plt.subplots(nrows=4, ncols=2) 
    #fig1.tight_layout()
    iter = 0
    for iter1 in xrange(4):      
        for i in range(2):
           for j in range(2):
               lmax = np.arange(len(hp.sphtfunc.anafast(mapp[iter, 0, :])))
               axes1[i, j].plot(lmax, lmax*(lmax +1.)*hp.sphtfunc.anafast(mapp[iter, 0, :])/2.0*np.pi, label = 'Stokes_I')
               axes1[i, j].plot(lmax, lmax*(lmax +1.)*hp.sphtfunc.anafast(mapp[iter, 1, :])/2.0*np.pi, label = 'Stokes_Q')
               axes1[i, j].plot(lmax, lmax*(lmax +1.)*hp.sphtfunc.anafast(mapp[iter, 2, :])/2.0*np.pi, label = 'Stokes_U')
               axes1[i, j].legend(loc = 'upper left') 
               axes1[i, j].set_xlabel('$lmax$')
	       axes1[i, j].set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
               axes1[i, j].set_title('Power Spectrum at %d MHz '%(iter + 1))
               fig1.tight_layout()
               plt.show()
               fig1.savefig('spectrum%d.png'%iter1) #dpi=100
               
               iter = iter + 1
            
#
'''
In [147]: %history
import numpy as np
from matplotlib import pylab as plt
plt.ion()
import healpy as hp
import pyfits
ls
hdu = pyfits.open('560_4_786432_K2JanskypPixel.fits')
dat = hdu[0].data
dat.shape
x1 = dat[0, 0:3, ...]
u = hp.sphtfunc.anafast(x1)
lmax = np.arange(u)
lmax = np.arange(len(u))
u.shape
lmax = np.arange(len(u[0]))
plt.plot(lmax,lmax*(lmax +1.)*u[0]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[0]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[0]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[1]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[2]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[3]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[4]/2.0*np.pi, label = 'Stokes_I')
plt.plot(lmax,lmax*(lmax +1.)*u[5]/2.0*np.pi, label = 'Stokes_I')
x1 = dat[0, 0:3, ...]
x1 = dat[0, 0, ...]
u = hp.sphtfunc.anafast(x1)
lmax = np.arange(len(u))
plt.plot(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
x1 = dat[0, 1, ...]
u = hp.sphtfunc.anafast(x1)
lmax = np.arange(len(u))
plt.plot(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
x1 = dat[0, 2, ...]
u = hp.sphtfunc.anafast(x1)
plt.plot(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
ax1 = plt.subplot2grid((3,3), (0,0))
u = hp.sphtfunc.anafast(dat[0, 0, ...])
lmax = np.arange(len(u))
ax1.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[0, 1, ...])
ax1.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[0, 2, ...])
ax1.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax1.legend(loc = 'best')
ax1.set_xlabel('$lmax$')
ax1.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax1.set_title('Power Spectrum at 1 MHz ')
ax2 = plt.subplot2grid((3,3), (0,1))
u = hp.sphtfunc.anafast(dat[1, 0, ...])
ax2.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[1, 1, ...])
ax2.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[1, 2, ...])
ax2.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax2.legend(loc = 'best')
ax2.set_title('Power Spectrum at 2 MHz ')
ax2.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax2.set_xlabel('$lmax$')
ax3 = plt.subplot2grid((3,3), (0,2))
u = hp.sphtfunc.anafast(dat[2, 0, ...])
ax3.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[2, 1, ...])
ax3.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[2, 2, ...])
ax3.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax3.set_xlabel('$lmax$')
ax3.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax3.set_title('Power Spectrum at 3 MHz ')
ax3.legend(loc = 'best')
ax4 = plt.subplot2grid((3,3), (1,0))
u = hp.sphtfunc.anafast(dat[3, 0, ...])
ax4.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[3, 1, ...])
ax4.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[3, 2, ...])
ax4.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax4.legend(loc = 'best')
ax4.set_title('Power Spectrum at 4 MHz ')
ax4.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax4.set_xlabel('$lmax$')
ax5 = plt.subplot2grid((3,3), (1,1))
u = hp.sphtfunc.anafast(dat[4, 0, ...])
ax5.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[4, 1, ...])
ax5.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[4, 2, ...])
ax5.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax5.set_xlabel('$lmax$')
ax5.legend(loc = 'best')
ax5.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax5.set_title('Power Spectrum at 5 MHz ')
ax6 = plt.subplot2grid((3,3), (1,2))
u = hp.sphtfunc.anafast(dat[5, 0, ...])
ax6.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[5, 1, ...])
ax6.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[5, 2, ...])
ax6.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax6.legend(loc = 'best')
ax6.set_title('Power Spectrum at 6 MHz ')
ax6.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax6.set_xlabel('$lmax$')
ax7 = plt.subplot2grid((3,3), (2,0))
u = hp.sphtfunc.anafast(dat[6, 0, ...])
ax7.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[6, 1, ...])
ax7.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[6, 2, ...])
ax7.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax7.legend(loc = 'best')
ax7.set_xlabel('$lmax$')
ax7.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax7.set_title('Power Spectrum at 7 MHz ')
ax8 = plt.subplot2grid((3,3), (2,1))
u = hp.sphtfunc.anafast(dat[7, 0, ...])
ax8.set_title('Power Spectrum at 8 MHz ')
ax8.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[7, 1, ...])
ax8.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[7, 2, ...])
ax8.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax8.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
ax8.legend(loc = 'best')
ax8.set_xlabel('$lmax$')
ax9 = plt.subplot2grid((3,3), (2,2))
u = hp.sphtfunc.anafast(dat[8, 0, ...])
ax9.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_I')
u = hp.sphtfunc.anafast(dat[8, 1, ...])
ax9.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_Q')
u = hp.sphtfunc.anafast(dat[8, 2, ...])
ax9.plot(lmax, lmax*(lmax +1.)*u/2.0*np.pi, label = 'Stokes_U')
ax9.set_xlabel('$lmax$')
ax9.set_ylabel('$lmax\, (lmax+1)\, cl/2\pi$ (Jy/pixel)')
plt.plot??
ax9.set_xticks??
ax9.set_yticks??
ax9.format_ydata??
ax9.yaxis.set_major_formatter(FormatStrFormatter('%0.0e'))
import matplotlib.ticker as mtick
ax9.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
ax9.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax9.set_xlabel('$lmax$')
ax9.legend(loc = 'best')
ax9.set_title('Power Spectrum at 9 MHz ')
ax9.set_yscale??
plt.subplots_adjust(hspace=.5)
%history


'''
if __name__ == '__main__':
	run_test()

