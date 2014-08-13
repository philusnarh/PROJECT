#!/usr/bin/env python
"""
kelvin2stokes.py
===========================
This program converts image data in Kelvin (K) into Stokes in Jansky (Jy)
"""
import healpy as hp
import pyfits as pf
import numpy as np
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
def kelvinConversion(temp):
    """ Converts Kelvin into Stokes by using Rayleigh-Jeans law which attempts 
    to describe the spectral radiance of electromagnetic radiation at all wavele    ngths

    Parameters
    ----------
    temp:image data in Kelvin
    K: Boltzmann constant
    c: speed of light
    v: frequency
    """
    k = 1.3806488e-23 
    c = 299792458
    I = []
    for v in range(560):
        Iv = 2*((v + 1)*1e6)**2 * k * temp[v]/c**2   # Rayleigh-Jeans law
        I.append(Iv)
    return np.array(I)/1e-26
#
def run_tests():
    """Tests conversion"""
    print "Opening FITS File for conversion ..."
    #load FITS File
    hdu = openFitsFile('newhd2fits.fits')
    temp = 1e-6*hdu[0].data.reshape(560, 4, 786432)
    print "Running conversion ..."
    M =  kelvinConversion(temp)
    print "Stokes: in Jy:", M
    print "Creating New FITS File ..."
    hdu1 = pf.PrimaryHDU(M)
    hdu1.writeto('Kel2Stok.fits')
    print "Conversion Don!! ..."
#
if __name__ == '__main__':
   run_tests()
