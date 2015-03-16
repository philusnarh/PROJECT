#!/usr/bin/env python
"""
hdf2fits.py
=========
Converts HDF to FITS format
"""
__author__ = "T. Ansah - Narh"
__version__ = "1.0"
import h5py
import pyfits
import numpy as np
from pyfits import getdata, getheader
print " ------------------------------------------------------"
print "\n >> Please enter HDF5 file name to convert (including extension):"
print " ------------------------------------------------------"
filename=raw_input() 	#Input .hdf5 file
print " >> Opening file and reading data"
print " ------------------------------------------------------"
hdf_data =h5py.File(filename,'r')
print " >> Extracting data from HDF file"
print " ------------------------------------------------------"
extract_data = hdf_data.values()[:] 
print " >> Converting HDF data into arrays"
print " ------------------------------------------------------"
a_ray = np.array(extract_data)
###############################################################
print " >> Converting HDF data into FITS"
print " ------------------------------------------------------"
hdu = pyfits.PrimaryHDU(a_ray)
hdu.writeto('newhd2fits.fits')
##############################################################
print " ------------------------------------------------------"
print "\n >> Please checking & fixing headers:"
print " ------------------------------------------------------"
filename1='newhd2fits.fits'  #Input .fits file

print " >> Opening FITS file and reading header"
print " ------------------------------------------------------"

hdulist = pyfits.open(filename1)
in_header=hdulist[0].header
in_header.update('BITPIX', -64)

print " >> BITPIX Fixed\n >> Now checking for all required fits header keywords"
print " ------------------------------------------------------"

try:
    in_header['COMMENT']
except KeyError:
    print "\t No COMMENT Found, adding entry for this"
    in_header.update('COMMENT', 'Standard WCS reduction:')
try:
    in_header['BZERO']
except KeyError:
    print "\t No BZERO Found, adding entry for this"
    in_header.update('BZERO',0.0,'PhysValue = BZERO + BSCALE * ArrayValue')
try:
    in_header['BSCALE']
except KeyError:
    print "\t No BSCALE Found, adding entry for this"
    in_header.update('BSCALE',1.0,'PhysValue = BZERO + BSCALE * ArrayValue')

try:
    in_header['CRVAL1']
except KeyError:
    print "\t No CRVAL1 Found, adding entry for this"
    in_header.update('CRVAL1',11.88798983,'WCS Ref value (RA in decimal degrees)')
try:
    in_header['CRVAL2']
except KeyError:
    print "\t No 'CRVAL2' Found, adding entry for this"
    in_header.update('CRVAL2',-25.33718377,'WCS Ref value (DEC in decimal degrees)')
try:
    in_header['EPOCH']
except KeyError:
    print "\t No 'EPOCH' Found, adding entry for this"
    in_header.update('EPOCH',1999.55329832624,'Epoch in Julain Year at start of observation')
try:
    in_header['CRPIX1']
except KeyError:
    print "\t No 'CRPIX1' Found, adding entry for this"
    in_header.update('CRPIX1',256.5,'WCS Coordinate reference pixel')
try:
    in_header['CRPIX2']
except KeyError:
    print "\t No 'CRPIX2' Found, adding entry for this"
    in_header.update('CRPIX2',356.5,'WCS Coordinate reference pixel')
try:
    in_header['CD1_1']
except KeyError:
    print "\t No 'CD1_1' Found, adding entry for this"
    in_header.update('CD1_1',-0.000277777783923599,'WCS Coordinate scale matrix')
try:
    in_header['CD1_2']
except KeyError:
    print "\t No 'CD1_2' Found, adding entry for this"
    in_header.update('CD1_2',1.78947724260645E-08, 'WCS Coordinate scale matrix')
try:
    in_header['CD2_1']
except KeyError:
    print "\t No CD2_1 Found, adding fake entry for this"
    in_header.update('CD2_1',1.78947724260645E-08,'WCS Coordinate scale matrix')
try:
    in_header['CD2_2']
except KeyError:
    print "\t No CD2_2 Found, adding fake entry for this"
    in_header.update('CD2_2',0.000277777783923599,'WCS Coordinate scale matrix')
try:
    in_header['CD3_1']
except KeyError:
    print "\t No 'CD3_1' Found, adding entry for this"
    in_header.update('CD3_1',1.78947724260645E-08, 'WCS Coordinate scale matrix')
try:
    in_header['CD3_2']
except KeyError:
    print "\t No 'CD3_2' Found, adding entry for this"
    in_header.update('CD3_2',1.78947724260645E-08, 'WCS Coordinate scale matrix')
try:
    in_header['CTYPE1']
except KeyError:
    print "\t No CTYPE1 Found, adding  entry for this"
    in_header.update('CTYPE1','RA---SIN', 'WCS Coordinate type') 
try:
    in_header['CTYPE2']
except KeyError:
    print "\t No CTYPE2 Found, adding  entry for this"
    in_header.update('CTYPE2','DEC---SIN', 'WCS Coordinate type')    
try:
    in_header['EQUINOX']
except KeyError:
    print "\t No EQUINOX Found, adding entry for this"
    in_header.update('EQUINOX',2000.0,' Equinox')
try:
    in_header['CDELT1']
except KeyError:
    print "\t No CDELT1 Found, adding entry for this"
    in_header.update('CDELT1',0.0002777777845,'WCS Coordinate scale matrix')
try:
    in_header['CDELT2']
except KeyError:
    print "\t No CDELT2 Found, adding entry for this"
    in_header.update('CDELT2',0.0002777777845,'WCS Coordinate scale matrix')
try:
    in_header['CTYPE3']
except KeyError:
    print "\t No CTYPE3 Found, adding fake entry for this"
    in_header.update('CTYPE3','STOKES', 'axes 3 is the spectra')
try:
    in_header['CNPIX1']
except KeyError:
    print "\t No CNPIX1 Found, adding fake entry for this"
    in_header.update('CNPIX1','0','New CNPIX1')
try:
    in_header['CNPIX2']
except KeyError:
    print "\t No CNPIX2 Found, adding fake entry for this"
    in_header.update('CNPIX2','0','New CNPIX2')
try:
    in_header['RADESYS']
except KeyError:
    print "\t No RADESYS Found, adding fake entry for this"
    in_header.update('RADESYS','FK5', 'Coordinate system')
try:
    in_header['MJD-OBS']
except KeyError:
    print "\t No MJD-OBS Found, adding fake entry for this"
    in_header.update('MJD-OBS',51381.3422136574,'modified Julian date at start of observation')
try:
    in_header['CRPIX3']
except KeyError:
    print "\t No CRPIX3 Found, adding fake entry for this"
    in_header.update('CRPIX3',1.0)
###############################
try:
    in_header['CUNIT1']
except KeyError:
    print "\t No CUNIT1 Found, adding entry for this"
    in_header.update('CUNIT1','DEG', 'RA coordinate')

try:
    in_header['CUNIT2']
except KeyError:
    print "\t No CUNIT2 Found, adding entry for this"
try:
    in_header['CDELT3']
except KeyError:
    print "\t No CDELT3 Found, adding entry for this"
    in_header.update('CDELT3', 1.0)
try:
    in_header['CRPIX1']
except KeyError:
    print "\t No CRPIX1 Found, adding entry for this"
    in_header.update('CRPIX1', 256.5)
try:
    in_header['CRPIX2']
except KeyError:
    print "\t No CRPIX2 Found, adding entry for this"
###################################################

print " >> All missing required fits header keywords added\n >> except NAXISn. I'll check in a second."
print " ------------------------------------------------------"
    
data=getdata(filename1) #Read in fits data
fl_data=np.float64(data) #convert in array to float64
hdu_data=pyfits.PrimaryHDU(fl_data) #create FITS type data 

print " >> Data array Fixed"
print " ------------------------------------------------------"

naxisvals=fl_data.shape
print naxisvals[0]
print naxisvals[1]

try:
    in_header['NAXIS1']
except KeyError:
    print "\t No NAXIS1 Found, adding entry for this"
    in_header.update('NAXIS1',naxisvals[0])

try:
    in_header['NAXIS2']
except KeyError:
    print "\t No NAXIS2 Found, adding entry for this"
    in_header.update('NAXIS2',naxisvals[0])
try:
    in_header['NAXIS3']
except KeyError:
    print "\t No NAXIS3 Found, adding entry for this"
    in_header.update('NAXIS3',naxisvals[0])

print " >> NAXIS values now checked."
print " ------------------------------------------------------"
hdulist.close()
pyfits.writeto('fixed_'+filename1,fl_data,in_header) #Write new fits file


print " >> Done!\n\n >> The corrected fits file is named fixed_"+filename1
print " ------------------------------------------------------"
