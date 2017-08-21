#!/usr/bin/env python
"""
convert_npy2fits.py
===================
	This script converts antenna beams of the form (channel,2,2,Nside,Nside)
	from .npy to .fits format
	
	*** TO RUN SCRIPT: python convert_npy2fits.py -f config.ini
"""

from astropy.io import fits as pyfits
from argparse import ArgumentParser
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import os, sys
import ConfigParser
import time
import datetime
import ephem



def azel2radec(az, el, lat, lon, alt, datetime):
	""" This function converts azimuth, elevation to 
		right ascension (RA) and declination (DEC)
		-------------------------------------------
		parameters::
		az: float, azimuth in degrees
		el: float, elevation in degrees
		lat: float, latitute in degrees
		lon: float, longitude in degrees
		alt: float, altitute in metres
		datetime: date & time, eg. "2017/02/23 03:28:04.182"
		return::
		RA & DEC in degrees
		"""
	

	meerkat_observer = ephem.Observer()
	meerkat_observer.lon = lon  
	meerkat_observer.lat = lat  
	meerkat_observer.elevation = alt
	meerkat_observer.date = datetime
	meerkat_observer.epoch = ephem.J2000
	
	print '\n >> RA & DEC :', meerkat_observer.radec_of(float(ephem.degrees(az)),
				  			    float(ephem.degrees(el)))

	return meerkat_observer.radec_of(float(ephem.degrees(az)), float(ephem.degrees(el)))



def create_FITSheader(CRVAL1, CRVAL2, CRVAL3, CRPIX1, CDELT1, 
		      CRPIX3 = None, 
		      CDELT3 = None, 
		      CUNIT1 = 'DEG', 
		      CUNIT3 = 'Hz ',
		      CTYPE1 = 'RA---SIN',
		      CTYPE2 = 'DEC--SIN', 
		      CTYPE3 = 'FREQ ' ):

	"""parameters::
		CRVAL1: float, RA in degrees
		CRVAL2: float, DEC in degrees
		CRVAL3: float, start frequency
		CRPIX1: float, reference pixel for x-axis
		CRPIX2: float, reference pixel for y-axis
		CDELT1: float, pixel size in degrees for x-axis
		CRPIX3: float, reference pixel step
		CDELT3: float, step frequency/ channel width
		return::
		FITs headers
		"""


		      
 	# Create header
   	hdr = pyfits.Header()

	hdr['CRVAL1'] = CRVAL1                                                  
	hdr['CRVAL2'] = CRVAL2	
	hdr['CRVAL3'] = CRVAL3                                                 
	hdr['CRPIX1'] = CRPIX1
	hdr['CRPIX2'] = CRPIX1
	if CRPIX3 is not None:
		hdr['CRPIX3'] = CRPIX3                                                               
                                                
	hdr['CDELT1']  =  -CDELT1                                                  
	hdr['CDELT2']  =  CDELT1 
	
	if CDELT3 is not None:                                                 
		hdr['CDELT3']  =  CDELT3 
		                                                  
	hdr['CUNIT1']  = CUNIT1                                                           
	hdr['CUNIT2']  = CUNIT1                                                             
	hdr['CUNIT3']  = CUNIT3                                                          
	hdr['CTYPE1']  = CTYPE1                                                            
	hdr['CTYPE2']  = CTYPE2                                                             
	hdr['CTYPE3']  = CTYPE3 
	
	return hdr 


def run_header_test():

	#	start-time
	start = time.time()
	startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Start at %s" % startime


	for i, arg in enumerate(sys.argv):		
        	if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
 
	parser = ArgumentParser(description="Write beams of the form (channel,2,2,Nside,Nside)\
						in .npy to {XX, XY, YX, YY} FITS format.")

	add = parser.add_argument
	
	add("-f", "--setup", type=str, required=True, help="add config file to run script")
		
	args = parser.parse_args() 
	
	# ++++ Extract Input Paramters +++++++++++++
		
	config = ConfigParser.ConfigParser()
	config.read('%s' %args.setup)
	azim = float(config.get('AZEL2RADEC', 'azimuth\in_degrees', 1))
	elev = float(config.get('AZEL2RADEC', 'elevation\in_degrees', 1))
	lat = float(config.get('AZEL2RADEC', 'latitude\in_degrees', 1))
	lon = float(config.get('AZEL2RADEC', 'longitude\in_degrees', 1))
	alt = float(config.get('AZEL2RADEC', 'altitute\in_metres', 1))
	date_time = config.get('AZEL2RADEC', 'datetime', 1)
	crval3_strtfreq = float(config.get('AZEL2RADEC', 'start_frequency', 1))
	cdelt3_stepfreq = float(config.get('AZEL2RADEC', 'step_frequency', 1))
	imsize = float(config.get('AZEL2RADEC', 'imagesize\in_degrees', 1))
	beam = config.get('AZEL2RADEC', 'jones_beam_filename', 1)
	output_direc = config.get('AZEL2RADEC', 'save_directory', 1)
	
	# 	
	
	# +++ load Jones beams ++++
	mp = np.load('%s' %beam)
	crpix1 = mp.shape[-1]/2. - 1
	cdelt1 = imsize/mp.shape[-1]
	crpix3 = 1.0
	
	
	# ++++ convert AZEL -->> RADEC ++++++++++++++++	
	
	ra_crval1, dec_crval2 = azel2radec(az = azim, el = elev , lat = lat, 
					   lon = lon, alt = alt, datetime = date_time)
					   
   	# +++++ create headers  +++++++++++++++++++++++

	header = create_FITSheader(CRVAL1 = ra_crval1, CRVAL2 = dec_crval2, 
				CRVAL3 = crval3_strtfreq, CRPIX1 = crpix1 , 
				CDELT1 = cdelt1 , CRPIX3 = crpix3 , CDELT3 = cdelt3_stepfreq)	

	# +++++ Writre to FITS  +++++++++++++++++++++++	
	Gxx = mp[:,0,0,...,...]
	Gxy = mp[:,0,1,...,...]
	Gyx = mp[:,1,0,...,...]
	Gyy = mp[:,1,1,...,...]

	pyfits.writeto('%s_XX_re.fits' %(os.path.basename(beam[:-4])), Gxx.real, header=header, clobber=True)
	pyfits.writeto('%s_XY_re.fits' %(os.path.basename(beam[:-4])), Gxy.real, header=header, clobber=True)
	pyfits.writeto('%s_YX_re.fits' %(os.path.basename(beam[:-4])), Gyx.real, header=header, clobber=True)
	pyfits.writeto('%s_YY_re.fits' %(os.path.basename(beam[:-4])), Gyy.real, header=header, clobber=True)
	pyfits.writeto('%s_XX_im.fits' %(os.path.basename(beam[:-4])), Gxx.imag, header=header, clobber=True)
	pyfits.writeto('%s_XY_im.fits' %(os.path.basename(beam[:-4])), Gxy.imag, header=header, clobber=True)
	pyfits.writeto('%s_YX_im.fits' %(os.path.basename(beam[:-4])), Gyx.imag, header=header, clobber=True)
	pyfits.writeto('%s_YY_im.fits' %(os.path.basename(beam[:-4])), Gyy.imag, header=header, clobber=True)
	
	if os.path.exists('%s' %output_direc):
		os.popen('rm -rf %s' %output_direc)	
	os.mkdir('%s' %output_direc)
	os.system('mv %s*.fits %s' %(os.path.basename(beam[:-4]), output_direc))

#	stop-time
	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0
	print "Total run time: %7.2f hours" % elasped_time
	

if __name__ == '__main__':
	run_header_test()
