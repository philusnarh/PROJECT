#!/usr/bin/env python
"""
beampattern_n_interferometer_sim.py
===================================

This script generates fully polarized primary beams and visibilities from OSKAR 2.6.1

"""
import glob
#import pyfits
from astropy.io import fits as pyfits
from optparse import OptionParser
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import os, sys
import shutil
import ConfigParser
import time
import datetime
import cPickle
import ephem
import pylab as plt
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
		print "\n >> Error: cannot open %s"%filename
		raise
#

def ampphase2complex(amplitude, phase):
    return amplitude * np.exp(1j*phase)
    


#def write_fits(dataset, beamcube, filename, freqs, diameter, scanner, tracker):
#    # Calculate parameters from the dataset
#    el = dataset.env_el[0]
#    az = np.mean(dataset.scanaz*(180./np.pi))
#    ID = dataset.filename.split('/')[-1][:-3]
#    
#    # Create header
#    hdr = fits.Header()
#    hdr['ID'] = ID
#    ctypes = ['AZIMUTH', 'ELEVATION', 'JONES0', 'JONES1', 'FREQ']
#    crvals = [az, el, 0,0, freqs[0]/1e3]
#    cdelts = [diameter/beamcube.shape[-2], diameter/beamcube.shape[-1], 1,1, freqs[1]/1e3-freqs[0]/1e3]
#    cunits = ['deg', 'deg', '', '', 'GHz']
#    for i in range(len(beamcube.shape)):
#        ii = str(i+1)
#        hdr['CTYPE'+ii] = ctypes[i]
#        hdr['CRVAL'+ii] = crvals[i]
#        hdr['CDELT'+ii] = cdelts[i]
#        hdr['CUNIT'+ii] = cunits[i]
#    hdr['TELESCOP'] = 'MeerKAT'
#    hdr['OBSSTART'] = dataset.rawtime[0]
#    hdr['OBSEND'] = dataset.rawtime[-1]
#    hdr['DURATION'] = (dataset.rawtime[-1]-dataset.rawtime[0])/3600.
#    hdr['SCANNER'] = scanner
#    hdr['TRACKER'] = tracker
#    hdr['TARGET'] = dataset.target.name
#    hdr['TARGRA'] = dataset.target.radec()[0]*(180./np.pi)
#    hdr['TARGDEC'] = dataset.target.radec()[1]*(180./np.pi)
#    hdr['DATE'] = str(datetime.datetime.now())
#    
#    # Write real and imag parts of data
#    hdr['PART'] = 'REAL'
#    hdu_r = fits.PrimaryHDU(beamcube.real, header=hdr)
#    hdr['PART'] = 'IMAG'
#    hdu_i = fits.PrimaryHDU(beamcube.imag, header=hdr)
#    hdu_r.writeto(filename+'_real.fits', clobber=True)
#hdu_i.writeto(filename+'_imag.fits', clobber=True)


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
	# az  = -35.000262416  #deg
	# el  =  61.0 #51.0910024663 #deg
	# lat = -30.7214 #deg
	# lon = 21.4106 #deg
	# alt = 1038  # m
	# #ut  = 2455822.20000367 #julian date
	# date =  "2017/02/23 03:28:04.182" # '2000/01/01 12:00:00'  
	# Which Julian Date does Ephem start its own count at?
	# J0 = ephem.julian_date(0)

	meerkat_observer = ephem.Observer()
	meerkat_observer.lon = lon  
	meerkat_observer.lat = lat  
	meerkat_observer.elevation = alt
	meerkat_observer.date = datetime
	meerkat_observer.epoch = ephem.J2000
	#meerkat_observer.pressure = 0
	#ut - J0

	#print observer.date
	print meerkat_observer.radec_of(float(ephem.degrees(az)),
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



def beam_model(xx_real, xx_imag, xy_real, xy_imag, yx_real, yx_imag, yy_real, yy_imag, header):
	    #
	    '''
	    Extracts the Jones components into XX, XY, YX, YY
	    to produce Mueller components'''
	    
	    #
	    #	Generating Jones Terms
	    #
	    xx = xx_real + 1j*xx_imag
	    xy = xy_real + 1j*xy_imag
	    yx = yx_real + 1j*yx_imag
	    yy = yy_real + 1j*yy_imag

	    #	    
	    #	Generating Mueller Terms
	    #	   
	    m_ii = 0.5*(xx*np.conjugate(xx) + xy*np.conjugate(xy) + yx*np.conjugate(yx) + yy*np.conjugate(yy))
	    m_iq = 0.5*(xx*np.conjugate(yx) - xy*np.conjugate(xx) + yx*np.conjugate(yy) + yy*np.conjugate(xy)) 
	    m_iu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) + yx*np.conjugate(yy)+ yy*np.conjugate(yx)) 
	    m_iv = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
	    m_qi = 0.5*(xx*np.conjugate(yx) - xy*np.conjugate(xx) + yx*np.conjugate(yy) + yy*np.conjugate(xy))   
	    m_qq = 0.5*(xx*np.conjugate(xx) - xy*np.conjugate(xy) - yx*np.conjugate(yx)+ yy*np.conjugate(yy)) 
	    m_qu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx)) 
	    m_qv = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
	    m_ui = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) + yx*np.conjugate(yy)+ yy*np.conjugate(yx))
	    m_uq = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx))  
	    m_uu = 0.5*(xx*np.conjugate(yy) + yy*np.conjugate(xx) + xy*np.conjugate(yx)+ yx*np.conjugate(xy))
	    #m_uv = 0.5*1j*(xx*np.conjugate(yy) + yx*np.conjugate(xy) - yy*np.conjugate(xx) - xy*np.conjugate(yx)) 
	    #m_uv = 0.5*1j*(xx*np.conjugate(yy) - yy*np.conjugate(xx) - xy*np.conjugate(yx) - xy*np.conjugate(xy))
	    m_uv = 0.5*1j*(xx*np.conjugate(yy) + yy*np.conjugate(xx) + xy*np.conjugate(yx) + xy*np.conjugate(xy))  
	    m_vi = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))  
	    m_vq = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
	    #m_vu = 0.5*1j*(xx*np.conjugate(yy) + yx*np.conjugate(xy) - yy*np.conjugate(xx) - xy*np.conjugate(yx))
	    m_vu = 0.5*1j*(xx*np.conjugate(yy) + yy*np.conjugate(xx) + xy*np.conjugate(yx) + xy*np.conjugate(xy))
	    #m_vu = 0.5*1j*(-xx*np.conjugate(yy) + yy*np.conjugate(xx) - xy*np.conjugate(yx) + xy*np.conjugate(xy))    
	    m_vv = 0.5*(xx*np.conjugate(yy) - yx*np.conjugate(xy) + yy*np.conjugate(xx) - xy*np.conjugate(yx)) 
	    #
	    M = []
	    M.append(m_ii.real) 
	    M.append(m_iq.real) 
	    M.append(m_iu.real)
	    M.append(m_iv.real)
	    M.append(m_qi.real) 
	    M.append(m_qq.real)
	    M.append(m_qu.real)
	    M.append(m_qv.real)
	    M.append(m_ui.real)
	    M.append(m_uq.real)
	    M.append(m_uu.real)
	    M.append(m_uv.real)
	    M.append(m_vi.real)
	    M.append(m_vq.real)
	    M.append(m_vu.real)
	    M.append(m_vv.real) 
	    M = np.array(M)
	    
	   	   
	    return   np.swapaxes(M,0,1)     # np.swapaxes(M,1,2)	
	    #

def headers(Fits_filename):
	'''Gets FITS header'''
	
    	hd = pyfits.getheader(Fits_filename)
    
	return hd 
	#
	
    
def run_sim_test():
	#	
	#	start-time
	start = time.time()
	startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Start at %s" % startime

	# ++++++++++++ INPUT PARAMETERS ++++++++

	crval3 = 900.0e6
	crpix1 = 512/2. - 1
	cdelt1 = 6.0/512.
	crpix3 = 1.0
	cdelt3 = 1.0

	ra_crval1, dec_crval2 = azel2radec(az = -35.000262416, el = 51.0910024663 , lat = -30.7214, 
									   lon = 21.4106, alt = 1038, datetime = '2017/02/23 03:28:04.182')

	header = create_FITSheader(CRVAL1 = ra_crval1, CRVAL2 = dec_crval2, 
			CRVAL3 = crval3, CRPIX1 = crpix1 , CDELT1 = cdelt1 , CRPIX3 = crpix3 , CDELT3 = cdelt3)
	#	+++++++++++++++ load ant-beams  +++++++++++++++++++++++++++++++++++++++++++++++++

	bm = ['/home/narh/meerkat-hbeam/1487813282_m000_900MHz_1MHz_500channels_Jones.npy',
		  '/home/narh/meerkat-hbeam/1487813282_m012_900MHz_1MHz_500channels_Jones.npy',
		  '/home/narh/meerkat-hbeam/1487813282_m017_900MHz_1MHz_500channels_Jones.npy']

 	
 	
	#m12 = np.load('/home/narh/meerkat-hbeam/1487813282_m012_900MHz_1MHz_500channels_Jones.npy')
	#m17 = np.load('/home/narh/meerkat-hbeam/1487813282_m017_900MHz_1MHz_500channels_Jones.npy')

	for iter in xrange(len(bm)):
		print '\n >> %s ' %(os.path.basename(bm[iter][:-4]))
		mp = np.load('%s' %bm[iter])
		Gxx = mp[:,0,0,...,...]
		Gxy = mp[:,0,1,...,...]
		Gyx = mp[:,1,0,...,...]
		Gyy = mp[:,1,1,...,...]

		pyfits.writeto('%s_XX_re.fits' %(os.path.basename(bm[iter][:-4])), Gxx.real, header=header, clobber=True)
		pyfits.writeto('%s_XY_re.fits' %(os.path.basename(bm[iter][:-4])), Gxy.real, header=header, clobber=True)
		pyfits.writeto('%s_YX_re.fits' %(os.path.basename(bm[iter][:-4])), Gyx.real, header=header, clobber=True)
		pyfits.writeto('%s_YY_re.fits' %(os.path.basename(bm[iter][:-4])), Gyy.real, header=header, clobber=True)
		pyfits.writeto('%s_XX_im.fits' %(os.path.basename(bm[iter][:-4])), Gxx.imag, header=header, clobber=True)
		pyfits.writeto('%s_XY_im.fits' %(os.path.basename(bm[iter][:-4])), Gxy.imag, header=header, clobber=True)
		pyfits.writeto('%s_YX_im.fits' %(os.path.basename(bm[iter][:-4])), Gyx.imag, header=header, clobber=True)
		pyfits.writeto('%s_YY_im.fits' %(os.path.basename(bm[iter][:-4])), Gyy.imag, header=header, clobber=True)



	
	"""
	Gx = m00[0,0,1].real

	
	option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Run OSKAR beam pattern and interferometer simulations')
	
	option.add_option('-f','--file1', dest='setup1', 
	                  default='beamsetup.ini', action="store", type="string",
	                  help='Enter config file for input parameters ')
	                  
  	option.add_option('-t','--file2', dest='setup2', 
	                  default='tru_oskar_config.ini', action="store", type="string", 
	                  help='Enter OSKAR setup file for simulation')
	
	options,args = option.parse_args()
	#
	
	# run OSKAR
	os.system('oskar_sim_beam_pattern %s' %options.setup2)
	
	XXamp = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_AMP_XX.txt'
	XYamp = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_AMP_XY.txt'
	YXamp = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_AMP_YX.txt'
	YYamp = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_AMP_YY.txt'	
	XXphse = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_PHASE_XX.txt'
	XYphse = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_PHASE_XY.txt'
	YXphse = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_PHASE_YX.txt'
	YYphse = 'beam_pattern_S0000_TIME_SEP_CHAN_SEP_PHASE_YY.txt'
	
	
	
	amp = np.loadtxt('%s' %XXamp)
	phs = np.loadtxt('%s' %XXphse)
	
	cmplx = ampphase2complex(amplitude=amp, phase=phs)
	xx = (cmplx.imag).reshape(128,128)
	xxr = (cmplx.real).reshape(128,128)
	
	plt.figure(1);plt.imshow(xx)
	plt.colorbar()
	plt.savefig('fig.png')
	plt.close()
	plt.figure(2);plt.imshow(xxr)
	plt.colorbar()
	plt.savefig('fig2.png')
	plt.close()
#	plt.show()
	print 'Done !!!!!!!!!!!'
	
	
	# ++++++++++  Run beam pattern simulation  +++++++++++
			
	if os.path.exists('%s' %options.setup1):
		config = ConfigParser.ConfigParser()
		config.read('%s' %options.setup1)
		start_save_file = config.get('simulations', 'start_saving_beam_names_as', 1)
		beamname = 'jones'
		
		if config.get('simulations', 'run_beam_pattern_sim', 1) == 'true':			
			
			if os.path.exists('%s' %options.setup2):
				ss = ConfigParser.ConfigParser()
				ss.read('%s' %(options.setup2))
				nchanns = int(ss.get('observation', 'num_channels', 1))
				start_freq = float(ss.get('observation', 'start_frequency_hz', 1))
				sf = start_freq
				freq_inc = float(ss.get('observation', 'frequency_inc_hz', 1))
				beampath = ss.get('beam_pattern', 'root_path', 1)
				
				
				for n in xrange(nchanns):
					# iterate channel and run beam simulation
					os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, 1))					
					os.system('oskar_sim_beam_pattern %s' %options.setup2)
					
					#	station beams					
			    		count_beams = len(glob.glob('*RAW_COMPLEX.txt'))
				    	#
				    	print "\n >> Opening Amplitude & Phase Beampattern:"
				    	for cc in xrange(count_beams):	    		
				    	
				    		complex_pattern = np.loadtxt('%s_S%s_TIME_SEP_CHAN_SEP_RAW_COMPLEX.txt'  %(beampath, str(cc).zfill(4)))	
				    		alpha = file('station_%s_chan_%d.save'%(str(cc).zfill(4), n), 'wb')		        
						cPickle.dump(complex_pattern, alpha, protocol = cPickle.HIGHEST_PROTOCOL)
						alpha.close()
					start_freq += freq_inc
					print 'start_freq', start_freq 
					os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, start_freq))	    	
			    	
			    	os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, sf))
			    	os.system('oskar_settings_set  %s observation/num_channels %s' %(options.setup2, nchanns))
			    	
			    	#	header
			    	
			    	
			    	header = headers('%s_S%s_TIME_SEP_CHAN_SEP_AMP_XX.fits' %(beampath, str(0).zfill(4)))
#			    	header['CRVAL3'] = sf 
#			    	header.update('CDELT3', freq_inc, after='CRVAL3')
			    	#header.update['CDELT3'] = freq_inc
			    	
			    	#	Checking for number of station beams			    	
			    	
			    	beam_size = int(ss.get('beam_pattern', 'beam_image\size', 1))
			    	
			    	for cc in xrange(count_beams):
			    		
			    		xx_re = []
				    	xx_im = []
				    	xy_re = []
				    	xy_im = []
				    	yx_re = []
				    	yx_im = []
				    	yy_re = []
				    	yy_im = []
				    	#	 
			    	
				    	for n in xrange(nchanns): 				    					    		
				    		
				    		complex_pattern = np.load('station_%s_chan_%d.save'%(str(cc).zfill(4), n))
					    	 
					    	if complex_pattern[:, 0].reshape(beam_size, beam_size).max() < 0.9:
					    		x1 = -complex_pattern[:, 0].reshape(beam_size, beam_size)
					    		xx_re.append(x1/x1.max())	# 
				    		else:
				    			xx_re.append(complex_pattern[:, 0].reshape(beam_size, beam_size))
					    	xx_im.append(complex_pattern[:, 1].reshape(beam_size, beam_size))
					    	xy_re.append(complex_pattern[:, 2].reshape(beam_size, beam_size))
					    	xy_im.append(complex_pattern[:, 3].reshape(beam_size, beam_size))
					    	yx_re.append(complex_pattern[:, 4].reshape(beam_size, beam_size))
					    	yx_im.append(complex_pattern[:, 5].reshape(beam_size, beam_size))
					    	if complex_pattern[:, 6].reshape(beam_size, beam_size).max() < 0.9:
					    		yy_re.append(-complex_pattern[:, 6].reshape(beam_size, beam_size))
				    		else: 
					    		yy_re.append(complex_pattern[:, 6].reshape(beam_size, beam_size))
					    	
					    	yy_im.append(complex_pattern[:, 7].reshape(beam_size, beam_size))
					    	
				        xx_re = np.array(xx_re)
				    	xx_im = np.array(xx_im)
				    	xy_re = np.array(xy_re)
				    	xy_im = np.array(xy_im)
				    	yx_re = np.array(yx_re)
				    	yx_im = np.array(yx_im)
				    	yy_re = np.array(yy_re)
				    	yy_im = np.array(yy_im)
				    	
				    	pyfits.writeto('%s_%s_S%s_XX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_re, header=header, clobber=True)
			   		pyfits.writeto('%s_%s_S%s_XY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_im, header=header, clobber=True)
#					break
#					mdata = []
#					for n in xrange(nchanns):
#						mueller_data = beam_model(np.sqrt(abs(xx_re[n])), xx_im[n], 
#									xy_re[n], xy_im[n], yx_re[n], yx_im[n],
#									np.sqrt(abs(yy_re[n])), yy_im[n], header)
#						mdata.append(mueller_data)
#						
#					print 'mdata shape; ', np.array(mdata).shape		    	
#				    		    	
#				    	#raise
#					pyfits.writeto('%s_S%s_mueller_beam.fits' %(start_save_file, str(cc).zfill(4)),
#					               np.array(mdata), header=header, clobber=True)
#					#					    	
					mueller_data = beam_model(np.sqrt(abs(xx_re)), xx_im, xy_re, xy_im, yx_re, yx_im, np.sqrt(abs(yy_re)), yy_im, header)
						    	
				    		    	
				    	pyfits.writeto('%s_%s_S%s_XX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_re, header=header, clobber=True)
			   		pyfits.writeto('%s_%s_S%s_XY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_im, header=header, clobber=True)
					pyfits.writeto('%s_S%s_mueller_beam.fits' %(start_save_file, str(cc).zfill(4)), mueller_data, header=header, clobber=True)
#					#
					
			else:
				
				print "\n >> Error: OSKAR setup file does not exist"       		
	       			raise				
		
		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		
		#	Creating Output directory
				
		filename = config.get('simulations', 'output_directory', 1)
		if os.path.isdir("%s" %filename):		
	       		os.popen('rm -rf %s' %filename)
		os.makedirs('%s' %filename)
		
		#	Saving Beam-patterns in output directory
		
		if os.path.exists('%s_%s_S%s_XX_re.fits' %(start_save_file, beamname, str(0).zfill(4))):
				
		 	os.system('mv -f %s* %s'  %(start_save_file, filename))
		 	
		 	#	removing temporary  directories
		 	
		        os.system('rm -f *.%s %s' %('fits', '*RAW_COMPLEX.txt'))
		        os.system('rm -f station*.save')
		        os.system('rm -f *.log')
		        
		#  +++++++++++ Run interferometer Simulation	+++++++++++
		
		if config.get('simulations', 'run_interferometer_sim', 1) == 'true':
			if os.path.exists('%s' %options.setup2):
				ss = ConfigParser.ConfigParser()
				ss.read('%s' %(options.setup2))	
						
				integration_time = int(ss.get('interferometer', 'time_average_sec', 1))
				beamtime_slots = int(ss.get('observation', 'num_time_steps', 1))
				beamform_len = ss.get('observation', 'length', 1)
				st = config.get('simulations', 'run_interferometer_sim\Synthesis_steps', 1) 
				msname = config.get('simulations', 'run_interferometer_sim\msname', 1)			
				x = time.strptime(st.split(',')[0],'%H:%M:%S.%f')
							
				date_time = datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()			
				os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, int(date_time/integration_time) + 1))
				os.system('oskar_settings_set  %s observation/length %s' %(options.setup2, st))
				nchanns = int(ss.get('observation', 'num_channels', 1))
				start_freq = float(ss.get('observation', 'start_frequency_hz', 1))
				sf = start_freq
				freq_inc = float(ss.get('observation', 'frequency_inc_hz', 1))
				visname = ss.get('interferometer', 'oskar_vis_filename', 1)
				
				if nchanns > 50:
					for iter in xrange(nchanns):
						os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, 1))					
						print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
						print "\n >> OSKAR Interferometer simulation starts channel %d:" %(iter +1) 
						print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
						
						os.system('oskar_sim_interferometer %s' %options.setup2)
						start_freq += freq_inc
						os.system('oskar_settings_set  %s observation/start_frequency_hz %d' %(options.setup2, start_freq))
						os.system('oskar_vis_to_ms %s -o %s_chan%s.MS' %(visname, msname, str(iter).zfill(4)))				
					
					os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, beamtime_slots))
					os.system('oskar_settings_set  %s observation/length %s' %(options.setup2, beamform_len))
					os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, sf))
					os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, nchanns))
					
				else:
					print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
					print "\n >> OSKAR Interferometer simulation starts " 
					print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
					
					os.system('oskar_sim_interferometer %s' %options.setup2)
					os.system('oskar_vis_to_ms %s -o %s.MS' %(visname, msname))
					os.system('mv -f *.MS %s' %filename )
				        os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, beamtime_slots))
				        os.system('oskar_settings_set  %s observation/length %s' %(options.setup2, beamform_len))
		        
		        else:
				
				print "\n >> Error: OSKAR setup file does not exist"       		
	       			raise					     					
				
		
	else:
       		print "\n >> Error: config file for input parameters does not exist"       		
       		raise
	
	#	stop-time
	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0
	print "Total run time: %7.2f hours" % elasped_time
	"""

if __name__ == '__main__':
	run_sim_test()
	
		   	
