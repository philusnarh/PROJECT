#
import pyfits as pf
import numpy as np
#from pyfits import getdata, getheader
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import healpy as hp
from scipy import interpolate
from scipy import ndimage
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.ndimage.interpolation import map_coordinates
import time
import copy
import astLib.astWCS
import cPickle
import ConfigParser
#
def openFitsFile(filename):
	""" Opens a FITS Files using PyFITS
		
        Parameters
        ----------
        filename:str
        path to file  """
        try:
        	hdu = pf.open(filename)        	
        	return hdu
        except:
        	print "Error: cannot open %s"%filename
		raise
#
def mask_beam(bm_data, fwhm, cdelt):
	shape = bm_data.shape
	rad = np.linspace(-shape[0]/2,shape[-1]/2,shape[0])
	rad2d =  np.sqrt(rad[np.newaxis,:]**2+rad[:,np.newaxis]**2)
	mask = rad2d <= fwhm/abs(cdelt)
	masked_beam = bm_data*mask	
	return masked_beam
#
def query_npix(nside, npix, radius):
	print 'searching for nearest pixels:......'
	theta, phi = hp.pix2ang(nside, np.arange(npix))
	vec = np.array(hp.pix2vec(nside, np.arange(npix)))	
	Q_indx = []
	pix_2_ang = []
	for jj in xrange(npix):
		query = hp.query_disc(nside, vec[:,jj].T, np.deg2rad(radius), inclusive=False)
		Q_indx.append(query)
		r1, r2 = hp.pix2ang(nside, query)		
		pix_2_ang.append([r1, r2])		
	print 'Done searching'		
	return 	np.array(Q_indx), np.array(pix_2_ang), theta, phi		
#
def beam_map2skycoord(header, theta_query, phi_query, theta, phi):
	#
	#print 'Centring Beam in Pixel Coordinate'	
	#hdulist = openFitsFile(primary_beam_fits)	
	#header = hdulist[0].header
	#b_data = hdulist[0].data
	#crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
        #crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
        """ 
        * parameter crval      Array of CRVAL keyword (reference) values for each axis.
        * parameter cdelt      Array of CDELT keyword (delta) values for each axis.
        * parameter crpix      Array of CRPIX keyword (reference pixel) values for each axis.  """
        
        
        #	Convert sky coord to pixel coord Using wcstools
	x3 = []	
	theta_dec = np.rad2deg(np.pi/2 - theta)
	phi_ra = np.rad2deg(phi) 
	query_dec = np.rad2deg(np.pi/2 - theta_query)
	query_ra = np.rad2deg(phi_query)	
	header['CRVAL1'] = phi_ra	
	header['CRVAL2'] = theta_dec	
	wcs1 = astLib.astWCS.WCS(header,mode="pyfits")	
	x3.append([wcs1.wcs2pix(i,j) for i, j in zip(query_ra.ravel(), query_dec.ravel())])	
	return np.array(x3)[0, ..., 0], np.array(x3)[0, ..., 1]
#
def fullSky_convolve(healpix_map, primary_beam_fits, ra_n_dec_query_indx, ra_theta, dec_phi, nearest_ipix_indx):
	""" Convolve foregrounds with the beam """
	print 'Full Sky Convolution: ..... !!!'
	#	       
	hdulist = openFitsFile(primary_beam_fits)
	beam_data = hdulist[0].data
	header = hdulist[0].header	
	npix = hp.get_map_size(healpix_map[0, ...])
	#npix = hp.get_map_size(healpix_map)	
	#
	#crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
        #crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
        """ 
        * parameter crval      Array of CRVAL keyword (reference) values for each axis.
        * parameter cdelt      Array of CDELT keyword (delta) values for each axis.
        * parameter crpix      Array of CRPIX keyword (reference pixel) values for each axis.  """
        #
        Xc = []
        Yc = []
	for j in xrange(4):
		print 'started Stokes: %d' %j		
		for iter in xrange(0 + j, 16, 4):		
			outpt = np.zeros(shape = npix, dtype=np.float64)			
			bm_data = beam_data[iter]			
			#			  
			print 'beam mapping:', iter
			if j < 3:			
				for itr in xrange(npix):				
					if iter == 0:
						Xc_beam, Yc_beam = beam_map2skycoord(header,  ra_n_dec_query_indx[itr][0],  ra_n_dec_query_indx[itr][1],\
						ra_theta[itr], dec_phi[itr])
						Xc.append(Xc_beam)
						Yc.append(Yc_beam)								
						bm_map = ndimage.map_coordinates(bm_data.T, [Xc_beam, Yc_beam], mode = 'nearest', order=3)
					else:
						bm_map = ndimage.map_coordinates(bm_data.T, [Xc[itr], Yc[itr]], mode = 'nearest', order=3)				
					weighted_healpix_map = healpix_map[j, ...][nearest_ipix_indx[itr]]				
					outpt[itr] = np.sum(weighted_healpix_map*(bm_map/abs(bm_map).sum()))
					#print itr
			else:				
				outpt = np.zeros(shape = npix, dtype=np.float64)				
			alpha = file('mapp_%d.save'%iter, 'wb')		        
			cPickle.dump(outpt, alpha, protocol = cPickle.HIGHEST_PROTOCOL)
			alpha.close()						
		        #
        measured_map = [] 		
	loaded_objects = []
	for itr4 in xrange(16):		
                alpha = file('mapp_%d.save'%itr4, 'rb')                
                loaded_objects.append(cPickle.load(alpha))
        alpha.close()
        measured_map.append(copy.deepcopy(loaded_objects))
        os.system('rm -f mapp_*.save')	                              
        return measured_map
#
def run_test():
	#
	option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Convolves full-sky map with fully polarized beams')	
	option.add_option('-f','--file1', dest='filename1', 
	                   default='oskar_beams_n_full_sky_conv_setup.ini', action="store", type="string", help='Enter setup file for input parameters')
   	'''option.add_option('-t','--file2', dest='filename2', 
	                  default='fully polarized beams.fits', action="store", type="string", help='Enter the filename of fully polarized beams in FITS format')
	option.add_option('-d','--file3', dest='filename3', 
	                  default='dssetup.ini', action="store", type="string", help='Enter setup file for distorted beams')'''
        #options,args = option.parse_args(sys.argv[1:])
	#path = args[0]
	options,args = option.parse_args()
	filename_input = options.filename1
	#filename_non_distorted = options.filename2
	#filename_distorted = options.filename3
	#
	start = time.time()
	startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Start at %s" % startime
	print 'this will take some time: ++++++++'
	if os.path.exists('%s' %options.filename1):		
		#
		config = ConfigParser.ConfigParser()
		config.read('%s' %options.filename1)
		beam_image_radius = float(config.get('full_sky_convolution', 'radius_of_beam_image', 1))
		mueller_4_x_4_beam_fits = config.get('full_sky_convolution', 'fits_file\input_beam_directory', 1)
		hdu_open = config.get('full_sky_convolution', 'fits_file\input_full_sky_map_directory', 1)
		start_save_as_name = config.get('full_sky_convolution', 'save\start_save_filenames_as', 1)
		results_directory = config.get('full_sky_convolution', 'output_directory', 1)		
		hdulist = openFitsFile('%s' %(hdu_open))
		shape = hdulist[0].data.shape
		if len(shape) == 3:		
			Nside = hp.get_nside(hdulist[0].data[0][0])
			Npix = hp.get_map_size(hdulist[0].data[0][0])	
			#final_conv_map = []	
			#radd = 2.5
			# *******************************	
			#	Query_Disk Function      
			# ********************************
			pix_indx, RA_n_DEC, t1, t2 = query_npix(Nside, Npix, beam_image_radius)		
			#				
			for itr in xrange(shape[0]):	
				foreground_map = hdulist[0].data[itr]				
				conv_foreground_data = fullSky_convolve(foreground_map, '%s' %(mueller_4_x_4_beam_fits), RA_n_DEC, t1, t2, pix_indx)
				#np.save('Disumap_'+ str(itr), conv_foreground_data) 
				#final_conv_map.append(conv_foreground_data)
				# ***************************************
				#	Pickling Convolve map
				# ***************************************
				alpha = file('healpix_map_%d.save' %itr, 'wb')		        
		        	cPickle.dump(conv_foreground_data, alpha, protocol = cPickle.HIGHEST_PROTOCOL)
		        	alpha.close()
		        	print 'Just dumped convolve map channel: %d ' %itr
			# ************************************************
			measured_map = []
			loaded_objects = []
			print 'Loading dumped files:-------' 			
			for iter in xrange(shape[0]):					
				alpha = file('healpix_map_%d.save' %iter, 'rb')				             
				loaded_objects.append(cPickle.load(alpha))				
			alpha.close()
			measured_map.append(copy.deepcopy(np.array(loaded_objects)[:,0,:,:])) 			
			pf.writeto('%s_with_NSIDE_%d_NPIX_%d.fits' %(start_save_as_name, Nside, Npix), np.array(measured_map)[0], clobber=True)
			#
			stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
			print "Stop at %s" % stoptime 
			end = time.time()
			elasped_time = (end - start)/3600.0
			print "Total run time: %7.2f hours" % elasped_time
			#
			print 'Done !!!!!!!!!!!'
		else:
			print "ShapeError: foregroung shape must be (frequency, stokes, npix)" 
			raise
	else:
       		print "Error: setup file for input parameters does not exist"       		
       		raise
       	os.makedirs('%s' %results_directory)
  	os.system('mv -f %s* %s' %(start_save_as_name, results_directory))
if __name__ == '__main__':
	run_test()


	

 

