#
import pyfits as pf
import numpy as np
#from pyfits import getdata, getheader
#from numpy import linalg as la
from matplotlib import pylab as plt
from optparse import OptionParser
import sys,os
import healpy as hp
from scipy import interpolate
from scipy import ndimage, signal
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.ndimage.interpolation import map_coordinates
#import pywcs
import time
import copy
import wcstools
import cPickle
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
def normalize(lst):
	s = np.sum(lst)
	return map(lambda v: float(v)/s, lst)
#
def cwt(data, wavelet, width):
	wavelet_data = wavelet(min(10*width, len(data)), width)
	print wavelet_data
	outcome = np.convole(data, wavelet_data, mode = 'sam')
	return outcome
#
def fullSky_convolve(healpix_map, radius, primary_beam_fits):
	""" Convolve foregrounds with the beam """
	#
	measured_map = []
        #temp = []
        #gauss = Gaussian1DKernel(stddev=2)
	hdulist = openFitsFile(primary_beam_fits)
	beam_data = hdulist[0].data
	header = hdulist[0].header
	#ra = np.deg2rad(header['CRVAL1'])
	#dec = np.deg2rad(header['CRVAL2'])
	nside = hp.get_nside(healpix_map[0, ...])
	npix = hp.get_map_size(healpix_map[0, ...])			# total number of pixels in the map must be  12 * nside^2 
	
	#	ipix = np.arange(npix)
	#
	crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
        crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
        """ 
        * parameter crval      Array of CRVAL keyword (reference) values for each axis.
        * parameter cdelt      Array of CDELT keyword (delta) values for each axis.
        * parameter crpix      Array of CRPIX keyword (reference pixel) values for each axis.  """
        #
	for j in xrange(4):
		print 'started Stokes: %d' %j	
		for iter in xrange(0 + j, 16, 4):
			outpt = np.zeros(shape = npix, dtype=np.float64)
			# mask beam
			bm_data = beam_data[iter]
			#masked_beam= beam_data[iter]
			shape = bm_data.shape
			rad = np.linspace(-shape[0]/2,shape[-1]/2,shape[0])
			rad2d =  np.sqrt(rad[np.newaxis,:]**2+rad[:,np.newaxis]**2)
			mask = rad2d <= radius/abs(cdelt2)
			masked_beam = bm_data*mask
			#        	        		
			for itr in xrange(npix):							#	npix
				theta, phi = hp.pix2ang(nside, itr)			#     coordinate at which the power is being calculated		
				d1_indx = hp.get_all_neighbours(nside, theta, phi)		
				d1_indx = np.where(d1_indx < 0, itr, d1_indx)
				r1, r2 = hp.pix2ang(nside, d1_indx)
				d2_indx = hp.get_all_neighbours(nside, r1, r2)		
				d2_indx = np.where(d2_indx < 0, itr, d2_indx)
				n1, n2 = hp.pix2ang(nside, d2_indx)
				disk_ipix_indx = hp.get_all_neighbours(nside, n1, n2).ravel()
				disk_ipix_indx = np.where(disk_ipix_indx < 0, itr, disk_ipix_indx)	
				#	++++++++++++++++++++++++++++++++++++
				#vec = hp.pix2vec(nside, itr)                        
				#disk_ipix_indx = hp.query_disc(nside, vec, np.deg2rad(radius), inclusive = True)      #	nearest pixels at specific power coordinate
				#print len(disk_ipix_indx)
				#h_map = np.zeros(shape = len(disk_ipix_indx), dtype=np.float64)
								
				t1, t2 = hp.pix2ang(nside, disk_ipix_indx)
				#d1_indx = hp.get_all_neighbours(nside, r1, r2)
				#d1_indx = np.where(d1_indx < 0, itr, d1_indx)
				#t1, t2 = hp.pix2ang(nside, d1_indx)
				#for itr1 in xrange(len(disk_ipix_indx)):					
				#s = interpolate.NearestNDInterpolator((t1[:,itr1], t2[:,itr1]), healpix_map[j, ...][d1_indx[:,itr1]])
				#h_map[itr1] = s(r1[itr1], r2[itr1])
				#s = interpolate.NearestNDInterpolator((t1.ravel('F'), t2.ravel('F')), healpix_map[j, ...][d1_indx.ravel('F')])
				#h_map = s(r1, r2)
				#print len(h_map)
				#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				# disk_ipix_indx = hp.get_all_neighbours(nside, theta, phi)
				# disk_ipix_indx = np.where(disk_ipix_indx < 0, itr, disk_ipix_indx)								
				#ra_theta, dec_phi = hp.pix2ang(nside, disk_ipix_indx)
				#weighted_healpix_map = hp.pixelfunc.get_interp_val(healpix_map[j, ...], ra_theta, dec_phi, nest=False)
				#widths = np.linspace(0,1, len(h_map)/2)
				#print widths.size
				t = np.linspace(-3.,3., len(disk_ipix_indx)/2)
				T = np.sinc(t)*np.sinc(t/3.)
				weighted_healpix_map = np.convolve(healpix_map[j, ...][disk_ipix_indx], T/T.sum(), mode='same')	
				#print weighted_healpix_map
				# weighted_healpix_map = np.convolve(healpix_map[j, ...][disk_ipix_indx], T/T.sum(), mode='same')	
				# weighted_healpix_map = np.convolve(h_map, T/T.sum(), mode='same')	
				# pixel coordinates
				#ra_theta, dec_phi = hp.pix2ang(nside, disk_ipix_indx)
				bx = theta - t1
				by = phi - t2
				mt = wcstools.wcs(primary_beam_fits, ext=0,rot=0,cd=False)
				wx, wy = mt.xy2sky(crpix1, crpix2)
				bx1 = np.rad2deg(bx)/180.  +  wx
				#print bx1				
				by1 = np.rad2deg(by)/360. +  wy
				#print by1
				s1, s2 = np.array(mt.sky2xy(bx1, by1))
				ra_indx = s1.astype(int)
				dec_indx = s2.astype(int)								
				#	Apply Convolution
				#print 'iter %d' % itr     	
				#outpt[itr] = np.sum(healpix_map[j, ...][disk_ipix_indx]*(beam_data[iter][beam_ipix_indx[0, :], beam_ipix_indx[1,:]]))
				#break 
				#outpt[itr] = np.sum(healpix_map[j, ...][disk_ipix_indx]*(beam_data[iter][ra_indx, dec_indx]/(beam_data[iter][ra_indx, dec_indx]).sum()))
				#z = convolve(masked_beam[ra_indx, dec_indx]/masked_beam[ra_indx, dec_indx].sum(),gauss)
				#hmap = ndimage.filters.gaussian_filter(healpix_map[j, ...][disk_ipix_indx], sigma = 1)
				#outpt[itr] = np.sum(healpix_map[j, ...][disk_ipix_indx]*(masked_beam[ra_indx, dec_indx]/masked_beam[ra_indx, dec_indx].sum()))	#
				outpt[itr] = np.sum(weighted_healpix_map*(masked_beam[ra_indx, dec_indx]/masked_beam[ra_indx, dec_indx].sum()))				
				#outpt[itr] = np.sum(hmap*(masked_beam[ra_indx, dec_indx]/masked_beam[ra_indx, dec_indx].sum()))
				#outpt[itr] = np.sum(healpix_map[j, ...][disk_ipix_indx]*(bm_data[ra_indx, dec_indx]/bm_data[ra_indx, dec_indx].sum()))		        
		        #z = convolve(outpt1,gauss)
		        alpha = file('disobjects%d.save'%iter, 'wb')
		        #h_map = ndimage.filters.gaussian_filter(outpt, sigma = 3.)
		        cPickle.dump(outpt, alpha, protocol = cPickle.HIGHEST_PROTOCOL)
		        alpha.close()
		        print 'Just dumped truobjects%d.save:-------'%iter
		        
	print 'Loading dumped files:-------' 	
	loaded_objects = []
	for itr4 in xrange(16):		
                alpha = file('disobjects%d.save'%itr4, 'rb')                
                loaded_objects.append(cPickle.load(alpha))
        alpha.close()
        measured_map.append(copy.deepcopy(loaded_objects)) 		
	#np.save('nearmap%d' %j, np.array(measured_map))                               
        return measured_map
#
def run_test():
	start = time.time()
        startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Start at %s" % startime
        print 'this will take some time: ++++++++'
	#hdu_gridd = openFitsFile('560_4_786432_K2JanskypPixel.fits')
	hdu_gridd = openFitsFile('h11-16.fits')
	#stokes_p = hdu_gridd[0].data[0]
	trumap = []
	#dismap = []
	radd = 150.		
	for itr in xrange(hdu_gridd[0].data.shape[0]):
		stokes_p = hdu_gridd[0].data[itr]
		#tru_image_data = fullSky_convolve(stokes_p,radd/60. ,'mueller_4_x_4_beamw.fits', 'chan_%d_cmap' %itr)
		tru_image_data = fullSky_convolve(stokes_p,radd/60. ,'dis_mueller_4_x_4_beam.fits')
		np.save('dschan_Trrmap'+ str(itr), tru_image_data) 
		trumap.append(tru_image_data)
	pf.writeto('dconv_map_560_x_16_786432_trubeam.fits', np.array(trumap), clobber=True)
	#radd = 150.
	#tru_image_data = fullSky_convolve(stokes_p,radd/60. ,'mueller_4_x_4_beamw.fits', 'tcnear%dmap' %radd) 		#	np.deg2rad(radd/60.)
        #np.save('tcfullsky%d' %radd, np.array(tru_image_data))
        #
	'''for radd in xrange(30, 160, 30):
		tru_image_data = fullSky_convolve(stokes_p,radd/60. ,'mueller_4_x_4_beamw.fits', 'mskcnear%dmaptru' %radd/60.) 		#	np.deg2rad(radd/60.)
                np.save('corfullsky%d' %radd, np.array(tru_image_data))'''
        #
	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Stop at %s" % stoptime 
        end = time.time()
        elasped_time = (end - start)/3600.0
        print "Total run time: %7.2f hours" % elasped_time
        #
        print 'Done !!!!!!!!!!!'
if __name__ == '__main__':
	run_test()


	

 

