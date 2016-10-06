#!/usr/bin/env python
"""
beampattern_n_interferometer_sim.py
===================================

This script generates fully polarized primary beams and visibilities from OSKAR 2.6.1

"""
import glob
import pyfits
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
	    m_uv = 0.5*1j*(xx*np.conjugate(yy) + yx*np.conjugate(xy) - yy*np.conjugate(xx) - xy*np.conjugate(yx))  
	    m_vi = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))  
	    m_vq = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
	    m_vu = 0.5*1j*(xx*np.conjugate(yy) + yx*np.conjugate(xy) - yy*np.conjugate(xx) - xy*np.conjugate(yx))   
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
	    
	   	   
	    return np.swapaxes(M,0,1)    	
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
	
	option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Run OSKAR beam pattern and interferometer simulations')
	
	option.add_option('-f','--file1', dest='setup1', 
	                  default='beamsetup.ini', action="store", type="string",
	                  help='Enter config file for input parameters ')
	                  
  	option.add_option('-t','--file2', dest='setup2', 
	                  default='oskar_config.ini', action="store", type="string", 
	                  help='Enter OSKAR setup file for simulation')
	
	options,args = option.parse_args()
	#
	
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
				
				print "\n >> OSKAR Beam pattern simulation starts:"
				os.system('oskar_sim_beam_pattern %s' %options.setup2)			
				
	
				#	Checking for number of station beams
			    	#
			    	
			    	count_beams = len(glob.glob('*AMP_XX.fits'))
			    	#
			    	print "\n >> Opening Amplitude & Phase Beampattern:"
			    	for cc in xrange(count_beams):	    		
				    	beam_filename = 'beam_pattern_S%s_TIME_SEP_CHAN_SEP_AMP_' %(str(cc).zfill(4))
				    	phase_filename = 'beam_pattern_S%s_TIME_SEP_CHAN_SEP_PHASE_' %(str(cc).zfill(4))				    	
				    	
				    	
				    	v = ['XX.fits', 'XY.fits', 'YX.fits', 'YY.fits']
				    	
			    		#	Extracting the Amplitude and Phase
			    		
				    	for iter in xrange(4):
				    		hduAmp = openFitsFile('%s%s' %(beam_filename, v[iter])) 
				    		hduPh = openFitsFile('%s%s' %(phase_filename, v[iter]))
				    		
				    		if iter == 0:
				    			j_xx = hduAmp[0].data[0]*np.exp(1j*hduPh[0].data[0])
				    		
				    		elif iter == 1:
				    			j_xy = hduAmp[0].data[0]*np.exp(1j*hduPh[0].data[0])
				    			
			    			elif iter == 2:
				    			j_yx = hduAmp[0].data[0]*np.exp(1j*hduPh[0].data[0])
				    			
			    			else:
				    			j_yy = hduAmp[0].data[0]*np.exp(1j*hduPh[0].data[0])
				    	
				    	xx_re = []
				    	yy_re = []
				    	
				    	for n in xrange(j_xx.shape[0]): 				    					    		
				    		
				    		#complex_pattern = np.load('station_%s_chan_%d.save'%(str(cc).zfill(4), n))
					    	 
					    	if np.any(j_xx.real[n] < 0):
					    		xx_re.append(abs(j_xx.real[n]))
					    		'''import pylab as plt
						    	print np.array(xx_re).min()
						    	plt.imshow(np.sqrt(abs(np.array(xx_re)[0])))
						    	plt.colorbar()
						    	plt.show()
					    		raise'''
				    		else:
				    			xx_re.append(j_xx.real[n])
				    			
			    			if np.any(j_yy.real[n] < 0):			    				
			    				
					    		yy_re.append(abs(j_yy.real[n]))
					    		
				    		else:
				    			yy_re.append(j_yy.real[n])					    	
					    	
				        	
		    			xx_re = np.array(xx_re)
				    	xx_im = j_xx.imag
				    	xy_re = j_xy.real
				    	xy_im = j_xy.imag
				    	yx_re = j_yx.real
				    	yx_im = j_yx.imag
				    	
				    	yy_re = np.array(yy_re)
				    	yy_im = j_yy.imag
				    	
				    	#
				    	
				    	#	header
				    	#
				    	header = headers('beam_pattern_S%s_TIME_SEP_CHAN_SEP_AMP_XX.fits' %(str(0).zfill(4)))
				    	
				    	pyfits.writeto('%s_%s_S%s_XX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_re, header=header, clobber=True)
			   		pyfits.writeto('%s_%s_S%s_XY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_re.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_re, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_XY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), xy_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YX_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yx_im, header=header, clobber=True)
					pyfits.writeto('%s_%s_S%s_YY_im.fits' %(start_save_file, beamname, str(cc).zfill(4)), yy_im, header=header, clobber=True)
					
					mueller_data = beam_model((xx_re), xx_im, xy_re, xy_im, yx_re, yx_im, (yy_re), yy_im, header)
					#mueller_data = beam_model(np.sqrt(abs(xx_re)), xx_im, xy_re, xy_im, yx_re, yx_im, np.sqrt(abs(yy_re)), yy_im, header)
					
					pyfits.writeto('%s_S%s_mueller_beam.fits' %(start_save_file, str(cc).zfill(4)), mueller_data, header=header, clobber=True)
				""" 
				
				for n in xrange(nchanns):
					os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, 1))	
					
					os.system('oskar_sim_beam_pattern %s' %options.setup2)					
			    		count_beams = len(glob.glob('*RAW_COMPLEX.txt'))
				    	#
				    	print "\n >> Opening Amplitude & Phase Beampattern:"
				    	for cc in xrange(count_beams):	    		
				    	
				    		complex_pattern = np.loadtxt('beam_pattern_S%s_TIME_SEP_CHAN_SEP_RAW_COMPLEX.txt'  %(str(cc).zfill(4)))	
				    		alpha = file('station_%s_chan_%d.save'%(str(cc).zfill(4), n), 'wb')		        
						cPickle.dump(complex_pattern, alpha, protocol = cPickle.HIGHEST_PROTOCOL)
						alpha.close()
					start_freq += freq_inc
					os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, start_freq))	    	
			    	
			    	os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, sf))
			    	os.system('oskar_settings_set  %s observation/num_channels %s' %(options.setup2, nchanns))
			    	
			    	#	header
			    	
			    	header = headers('beam_pattern_S%s_TIME_SEP_CHAN_SEP_AMP_XX.fits' %(str(0).zfill(4)))
			    	
			    	
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
					    		xx_re.append(-complex_pattern[:, 0].reshape(beam_size, beam_size))
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
					#
					""" 
					
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
	

if __name__ == '__main__':
	run_sim_test()
	
		   	
