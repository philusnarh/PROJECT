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
from pyrap.tables import table
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
	

def AddImagingColumn(msname,modeldata = None):
	
	ms = table(msname,readonly=False)
	# This function returns the names of all the columns in the MS.
	desc = ms.getcoldesc("DATA")
	
	if modeldata == None:
		print "\n >> Add Imaging data column with zeros"
		msData = ms.getcol("DATA")
		modeldata = np.zeros_like(msData, dtype=None, order='K', subok=True)
		
	if 'MODEL_DATA' in ms.colnames():
		pass
	else:		
		
		desc["name"]="MODEL_DATA"
		ms.addcols(desc)
		ms.putcol("MODEL_DATA",modeldata)
		
	if 'CORRECTED_DATA' in ms.colnames():
		pass
	else:		
		
		desc["name"]="CORRECTED_DATA"
		ms.addcols(desc)
		ms.putcol("CORRECTED_DATA",modeldata)
	ms.close()
	
	return None
#
	
def vis_noise(msname, sefd):

	tab = table(msname)
	spwtab = table("%s/SPECTRAL_WINDOW"%msname)
	freq0 = spwtab.getcol("CHAN_FREQ")[0, 0]
	wavelength = 300e+6/freq0
	bw = spwtab.getcol("CHAN_WIDTH")[0, 0]
	dt = tab.getcol("EXPOSURE", 0, 1)[0]
	dtf = (tab.getcol("TIME", tab.nrows()-1, 1)-tab.getcol("TIME", 0, 1))[0]

	tab.close()
	spwtab.close()

	noise = sefd/numpy.sqrt(abs(2*bw*dt))

	return noise


def addnoise(msname, mean = 0.0, sigma = 1.0, sefd=None, column="DATA"):
    
	tab = table(msname, readonly=False)

	if sefd:
		noise = vis_noise(msname, sefd)

	nrows = tab.nrows()
	data = tab.getcol(column)
	dshape = list(data.shape)

	rowchunk = nrows/10 if nrows > 100 else nrows

	for row0 in range(0, nrows, rowchunk):
		nr = min(rowchunk, nrows-row0)
		dshape[0] = nr

		data = noise*sigma*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape) + np.complex(mean,mean))

		print("Adding noise to %s (rows %d to %d)"%(column ,row0, row0+nr-1))

		tab.putcol(column, data, row0, nr) 
	tab.close()
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
	                  default='auto_setup.ini', action="store", type="string", 
	                  help='Enter OSKAR setup file for simulation')
        #  skar_config.ini
	
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
				visname = ss.get('interferometer', 'oskar_vis_filename', 1)
				msname = config.get('simulations', 'run_interferometer_sim\msname', 1)
				nchanns = int(ss.get('observation', 'num_channels', 1))
				phase_centre_ra_deg = float(ss.get('observation', 'phase_centre_ra_deg', 1))
				phase_centre_dec_deg = float(ss.get('observation', 'phase_centre_dec_deg', 1))
				num_time_steps = int(ss.get('observation', 'num_time_steps', 1))
				start_freq = float(ss.get('observation', 'start_frequency_hz', 1))
				sf = start_freq
				freq_inc = float(ss.get('observation', 'frequency_inc_hz', 1))
				#
				os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, 1))
				os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, 1))
				#
				pt = np.loadtxt('pointing_file.txt') 
				os.mkdir('visdir')
				#cwd = os.getcwd()                            # Current directory path 
				for iter in xrange(nchanns):
										
					print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
					print "\n >> OSKAR Interferometer simulation starts channel %d:" %(iter + 1) 
					print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
					os.system('oskar_settings_set  %s observation/start_frequency_hz %d' %(options.setup2, start_freq))
					
					for nn in xrange(num_time_steps):						
						ra_deg = pt[:,0][nn]
						dec_deg = pt[:,1][nn]
						os.system('oskar_settings_set  %s observation/phase_centre_ra_deg %f' %(options.setup2, ra_deg))
						os.system('oskar_settings_set  %s observation/phase_centre_dec_deg %f' %(options.setup2, dec_deg))
						os.system('oskar_sim_interferometer %s' %options.setup2)					
						#os.system('oskar_vis_to_ms %s -o %s_chan%s.MS' %(visname, msname, str(iter).zfill(4)))
						#msname = msname[:-3].strip().replace(" ","")
						
						os.system('cp %s  %s_nstep%s.vis' %(visname, visname[:-3].strip().replace(" ",""), str(nn).zfill(4)))
						#os.system('oskar_vis_to_ms %s -o %s_nstep%s.ms ' %(visname, msname[:-3].strip().replace(" ",""), str(nn).zfill(4))) 
						os.system('mv %s_nstep%s.vis %s' %(visname[:-3].strip().replace(" ",""), str(nn).zfill(4), 'visdir'))
						#os.system('mv *.ms %s' %('msdir'))
						#os.system('cp %s  %s_chan%s.vis' %(visname, msname, str(iter).zfill(4)))
					#os.chdir("%s" %'msdir')
					os.chdir("%s" %'visdir')
					#ms_sort = sorted(glob.glob('./*.ms'))
					os.system('oskar_vis_add *.vis -o combined_vis_chan%s.vis' %(str(iter).zfill(4)))
					
					#os.system('pyxis %s output=concat_chan%s.MS' %(ms_sort, str(iter).zfill(4)))
					os.chdir("../")	
					start_freq += freq_inc
					#start_freq += freq_inc
					#oskar_vis_add *.vis
					
				"""		
				integration_time = int(ss.get('interferometer', 'time_average_sec', 1))
				beamtime_slots = int(ss.get('observation', 'num_time_steps', 1))
				beamform_len = ss.get('observation', 'length', 1)
				st = config.get('simulations', 'run_interferometer_sim\Synthesis_steps', 1) 
				
				if msname.endswith(('.ms', '.MS')):
					msname = msname[:-3].strip().replace(" ","")			
				x = time.strptime(st.split(',')[0],'%H:%M:%S.%f')
							
				date_time = datetime.timedelta(hours=x.tm_hour,minutes=x.tm_min,seconds=x.tm_sec).total_seconds()			
				os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, int(date_time/integration_time) + 1))
				os.system('oskar_settings_set  %s observation/length %s' %(options.setup2, st))
				nchanns = int(ss.get('observation', 'num_channels', 1))
				start_freq = float(ss.get('observation', 'start_frequency_hz', 1))
				sf = start_freq
				freq_inc = float(ss.get('observation', 'frequency_inc_hz', 1))
				
				
				if nchanns > 1000:
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
				
					print "\n >> Checking for Imaging Columns"
					AddImagingColumn(msname = '%s.MS' %msname,modeldata = None)
					os.system('mv -f *.MS %s' %filename )
					os.system('rm -f %s' %visname)
				"""  	
					
				        #os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, beamtime_slots))
				        #os.system('oskar_settings_set  %s observation/length %s' %(options.setup2, beamform_len))
		        
				os.system('oskar_settings_set  %s observation/start_frequency_hz %f' %(options.setup2, sf))
				os.system('oskar_settings_set  %s observation/num_channels %d' %(options.setup2, nchanns))
				os.system('oskar_settings_set  %s observation/num_time_steps %d' %(options.setup2, num_time_steps))
		        
		        else:
				
				print "\n >> Error: OSKAR setup file does not exist"       		
	       			raise	
	       			
       			# ++++++++++++ ADD NOISE  +++++++++++++++
	       		
			if config.get('add_noise_to_vis', 'add_noise', 1) == 'yes':
			
				print "\n >> adding noise to visibilities"
				
				t_sys = int(config.get('add_noise_to_vis', 'system_temperature\Kelvin_units', 1))
				ant_eff = float(config.get('add_noise_to_vis', 'antenna_efficiency', 1))
				
				aper_size = int(config.get('antenna-layout', 'dish_diameter', 1))
				data_col = config.get('antenna-layout', 'msname\data_column', 1)
				mu_noise = float(config.get('add_noise_to_vis', 'mean_noise\distribution', 1))
				scale_noise = float(config.get('add_noise_to_vis', 'std_deviation_noise\distribution', 1))
				ant_area = np.pi*aper_size**2
				k_bzm = 1.38e-38		# Boltzmann Constant
				sefd = 	2*k_bzm*t_sys/(ant_eff * ant_area)
					
				if config.get('add_noise_to_vis', 'overwrite_msname', 1) == 'no':
					print "\n >> creating new msname  "
					os.system('mv %s.MS %s_%s.MS' %(msname, msname, 'noise_data'))
					
				new_msname = glob.glob('*noise_data.MS')					
					
									
				addnoise(msname = new_msname, mean = mu_noise, sigma = scale_noise, sefd = sefd, column = data_col)
						     					
				
		
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
	
		   	
