#
"""
	 This script generates 2x2 Jones & 4x4 Mueller Beam Plots from OSKAR Beampattern Simulations. It also 
	 produces Beam Errors between 2 Beampatterns.
	  
	 To run script use:
	 eg: python beampattern_plot.py -f <configFile.ini>
	 where, <configFile.ini> is the input parameters in .ini file
  """
#

import pyfits 
import pywcsgrid2
import pywcs
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
from matplotlib import pylab as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import os, sys
import time
import ConfigParser
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
        	print "Error: cannot open %s"%filename
        	raise
#


def setup_axes(fig, header, nrows, ncols, i):
	    """ Create an inset axes with a given width and height """
	    gh = pywcsgrid2.GridHelper(wcs=header)
	    gh.locator_params(nbins=3)
	    g = axes_grid.ImageGrid(fig, 111,
		                    nrows_ncols=(nrows, ncols),
		                    ngrids=None,
		                    direction='row',
		                    axes_pad=0.02, add_all=True,
		                    share_all=True, aspect=True,
		                    label_mode='L', cbar_mode=None,
		                    cbar_location='right',
		                    axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh))) 
	    # make colorbar
	    
	    cax = []
	    for iter in xrange(i):		
		ax = inset_axes(g[iter],
		             width="45%", 
		             height="5%", 
		             loc=9)
		cax.append(ax)
		
	    return g, cax
	    

#
def Mueller_image(fig, header, nrows, ncols, i , Mueller_matrix, figname):
	    """ 4x4 Mueller matrix plot """
	    
	    g2, cax2 = setup_axes(fig, header, nrows, ncols, i)
	    
	    
	    cmap = plt.get_cmap('jet')	
	    images = []
	    
	    for i, ax in enumerate(g2):
		channel = Mueller_matrix[i]			
		im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos")
		if i in [0, 5, 10, 15]:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel),  np.max(channel)], \
		   orientation = "horizontal") 
		    
		   cax2[i].xaxis.set_ticks_position("bottom") 
		    
		else:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel), np.max(channel)],\
		   orientation = "horizontal")
		   
		   cax2[i].xaxis.set_ticks_position("bottom")  
		images.append(im)
	    plt.savefig(figname)
	    
	    return None 
#  

def Jones_image(fig, header, nrows, ncols, i , Jones_matrix, figname):
	    """ 2x2 Jones matrix plot """
	    g1, cax1 = setup_axes(fig, header, nrows, ncols, i)		
	    
	    cmap = plt.get_cmap('jet')		
	    images = []
	    
	    for i, ax in enumerate(g1):
		channel =  Jones_matrix[i]		
		im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos")
		plt.colorbar(im, cax=cax1[i], format = '%.0e', ticks = [np.min(channel),  np.max(channel)],\
		orientation = "horizontal")
		
		cax1[i].xaxis.set_ticks_position("bottom")
		images.append(im)
	    plt.savefig(figname)
    
    	    return None


def figur(n1, n2, n3, n4):
	    fx = plt.figure(n1, figsize=(n2, n3), dpi= n4)
	    
	    return fx
	    
#    
def headers(Fits_filename):
	    hd = pyfits.getheader(Fits_filename)
	    
	    return hd
	    
def jones_filename(beam_dir1, beam1_start_filename, station_num, chan_num):

	jns = ['XX', 'XY', 'YX', 'YY'] 							
	ty = ['re', 'im']
	k = ['Eampl.fits', 'Ephase.fits']
	
	for ii in xrange(len(ty)):
		m = []              	    
		for jj in xrange(len(jns)):			
			beam_file = '%s/%s_jones_S%s_%s_%s.fits' %(beam_dir1, beam1_start_filename, \
			str(station_num).zfill(4), jns[jj], ty[ii]) 
			
			hdu = openFitsFile('%s' %(beam_file))
			beam_data = hdu[0].data[chan_num]
			m.append(beam_data)
		beam_header = hdu[0].header
		pyfits.writeto('%s' %(k[ii]), np.array(m), header=beam_header, clobber=True)
	
	return None

#
def run_test():
	
	start = time.time()
        startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Start at %s" % startime
        
        #	Accept standard options by all instances
        
        option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Generates Jones 2 X 2 & 4 X 4 Mueller Matrices Plots from OSKAR Beampattern simulations')
	option.add_option('-f','--file1', dest='configFile', 
	                  default='beamsetup.ini', action="store", type="string",
	                  help='Enter the input parameters')
	
	options,args = option.parse_args()
	#
		
	if os.path.exists('%s' %options.configFile):
		config = ConfigParser.ConfigParser()
		config.read('%s' %options.configFile)
		beam_dir1 = config.get('plot_beam_image', 'beam_1\input_directoy', 1)
		beam_image1 = config.get('plot_beam_image', 'beam_1\image_name', 1)
		
		if config.get('plot_beam_image', 'all_stations', 1) == 'yes':
			total_num_stations = int(config.get('antenna-layout', 'num_of_ant_stations', 1))
			station_num = 0
			
		else:
		
			station_num = int(config.get('plot_beam_image', 'station_number', 1))
			total_num_stations = 1
			
		chan_num = int(config.get('plot_beam_image', 'channel_number', 1))
		beam1_start_filename = config.get('plot_beam_image', 'beam_1_start_name', 1)
		save_fits = config.get('plot_beam_image', 'save_mueller_image_as_FITS', 1)
		
		
		jj = 0
		
		while jj < total_num_stations:
		
			#  ++++++++++++++++++	FIRST PART	++++++++++++++++++++++++++++++++++++
		
			#	Extract Jones terms in amplitude & phase
			jones_filename(beam_dir1, beam1_start_filename, station_num, chan_num)		
			beam_amp = 'Eampl.fits' 			# amplitude
			hdu = openFitsFile('%s' %(beam_amp))
			beam_data_amp = hdu[0].data
			beam_ph = 'Ephase.fits'				# phase 
			hdu = openFitsFile('%s' %(beam_ph))
			beam_data_ph = hdu[0].data
			beam_header = headers('%s' %beam_ph) 
			
		
			#	Plot Jones terms
			fv1 = figur(1, 12, 12, 70)
	    		Jones_image(fv1, beam_header, 2, 2, 4 ,  beam_data_amp, '%s_S%s_2x2_Jones1_Real_Image.png'\
	    		%(beam_image1, str(station_num)))
	    		 		    	
	    		fv2 = figur(2, 12, 12, 70)
	    		Jones_image(fv2, beam_header, 2, 2, 4 , beam_data_ph, '%s_S%s_2x2_Jones1_Imag_Image.png'\
	    		%(beam_image1, str(station_num)))
		
			beam_file = '%s/%s_S%s_mueller_beam.fits' %(beam_dir1, beam1_start_filename, \
			str(station_num).zfill(4)) 
			
			hdu1 = openFitsFile('%s' %(beam_file)) 			
			mueller_beam1 = hdu1[0].data[chan_num]      
			fv1 = figur(3, 12, 12, 80)
		
			if save_fits == 'yes':
			
				pyfits.writeto('%s_S%s_mueller_4_x_4_beam1.fits' %(beam_image1, str(station_num).zfill(4)), \
				mueller_beam1, header = beam_header, clobber=True)
			
		    	Mueller_image(fv1, beam_header, 4, 4, 16 , mueller_beam1, '%s_S%s_Meuller_4_X_4_Images1.png' \
		    	%(beam_image1, str(station_num).zfill(4)))
		    	
		    	
		    	# +++++++++++++++++++++++	SECOND PART 	+++++++++++++++++++++++++++++++++++++++
		    	
		    	if config.get('plot_beam_image', 'compare_beam', 1) == 'true':
		    	
			    	beam_dir2 = config.get('plot_beam_image', 'beam_2\input_directoy', 1)
			    	beam2_start_filename = config.get('plot_beam_image', 'beam_2_start_name', 1)
			    	beam_diffname = config.get('plot_beam_image', 'beam_error\image_nam', 1)
			    	
			    	#	Extract Jones terms in amplitude & phase
				jones_filename(beam_dir2, beam2_start_filename, station_num, chan_num)		
				beam_amp = 'Eampl.fits' 			# amplitude
				hdu = openFitsFile('%s' %(beam_amp))
				beam_data_amp = hdu[0].data
				beam_ph = 'Ephase.fits'				# phase 
				hdu = openFitsFile('%s' %(beam_ph))
				beam_data_ph = hdu[0].data
				beam_header = headers('%s' %beam_ph)		
		
				#	Plot Jones terms
				fv1 = figur(1, 12, 12, 70)
		    		Jones_image(fv1, beam_header, 2, 2, 4 ,  beam_data_amp, '%s_S%s_2x2_Jones2_Real_Images.png'\
		    		%(beam_image2, str(station_num).zfill(4)))
		    				    	
		    		fv2 = figur(2, 12, 12, 70)
		    		Jones_image(fv2, beam_header, 2, 2, 4 , beam_data_ph, '%s_S%s_2x2_Jones2_Imag_Images2.png' %(beam_image2, \
		    		str(station_num).zfill(4)))
			 	
			 	#	Plot Jones terms
				beam_file = '%s/%s_S%s_mueller_beam.fits' %(beam_dir2, beam2_start_filename, str(station_num).zfill(4)) 
				hdu_1 = openFitsFile('%s' %(beam_file))  
				mueller_beam2 = hdu[0].data[chan_num]      
				fv1 = figur(3, 12, 12, 80)
			
				if save_fits == 'yes':
			
					pyfits.writeto('%s_S%s_mueller_4_x_4_beam2.fits' %(beam2_start_filename, str(station_num).zfill(4)),\
					 mueller_beam2, header = beam_header, clobber=True)
					 
			    	Mueller_image(fv1, beam_header, 4, 4, 16 , mueller_beam2, '%s_S%s_4x4_Meuller_Images.png' \
			    	%(beam_image2, str(station_num).zfill(4)))
			    	
			    	diff_beam = mueller_beam1 - mueller_beam2
			    	
			    	#	Plot Beam Errors
			    	
			    	if save_fits == 'yes':
			    	
			    		pyfits.writeto('%s_S%s_mueller_4_x_4_beam.fits' %(beam_diffname, \
			    		str(station_num).zfill(4)), diff_beam, header = beam_header, clobber=True)
			    		
			    	fv1 = figur(4, 12, 12, 80)
			    	Mueller_image(fv1, beam_header, 4, 4, 16 , diff_beam, '%s_S%s_4x4_Meuller_Images.png' \
			    	%(beam_image2, str(station_num).zfill(4)))
			    	
		    	
		    	os.system('rm -f Eampl.fits Ephase.fits')
		    	
		    	#	Creating Output directory
		    			
			filename = config.get('plot_beam_image', 'output_directory', 1)
			if os.path.isdir("%s" %filename):		
		       		#os.popen('rm -rf %s' %filename)
		       		os.system('mv -f *png %s '  %(filename))
		       		
		       		if os.path.exists('%s_S%s_mueller_4_x_4_beam1.fits' %(beam_image1, str(station_num).zfill(4))):
		       			os.system('mv -f *fits %s '  %(filename))
		
			else:		        
				os.makedirs('%s' %filename)
				os.system('mv -f *png %s '  %(filename))
				
		       		if os.path.exists('%s_S%s_mueller_4_x_4_beam1.fits' %(beam_image1, str(station_num).zfill(4))):
		       			os.system('mv -f *fits %s '  %(filename))
       			jj += 1
       			station_num += 1
       			
		
    	else:
				
		print "Error: OSKAR setup file does not exist"       		
		raise		
    	#        
     
        stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0		
        print "Total run time: %7.2f hours" % elasped_time


if __name__ == '__main__':
	run_test()
	
	

