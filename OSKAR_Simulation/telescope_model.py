#!/usr/bin/env python
"""
telescope_model.py
====================

This program models a cassagrain dish as a collection of dipoles

"""
#

import pyfits 
from optparse import OptionParser
from matplotlib import pylab as plt
import numpy as np
import scipy.interpolate as interpolate
import warnings
warnings.filterwarnings("ignore")
import os, sys
import ConfigParser
import time
#

def openFitsFile(filename):
	""" 
	Opens a FITS Filen using PyFITS
	
	Parameter
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
	
def inverse_transform_sampling(bin_edges, p_x, sampl_sz):
	""" 
	Sample uniform random numbers and transform them to your target distribution.
	
	Parameters:
	----------
	
	bin_edges: array-like
		uniform data
	
	p_x: array-like
		probability distribution	
	 
	sampl_sz: int
		number of dipoles
		
 	Returns:
 	--------
 	out array-like
 		inverse CDF 	
 	"""
 	
	cdf = np.cumsum(p_x)	
	cdf = cdf/sum(cdf)	
	cdf = cdf/max(cdf)		
	inv_cdf1 = interpolate.InterpolatedUnivariateSpline(cdf, bin_edges, k=5)	
        pdf = np.arange(0,1, 1./sampl_sz)
        
        return inv_cdf1(pdf)
#


def super_gaussian_fxn(bin_edges, fuzz_factor = 12.0, center = 0.0, width = 0.82):
	""" Generates a super-gaussian radial distribution 
	
	Parameters:
	----------
	
	fuzz_factor: float or int
		peak-shape
		
	center: float or int
		mean
	
	width: float or int
		standard deviation
	
	Returns:
 	--------
 	out array-like
 		radial distribution of dipoles 
	"""
	
	from scipy import special
						
	radd = (fuzz_factor/(2*width*special.gamma(1./width)))*(np.exp(- abs((bin_edges - center)/width)**fuzz_factor/2.))
	
	return radd
#

def grid_layout_enu(ants_x, ants_y, step_x, step_y, lonlat=None, filename=None, dish_diameter=13.5, mount="ALT-AZ", savefig=None):
	    """
	    ants_x : Antennas along x-axis (East)
	    ants_y : Antennas along y-axis (North)
	    step_x : Step size along x-axis (metres)
	    step_y : Step size along y-axis (metres)
	    lonlat : Telescope location (degrees)
	    """
	    
	    nants = ants_x * ants_y
	    enu = []
	    
	    if lonlat in [None, [], ""]:	    	
	    	lonlat = 0,0	    
	    zero_x = ants_x/2 - 0.5 if ants_x%2==0 else ants_x/2
	    zero_y = ants_y/2 - 0.5 if ants_y%2==0 else ants_y/2
	    
	    for i in xrange(ants_x):
		for j in xrange(ants_y):
		    enu.append( ( (i-zero_x)*step_x + lonlat[0], 
		                  (j-zero_y)*step_y + lonlat[1]))
	    
	    if filename:
		with open(filename, "w") as wstd:
		    wstd.write("# East North Up Dish_Diameter Station Mount\n")
		    for i,(x,y) in enumerate(enu):		        
		        wstd.write("%.6f %.6f 0.0 %.3f ST-%d %s\n"%(x, y, 
		                   dish_diameter, i, mount))		        
	    enu = np.array(enu).T                

	    if savefig:
		plt.figure(figsize=(12,12))
		plt.scatter(enu[0], enu[1])
		plt.xlabel("East [m]")
		plt.ylabel("North [m]")
		plt.grid()
		plt.savefig(savefig)
		plt.close()
		
	    return enu
#

def save_fig(prefix):
	"Save current figure in extended postscript and PNG formats."
	plt.savefig('%s.png' % prefix, format='PNG', dpi = 200)
	
	return None
	#


def feed_angle_displacement(input_directory, sub_directory, uniform_error, alpha_x, alpha_y, n_antenna):
	
	ss = np.loadtxt('%s/%s%s/layout.txt' %(input_directory, sub_directory, str(0).zfill(4)), delimiter=',')
	shape = ss.shape
	
	
	if uniform_error == 'false':
		print "\n >> Introducing random feed angle displacement per station:"		
		for iter in xrange(n_antenna):	
			feed_x =  np.random.uniform(0, alpha_x, shape[0]) 			
			feed_y =  np.random.uniform(0, alpha_y, shape[0]) 
			save_error_xy = 'feed_angle.txt'
			np.savetxt('%s' %(save_error_xy), np.array([feed_x, feed_y]).T, delimiter = ',')			
			os.system('mv -f %s  %s/%s%s ' %(save_error_xy, input_directory, sub_directory, str(iter).zfill(4)))
			#
	else:
		print "\n >> Introducing identical feed angle displacement for all stations:"	
		feed_x =  np.random.uniform(0, alpha_x, shape[0]) 			
		feed_y =  np.random.uniform(0, alpha_y, shape[0]) 
		save_error_xy = 'feed_angle.txt'	
		np.savetxt('%s' %(save_error_xy), np.array([feed_x, feed_y]).T, delimiter = ',')
		os.system('echo %s/%s* | xargs -n 1 cp %s -f' %(input_directory, sub_directory, save_error_xy))
		#
	
	return None
	#
	
		
def add_gain_phase_error(input_directory, sub_directory, uniform_error, gain_f, phase_off, time_gain_std, time_phase_std, n_antenna):
	
	ss = np.loadtxt('%s/%s%s/layout.txt' %(input_directory, sub_directory, str(0).zfill(4)), delimiter=',')
	shape = ss.shape
	
	
	if uniform_error == 'false':
		print "\n >> Introducing random gain_phase errors per station:"	
		for iter in xrange(n_antenna):	
			G_0 =  np.random.uniform(0, gain_f, shape[0]) 	
			phi_deg = np.random.uniform(0, phase_off, shape[0])
			G_std = np.random.uniform(0, time_gain_std, shape[0])
			phi_std = np.random.uniform(0, time_phase_std, shape[0])
			save_error_as = 'gain_phase.txt'
			np.savetxt('%s' %(save_error_as), np.array([G_0, phi_deg, G_std, phi_std]).T, delimiter = ',')
			os.system('cp %s -f %s/%s%s ' %( save_error_as, input_directory, sub_directory, str(iter).zfill(4)))
	
	else:
		print "\n >> Introducing identical gain_phase errors for all stations:"
		G_0 =  np.random.uniform(0, gain_f, shape[0]) 	
		phi_deg = np.random.uniform(0, phase_off, shape[0])
		G_std = np.random.uniform(0, time_gain_std, shape[0])
		phi_std = np.random.uniform(0, time_phase_std, shape[0])		
		save_error_as = 'gain_phase.txt'
		np.savetxt('%s' %(save_error_as), np.array([G_0, phi_deg, G_std, phi_std]).T, delimiter = ',')
		os.system('echo %s/%s* | xargs -n 1 cp %s -f' %(input_directory, sub_directory, save_error_as))	
	
	return None
	#
	
	
def run_dish_sim_model():
	#
	#	start-time
	start = time.time()
	startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Start at %s" % startime
	
	option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Generates telescope model to run OSKAR simulations')
	option.add_option('-f','--file1', dest='config', 
	                  default='script_config.ini', action="store", type="string",
	                  help='Enter the input parameters in sections 1 & 2 of script_config.ini')
	
	options,args = option.parse_args()
	#
	
	if os.path.exists('%s' %options.config):		
		#
		config = ConfigParser.ConfigParser()
		config.read('%s' %options.config)
	    	num_dipoles = int(config.get('antenna-layout', 'num_dipoles', 1))	    		    	
	    	dish_aperture = int(config.get('antenna-layout', 'dish_diameter', 1))
	    	num_ant = int(config.get('antenna-layout', 'num_of_ant_stations', 1))
	    	ant_directory = config.get('antenna-layout', 'telescope\input_directory', 1)
	    	ant_sub_directory  = config.get('antenna-layout', 'start_each_station_name_as', 1)	    	
	    	
      
       		fig = plt.figure()
		bins = np.arange(0.0,1,0.01) 
		pdf = super_gaussian_fxn(bin_edges = bins)
		radial_distribution = inverse_transform_sampling(bin_edges = bins, p_x = pdf, sampl_sz = num_dipoles)		
		radial_distribution = radial_distribution/radial_distribution.max() 		
		plt.hist(dish_aperture/2.0*radial_distribution,100, alpha = 0.9, weights = np.zeros_like(radial_distribution) + 1./ radial_distribution.size)
		plt.xlim(0, dish_aperture/2.0)
		plt.xlabel('Radius / m') 
		plt.ylabel('Distribution of Aperture Illumination')
		plt.title('Flat Top Gaussian Distribution')
		save_fig('radial-plot')
		plt.close()
		#
	
		fig = plt.figure()
		#
		# 	Model X & Y Orientation of Dipoles
		print " ------------------------------------------------------"
    		print "\n >> Modelling Cassegrain Dish as a collection of dipoles:"
    		print " ------------------------------------------------------"
    		
		#	Without Dish Support		
		if config.get('antenna-layout', 'add_beam_support', 1) == 'false':
			
			a = radial_distribution[np.where(radial_distribution <= 1)]
			U = np.random.uniform(0, 2*np.pi, len(a))
			x1 = dish_aperture/2*a*np.cos(U)			
			y1 = dish_aperture/2*a*np.sin(U)			
		
		#	Apply Dish support
		else:
			beam_supt = float(config.get('antenna-layout', 'add_beam_support\subflector_size', 1))	
			r1 = beam_supt/6.0
			a = radial_distribution[np.where((radial_distribution<=1) & (radial_distribution> r1))]
			U = np.random.uniform(0, 2*np.pi, len(a))	   
			x = dish_aperture/2*a*np.cos(U)			
			y = dish_aperture/2*a*np.sin(U)			
			s1 = np.where(np.logical_and(x >= - r1, x <= r1)) 
			s1 = np.take(s1, range(len(s1[0])- len(s1[0])/8))
			x1 = np.delete(x, s1)
			y1 = np.delete(y, s1)
			s2 = np.where(np.logical_and(y1 >= - r1, y1 <= r1))
			s2 = np.take(s2, range(len(s2[0]) - len(s2[0])/8))
			x1 = np.delete(x1, s2)
			y1 = np.delete(y1, s2)		    
		#
								
		plt.plot(x1, y1, 'b+') 
		plt.xlabel('X/ m', fontsize="15")
		plt.ylabel('Y/ m', fontsize="15")
		plt.title('Station Setup')
		save_fig('dish-plot_wt_%d_dipoles' %num_dipoles)
		plt.close()
	#
	
		#	Save Dipoles
	    	#
	    	
	    	for num in xrange(num_ant):
	    		         
			os.makedirs('folder/%s%s' %(ant_sub_directory, str(num).zfill(4)))
	    	os.chdir('folder/')
	    	np.savetxt('layout.txt' ,np.array([x1,y1, np.zeros(len(x1))]).T, delimiter = ',')
	    	 
	    	#	copy layout file into various stations
	    	os.system('echo %s* | xargs -n 1 cp layout.txt -f' %(ant_sub_directory))	    	
	    	os.chdir('../')
	    	
	       	if os.path.isdir("%s" %ant_directory):	
	       		os.system('rm -rf %s' %ant_directory)
	       		
       		os.makedirs('%s' %ant_directory)
	    	os.system('mv -f folder/%s* %s' %(ant_sub_directory, ant_directory))
	    	os.system('rm -rf folder')
	    	#
	    	prefix = '%s' %'layout'
	    	
	    	if config.get('antenna-layout', 'add_regular_grid', 1) == 'true':	    		 
	    		nx = int(config.get('antenna-layout', 'add_regular_grid\Xnumber_antennas', 1))
			ny = int(config.get('antenna-layout', 'add_regular_grid\Ynumber_antennas', 1))
			dx = int(config.get('antenna-layout', 'add_regular_grid\X_interval', 1))
			dy = int(config.get('antenna-layout', 'add_regular_grid\Y_interval', 1))
			lon = float(config.get('antenna-layout', 'add_regular_grid\longitude_deg', 1))
			lat = float(config.get('antenna-layout', 'add_regular_grid\latitude_deg', 1))
			if num_ant != nx*ny:
				print "InconsistentError: input number of antenna stations must be equal to antennas on regular-grid"       		
       				raise			
			if os.path.exists('%s' %prefix+".txt"):
				os.system('rm -rf %s' %prefix+".txt")			       	
		       	
			grid_layout_enu(ants_x = nx, ants_y = ny, step_x = dx, step_y = dy, lonlat = [lon, lat], dish_diameter = dish_aperture,\
					filename=prefix+".txt", savefig=prefix+"_ENU-%dx%d-grid.png" %(ny,nx))					
		
		os.system('mv -f  %s  %s' %(prefix+".txt", ant_directory))
			
	
	else:
       		print "Error: config file for input parameters does not exist"       		
       		raise
       		
	
	if config.get('distort_antenna-layout', 'include_distortion\per_station', 1) == 'true':
		#
		if config.get('distort_antenna-layout', 'add_gain_phase_error', 1) == 'true':
			xx_0 = float(config.get('distort_antenna-layout', 'sysmatic_gain_factor', 1))		
			xx_1 = float(config.get('distort_antenna-layout', 'sysmatic_phase_offset_deg', 1))		
			xx_2 = float(config.get('distort_antenna-layout', 'time_variable_gain_factor\std_deviation', 1))		
			xx_3 = float(config.get('distort_antenna-layout', 'time_variable_phase_factor\std_deviation', 1))
			
			err_type = config.get('distort_antenna-layout', 'add_uniform_gain_phase_error', 1)
					
			add_gain_phase_error(input_directory = ant_directory, sub_directory = ant_sub_directory, uniform_error = err_type, gain_f = xx_0, phase_off = xx_1 , \
			time_gain_std = xx_2, time_phase_std = xx_3, n_antenna = num_ant)
	
		if config.get('distort_antenna-layout', 'introduce_feed_displacement', 1) == 'true':
			f_x = float(config.get('distort_antenna-layout', 'feed_angle_X_deg', 1))		
			f_y = float(config.get('distort_antenna-layout', 'feed_angle_Y_deg', 1))
			
			err_type = config.get('distort_antenna-layout', 'add_uniform_feed_displacement', 1)
			
			feed_angle_displacement(input_directory = ant_directory, sub_directory = ant_sub_directory, uniform_error = err_type, alpha_x = f_x, alpha_y = f_y, n_antenna = num_ant)
	
	#	move plots to telescope directory
	os.system('mv -f  *.png  %s' %ant_directory)
		
	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0
	print "Total run time: %7.2f hours" % elasped_time
	

if __name__=="__main__":
	#
	run_dish_sim_model()


