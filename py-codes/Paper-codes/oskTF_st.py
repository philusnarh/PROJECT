import pyfits as pf
import numpy as np
import os, sys
from matplotlib import pylab as plt
import healpy as hp
from itertools import permutations, chain
from scipy.signal import savgol_filter 		
import warnings
warnings.filterwarnings("ignore")
#import os, sys
import ConfigParser
import time
from optparse import OptionParser

#
def figur(fig_num, w_size, h_size, DPI):
    fx = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
    
    return fx
    
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
def f_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U' , 'V']
	#k1 =  sorted(k)
	for i in xrange(conv_map.shape[0] ):	
		
		if i in xrange(1):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, cbar=True,\
			unit = 'Jy/pixel', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (1,4 , i + 1))		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, cbar=True, \
			unit = 'Jy/pixel',  margins=(0.01, 0.01, 0.0, 0.02),\
			notext = True, sub = (1, 4, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize=20, fontweight="bold") 
	#plt.show()
	fig.savefig(figname)
	plt.close(fig_num)
	
    	return None
    	
    	
#
def conv_healpix_plott_all(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['$Measured \, Stokes \, I^{MB_{T}}$', 
	     '$Measured \, Stokes \, Q^{MB_{T}}$', 
	     '$Measured \, Stokes \, U^{MB_{T}}$', 
	     '$Measured \, Stokes \, I^{MB_{XY}}$',
	      '$Measured \, Stokes \, Q^{MB_{XY}}$', 
	      '$Measured \, Stokes \, U^{MB_{XY}}$', 
	     '$ Error \, in \, Stokes \, I$', 
	     '$ Error \, in \, Stokes \, Q$', 
	     '$ Error \, in \, Stokes \, U$' ]
	#k1 =  sorted(k)
	
	for i in xrange(conv_map.shape[0]):	
		#print conv_map.shape[0]
		if i in xrange(6):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], 
			                     norm = 'log', xsize = 3200, cbar =True,
			                     unit = 'Jy/beam', margins = (0.01, 0.01, 0.0, 0.02), 
			                     notext = True, sub = (3,3 , i + 1))
			                     #		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear',
			                     xsize = 3200, cbar = True, unit = 'Jy/beam',
			                     margins = (0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize = 20, fontweight="bold") 
	fig.savefig(figname)	
	plt.close(fig_num)
	
    	return None
#



#
def conv_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U', 'V']*4
	k1 =  sorted(k)
	j = 0
	for i in xrange(conv_map[0].shape[0]):	         # [0, 1, 2, 4, 5, 6, 8, 9, 10]:
		if i in [0, 4, 8]:   		#[0]
		
			j = j + 1
			hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, \
			cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
			plt.title(r'$ %s \longrightarrow %s $' %( k[i],k1[i]), fontsize=20, fontweight="bold")
		elif i in [3, 7, 11, 12, 13, 14, 15]:
			continue
		
		else:
			j = j + 1
			hp.visufunc.mollview(conv_map[0, i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, \
			cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, j )) 
			plt.title(r'$ %s \longrightarrow %s $' %(k[i], k1[i]), fontsize=20, fontweight="bold")
	fig.savefig(figname)
	plt.close(fig_num)
	
    	return None
    	
    	
    	
#
def error_conv_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U']*3
	k1 =  sorted(k)
	j = 0
	for i in xrange(conv_map.shape[0]):	       
		'''if i in [0, 4, 8]:   		
		
			j = j + 1
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, \
			cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
			plt.title(r'$ E_{re}^{%s \longrightarrow %s}$' %( k[i],k1[i]), fontsize=20, fontweight="bold")
			'''
		#elif i in [3, 7, 11, 12, 13, 14, 15]:
			#continue
		
		#else:
		j = j + 1
#		if j in [1, 5, 9]:
#			ss = '$ E_{ab}^{%s \longrightarrow %s}$' %(k[i], k1[i])
#		else:
#			ss = '$ E_{ab}^{%s \longrightarrow %s}$' %(k[i], k1[i])
		hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, \
		cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, j )) 
		plt.title(r'$ E_{re}^{%s \longrightarrow %s}$' %(k[i], k1[i]), fontsize=20, fontweight="bold")
	fig.savefig(figname)
	plt.close(fig_num)
	
    	return None
    	
    	
#
def run_test():
	#       start-time
        start = time.time()
        startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Start at %s" % startime

        option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
        option.set_description('Plots convmaps for corrupted and uncorrupted foregrounds simulations')
        option.add_option('-f','--file1', dest='setup',
                          default='newsetup_XY.ini', action="store", type="string",
                          help='Enter config file for input parameters ')
        
        options,args = option.parse_args()
        #
        
        config = ConfigParser.ConfigParser()
        config.read('%s' %options.setup)
        true_skymap= config.get('parameters', 'healpix_map', 1)
        chan_map= int(config.get('parameters', 'healpix_map/channel', 1))
        uncorrupt_skymap = config.get('parameters', 'uncorrupted_healpix_map', 1)        
        corrupt_skymap = config.get('parameters', 'corrupted_healpix_map', 1)
        strt_sav_file = config.get('parameters', 'start_save_filenames_as', 1)
        result_directory = config.get('parameters', 'output_directory', 1)
        
	hdu_gridd = openFitsFile('%s' %(true_skymap))
	dat = hdu_gridd[0].data
	
	# dat.astype()
	#raise
	#a = a.astype(numpy.float32).
	#conv_healpix_plott(conv_map = v, figname = 'tt.png', fig_num = 1, w_size = 18, h_size = 10, DPI = 500)
	
	for iter in [chan_map]: #xrange(1):
		print 'Started plotting chan:++++	', iter
		
		hdu_gridd = openFitsFile('%s' %(corrupt_skymap))
		datt = hdu_gridd[0].data
		v1 = datt.copy()	
		hdu_gridd = openFitsFile('%s' %(uncorrupt_skymap))
		datt = hdu_gridd[0].data
		v = datt.copy()		
		
		
		conv_healpix_plott(conv_map = v, figname = '%s_Tchan_%d.png' %(strt_sav_file, iter), fig_num = 1 ,
		                   w_size = 18, h_size = 10, DPI = 500)
		conv_healpix_plott(conv_map = v1, figname = '%s_Dchan_%d.png' %(strt_sav_file, iter), fig_num = 2 , 
		                   w_size = 18, h_size = 10, DPI = 500)
		f_healpix_plott(conv_map = dat[iter], figname = '%s_fchan_%d' %(strt_sav_file, iter), fig_num = 3 , 
		                   w_size = 18, h_size = 3, DPI = 500)
		#plt.show()                  
                #break
		                   
		print '\n >>> Plotting Absolute difference btn the conv maps'
		#
		# Relative Error
		conv_errors = [] 
		from scipy.stats import sem
		for jj in xrange(16):
			
			if jj in [3, 7, 11, 12, 13, 14, 15]:
				continue
			else:
			
				m1 = v1[0, jj, ...]**2 #abs(v1[0, jj, ...]) 
				m2 = v[0, jj, ...]**2 #abs(v[0, jj, ...])
				cc = (m1 - m2)
				#xg = cc/cc.mean()	#sem(cc/cc.mean(), axis=None, ddof=0)
				#print xg
				
				#cc = (m1 - m2)**2
				conv_errors.append(cc)
				print jj
		
		np.save('fdhase_error_%d' %iter, np.array(conv_errors)) 
		error_conv_healpix_plott(conv_map = np.array(conv_errors), figname = '%s_Echan_%d.png' %(strt_sav_file, iter),
		 			fig_num = 15, w_size = 18, h_size = 10, DPI = 500)
		print '\n >>> done !!!!!!!!!!'
		
		#raise
#		
		T_arr = np.array([v[0, 0, ...], v[0, 1, ...], v[0, 2, ...], 
				 v[0, 4, ...] , v[0, 5, ...], v[0, 6, ...],
				 v[0, 8, ...] , v[0, 9, ...] , v[0, 10, ...]])
		np.save('Tconv_map', T_arr)
		lst = []
		ti = v[0, 0, ...] + v[0, 1, ...] + v[0, 2, ...] 
		tq = v[0, 4, ...] + v[0, 5, ...] + v[0, 6, ...] 
		tu = v[0, 8, ...] + v[0, 9, ...] + v[0, 10, ...] 
		
		di = v1[0, 0, ...] + v1[0, 1, ...] + v1[0, 2, ...] 
		dq = v1[0, 4, ...] + v1[0, 5, ...] + v1[0, 6, ...] 
		du = v1[0, 8, ...] + v1[0, 9, ...] + v1[0, 10, ...]
		
#		
		lst.append(ti)
		lst.append(tq)
		lst.append(tu)
		lst.append(di)
		lst.append(dq)
		lst.append(du)
		lst.append(ti - di)
		lst.append(tq - dq)
		lst.append(tu - du)
		#
		
		
		lst = np.array(lst)
		
		
		conv_healpix_plott_all(conv_map = lst, figname = '%s_All_chan_%d.png' %(strt_sav_file, iter), fig_num = 4 ,
		                       w_size = 18, h_size = 10, DPI = 500)
		
		
		#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Stokes I
		
		nm = int(420)	#720
		u1 = hp.sphtfunc.anafast(dat[iter,0,...], lmax = nm) 		# True Sky
		m_IT = v[0, 0, ...] + savgol_filter(v[0, 1, ...] + v[0, 2, ...], 21, 3, deriv = 1) 
		m_ID = v1[0, 0, ...] + savgol_filter(v1[0, 1, ...] + v1[0, 2, ...], 21, 3, deriv = 1)
		u = hp.sphtfunc.anafast(m_IT,lmax = nm) 				# Conv I^T
		ud = hp.sphtfunc.anafast(m_ID, lmax = nm)				# Conv. I^D
		
		#Leakage term  
		#ul = hp.sphtfunc.anafast(m_IT - m_ID, lmax = nm) 
		#
		lmax = np.arange(len(u1))/3.0
		#np.save('lmax_%d' %iter, np.array(lmax)) 
		fig = figur(fig_num = 5 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $I$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $I$ with ${MB_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$I$ with ${MB_{XY}}$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $I$')
		#		
		ax1.set_xlabel(r" $l$ (multipole)", fontsize="20")
		ax1.set_ylabel(r" Intensity $[Jy]$", fontsize="20")
				
		ax2 = ax1.twinx()		
		
		ax2.semilogy(lmax,u/u1, '-y', label = '${MB_{T}}$ ')
		print (u/u1).max()
		
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${MB_{XY}}$ ') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize="20")		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='upper right', framealpha = 0.5)		

		fig.savefig('iipwchan_%d.png' %iter)
		fig.clear()
				
		
		# 
		#	Dictionary
		#		I
		#if config.get('parameters', 'sav_percentage_leakage', 1) == 'yes' :
		
		#	stokes_leakage = []
		#	stokes_leakage.append(savgol_filter(ul/u1 * 100, 51, 3))
			
		#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Stokes Q
		
		u1 = hp.sphtfunc.anafast(dat[iter, 1,...], lmax = nm) 
		m_QT =	v[0, 5, ...] + savgol_filter(v[0, 6, ...], 21, 3, deriv = 2) 		
		m_QD =	v1[0, 5, ...] + savgol_filter(v1[0, 6, ...], 21, 3, deriv = 2) 
		u = hp.sphtfunc.anafast(m_QT, lmax = nm)
		ud = hp.sphtfunc.anafast(m_QD, lmax = nm) 
		
		#Leakage term     
		ul = hp.sphtfunc.anafast(m_QT - m_QD, lmax = nm) 
		# 
		
		fig = figur(fig_num = 6 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $Q$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $Q$ with ${MB_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$Q$ with ${MB_{XY}}$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in  $Q$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",fontsize="20")
		ax1.set_ylabel(r" Intensity $[Jy]$",fontsize="20")
				
		ax2 = ax1.twinx()		
		ax2.semilogy(lmax,u/u1, '-y', label = '${MB_{T}}$')
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${MB_{XY}}$') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize="20")		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='upper right', framealpha = 0.5)		

		fig.savefig('qqpwchan_%d.png' %iter)
		fig.clear()
		
		#
		
		#if config.get('parameters', 'sav_percentage_leakage', 1) == 'yes' :		
			
		#	stokes_leakage.append(savgol_filter(ul/u1 * 100, 21, 3))
		#
		# Stokes U
		
		u1 = hp.sphtfunc.anafast(dat[iter, 2,...], lmax = nm) 
		m_UT =	v[0, 10, ...] + savgol_filter(v[0, 9, ...], 21, 3, deriv = 2) 		
		m_UD =	v1[0, 10, ...] + savgol_filter(v1[0, 9, ...], 21, 3, deriv = 2) 
		u = hp.sphtfunc.anafast(m_UT, lmax = nm) 
		ud = hp.sphtfunc.anafast(m_UD, lmax = nm)
		
		#Leakage term  
		ul = hp.sphtfunc.anafast(m_UT - m_UD, lmax = nm)
		#
		fig = figur(fig_num = 7 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $U$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $U$ with ${MB_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$U$ with ${MB_{XY}}$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $U$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",fontsize="20")
		ax1.set_ylabel(r" Intensity $[Jy]$",fontsize="20")
				
		ax2 = ax1.twinx()		
		ax2.semilogy(lmax,u/u1, '-y', label = '${MB_{T}}$ ')
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${MB_{XY}}$') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize="20")		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='upper right', framealpha = 0.5)		

		fig.savefig('uupwchan_%d.png' %iter)
		fig.clear()
		
		
		#
		
#		if config.get('parameters', 'sav_percentage_leakage', 1) == 'yes' :		
#			
#			stokes_leakage.append(savgol_filter(ul/u1 * 100, 21, 3))
#			d1 =  { 'l_max' : np.array(lmax), 'IQU' : np.array(stokes_leakage)} 
#			np.savez('%s_fsky_percentage_leakage' %(strt_sav_file), **d1) 
			
	if os.path.isdir("%s" %(result_directory)):
		os.popen('rm -r %s' %result_directory)
	os.makedirs('%s' %(result_directory))
	os.system('mv -f *.%s %s' %('png', result_directory))


	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Stop at %s" % stoptime
	end = time.time()
	elasped_time = (end - start)/3600.0
	print "Total run time: %7.2f hours" % elasped_time
        
	print 'done wai !!!!!!!!!!!!'
if __name__ == '__main__':
	run_test()

