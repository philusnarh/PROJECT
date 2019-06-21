#import pyfits as pf
import astropy.io.fits as pf
import numpy as np
import os, sys
from matplotlib import pylab as plt
from matplotlib import cm
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
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
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
	k = ['I', 'Q', 'U', 'V' ]
	#k1 =  sorted(k)
	for i in xrange(conv_map.shape[0]):	   # (conv_map.shape[0] - 1)
		
		if i in xrange(1):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, cbar=True,\
			unit = 'k', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (2,2 , i + 1))		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, cbar=True, \
			unit = 'k',  margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (2, 2, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize=20, fontweight="bold") 
	fig.savefig(figname)
	plt.close(fig_num)
	
    	return None
    	
    	
#
#
def conv_healpix_plott(conv_map, Tmap, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U', 'V']*4
	k1 =  sorted(k)
	j = 0
	for i in xrange(conv_map[0].shape[0]):	         # [0, 1, 2, 4, 5, 6, 8, 9, 10]:
		if i in [0, 4, 8]:   		#[0]
		
			j = j + 1
			if i == 0:
				hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, min=Tmap[0].min(),\
				max=Tmap[i].max(), cbar=True, unit = 'K/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
				plt.title(r'$ %s \longrightarrow %s $' %( k[i],k1[i]), fontsize=20, fontweight="bold")
			else:
				hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, \
				cbar=True, unit = 'K/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
				plt.title(r'$ %s \longrightarrow %s $' %( k[i],k1[i]), fontsize=20, fontweight="bold")
		elif i in [3, 7, 11, 12, 13, 14, 15]:
			continue
		
		else:
			j = j + 1
			if i in [5,10]:
				if i == 5: ck = 1
				else: ck = 2
				hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'linear', xsize = 3200, min=Tmap[ck].min(),\
				max=Tmap[ck].max(), cbar=True, unit = 'K/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
				plt.title(r'$ %s \longrightarrow %s $' %( k[i],k1[i]), fontsize=20, fontweight="bold")
			else:
				hp.visufunc.mollview(conv_map[0, i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, \
				cbar=True, unit = 'k/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, j )) 
				plt.title(r'$ %s \longrightarrow %s $' %(k[i], k1[i]), fontsize=20, fontweight="bold")
	fig.savefig(figname)
	plt.close(fig_num)
	
    	return None
    	
def conv_healpix_plott_all(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['$Measured \, Stokes \, I^{T}$', 
	     '$Measured \, Stokes \, Q^{T}$', 
	     '$Measured \, Stokes \, U^{T}$', 
	     '$Measured \, Stokes \, I_{D}^{GP}$',
	      '$Measured \, Stokes \, Q_{D}^{GP}$', 
	      '$Measured \, Stokes \, U_{D}^{GP}$', 
	     '$ Error \, in \, Stokes \, I$', 
	     '$ Error \, in \, Stokes \, Q$', 
	     '$ Error \, in \, Stokes \, U$' ]
	#k1 =  sorted(k)
	
	for i in xrange(conv_map.shape[0]):	
		#print conv_map.shape[0]
		if i in xrange(6):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], 
			                     norm = 'log', xsize = 3200, cbar =True,
			                     unit = 'K/beam', margins = (0.01, 0.01, 0.0, 0.02), 
			                     notext = True, sub = (3,3 , i + 1))
			                     #		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear',
			                     xsize = 3200, cbar = True, unit = 'K/beam',
			                     margins = (0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize = 20, fontweight="bold") 
	fig.savefig(figname)	
	plt.close(fig_num)
	
    	return None
#



##
#def conv_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI):
#	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
#	k = ['I', 'Q', 'U', 'V']*4
#	k1 =  sorted(k)
#	j = 0
#	for i in xrange(conv_map[0].shape[0]):	         # [0, 1, 2, 4, 5, 6, 8, 9, 10]:
#		if i in [0, 4, 8]:   		#[0]
#		
#			j = j + 1
#			hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, \
#			cbar=True, unit = 'K/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3,3 , j))
#			plt.title(r'$ %s \longrightarrow %s $' %( k[i],k1[i]), fontsize=20, fontweight="bold")
#		elif i in [3, 7, 11, 12, 13, 14, 15]:
#			continue
#		
#		else:
#			j = j + 1
#			hp.visufunc.mollview(conv_map[0, i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, \
#			cbar=True, unit = 'k/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, j )) 
#			plt.title(r'$ %s \longrightarrow %s $' %(k[i], k1[i]), fontsize=20, fontweight="bold")
#	fig.savefig(figname)
#	plt.close(fig_num)
#	
#    	return None
#    	
#    	
    	
#
def error_conv_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI, mask=None):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U']*3
	k1 =  sorted(k)
	j = 0
#	mask_indx = []
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
		if mask is not None:
			print 'gakpo !!!'
			print 'yeh'
##			mp = mask_healmap(NSIDE = hp.npix2nside(len(conv_map[i,...])), data=conv_map[i,...])
#			ind1 = np.where(abs(conv_map[i,...]) >= mask[i])
##			ind1 = np.where(abs(conv_map[i,...]) >= mask)
#			mask_indx.append(ind1)
#			conv_map[i,...][ind1] = np.nan
#			
			conv_map[i,...][mask[0,0]] = np.nan     # mask[i,0]
			print conv_map[i,...]
			
			
		hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, \
		cbar=True, unit = 'K/beam', margins=(0.01, 0.01, 0.0, 0.02), notext = True, sub = (3, 3, j)) #, min=0, max=1000) 
		plt.title(r'$ E_{re}^{%s \longrightarrow %s}$' %(k[i], k1[i]), fontsize=20, fontweight="bold")
#	np.save('mask', mask_indx)
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
                          default='set00.ini', action="store", type="string",
                          help='Enter config file for input parameters ')
        
        options,args = option.parse_args()
        #
        
        config = ConfigParser.ConfigParser()
        config.read('%s' %options.setup)
        true_skymap= config.get('parameters', 'healpix_map', 1)
        chan_map= int(config.get('parameters', 'healpix_map/channel', 1))
        uncorrupt_skymap = config.get('parameters', 'uncorrupted_healpix_map', 1)
        corrupt_skymap = config.get('parameters', 'corrupted_healpix_map', 1)
        xycorrupt_skymap = config.get('parameters', 'xy_corrupted_healpix_map', 1)
        vcorrupt_skymap = config.get('parameters', 'v_corrupted_healpix_map', 1)
        strt_sav_file = config.get('parameters', 'start_save_filenames_as', 1)
        result_directory = config.get('parameters', 'output_directory', 1)
        tp = 1#  e5
	hdu_gridd = openFitsFile('%s' %(true_skymap))
	dat = hdu_gridd[0].data*tp
	fsz = 15
	
	# dat.astype()
	#raise
	#a = a.astype(numpy.float32).
	#conv_healpix_plott(conv_map = v, figname = 'tt.png', fig_num = 1, w_size = 18, h_size = 10, DPI = 500)
	
	for iter in [chan_map]:
		print 'Started plotting chan:++++ ', iter
		
		hdu_gridd = openFitsFile('%s' %(corrupt_skymap))
		datt = hdu_gridd[0].data*tp
		v1 = datt.copy()	
		hdu_gridd = openFitsFile('%s' %(uncorrupt_skymap))
		datt = hdu_gridd[0].data*tp
		v = datt.copy()
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# plotting polarisation leakage
		hdu_gridd = openFitsFile('%s' %(xycorrupt_skymap))
		datt = hdu_gridd[0].data*tp
		xy = datt.copy()	    # TRUE
		
		hdu_gridd = openFitsFile('%s' %(vcorrupt_skymap))
		datt = hdu_gridd[0].data*tp
		vl = datt.copy()
#		for k1 in [0,1,2]:
		
		# *****  5z ********
		i =  v1[0, 0, ...]
		mk = i.copy()
		q =  v1[0, 1, ...]
		u =  v1[0, 2, ...]
		st = np.sqrt(q**2 + u**2)*1e1/i
		# *****  10z ********
		i =  v[0, 0, ...]
		q =  v[0, 1, ...]
		u =  v[0, 2, ...]		
		st1 = np.sqrt(q**2 + u**2)*1e1/i
		# *****  5z ********
		i =  xy[0, 0, ...]
		q =  xy[0, 1, ...]
		u =  xy[0, 2, ...]		
		xyst = np.sqrt(q**2 + u**2)*1e1/i
		# *****  10z ********
		i =  vl[0, 0, ...]
		mv = i.copy()
		q = vl[0, 1, ...]
		u =  vl[0, 2, ...]		
		vst = np.sqrt(q**2 + u**2)*1e1/i
		
		
#		hp.zoomtool.mollzoom(st, coord = ['G', 'C'], norm= 'linear', 
#				     xsize = 3200, unit = 'mK/beam', 
#				     margins=(0.01, 0.01, 0.0, 0.02), 
#				    ) 
		
#		hp.visufunc.mollview(st, coord = ['G', 'C'], norm= 'linear', 
#				     xsize = 3200, cbar=True, unit = 'mK/beam', 
#				     margins=(0.01, 0.01, 0.0, 0.02), notext = True) 
		nm = 500 #int(420)	#720
#		uk = savgol_filter(hp.sphtfunc.anafast(mk, lmax = nm), 5, 1)  
#		vk = savgol_filter(hp.sphtfunc.anafast(mv, lmax = nm), 5, 1)  
		u1 = savgol_filter(hp.sphtfunc.anafast(st1, lmax = nm), 5, 1)  #TRUE 
		u = savgol_filter(hp.sphtfunc.anafast(st, lmax = nm), 5, 1)   # GP  
		uxy = savgol_filter(hp.sphtfunc.anafast(xyst, lmax = nm), 5, 1)  
		uv = savgol_filter(hp.sphtfunc.anafast(vst, lmax = nm), 5, 1)  
		lmax = 2.0*np.arange(len(u1))/4.8 #np.arange(len(u1))/3.0
#		x = np.linspace(0, len(u1), len(u1))
#		plt.semilogy(lmax,100*lmax*(lmax +1.)*u1/2.0*np.pi, label = '${QU \longrightarrow I_{5z}^{ANT1}} $ ')
#		plt.semilogy(lmax,100*lmax*(lmax +1.)*u/2.0*np.pi, linestyle='--',linewidth=3, 
#		             label = ' ${QU \longrightarrow I_{10z}^{ANT1}} = %4.4g $ ' %(abs(u1.std()-u.std())/u1.std()))
#	        plt.semilogy(lmax, abs(u1 - u),     
#		plt.semilogy(lmax,100*lmax*(lmax +1.)*uxy/2.0*np.pi,  label = ' ${QU \longrightarrow I_{5z}^{ANT5}}$' )
#		plt.semilogy(lmax,100*lmax*(lmax +1.)*uv/2.0*np.pi, linestyle='--',  linewidth=3, 
#		             label = ' ${QU \longrightarrow I_{10z}^{ANT5}} = %4.4g$' %(abs(uxy.std()-uv.std())/uv.std()))
#		
		
#		# ---- Linear Leakage terms ----------
#		plt.semilogy(lmax,100*lmax*(lmax +1.)*u1/2.0*np.pi, label = r'$QU \longrightarrow I_{T} $ ')
##		plt.semilogy(lmax,100*lmax*(lmax +1.)*uv/2.0*np.pi, label = ' $QU \longrightarrow I_{H}$' )
#		plt.semilogy(lmax,1000*lmax*(lmax +1.)*u/2.0*np.pi, linestyle='--',linewidth=3, label = r' $QU \longrightarrow I_{GP} $ ' )
#	        plt.semilogy(lmax,1000*lmax*(lmax +1.)*uxy/2.0*np.pi, linestyle='--',linewidth=3, label = r' $QU \longrightarrow I_{XY}$' )  
##		plt.semilogy(lmax,1000*lmax*(lmax +1.)*uxy/2.0*np.pi,  label = ' $QU \longrightarrow I_{XY}$' )
##		plt.semilogy(lmax,100*lmax*(lmax +1.)*uv/2.0*np.pi, linestyle='--',  linewidth=3, label = ' $QU \longrightarrow I_{H}$' )
##                plt.semilogy(lmax, abs(uxy - uv),  label = r'$Er_{L}^{ANT5}$')
#		plt.semilogy(lmax, abs(u1 - u),  label = r'$Er_{L}^{GP}$') 
#                plt.semilogy(lmax, abs(u1 - uxy),  label = r'$Er_{L}^{XY}$')		
#		
#		plt.xlim(xmin = 1)		
#		plt.xlabel(r" $l$ (multipole)",  fontsize="15")
#		plt.ylabel(r" Percentage Leakage", fontsize="15")
#		plt.legend(loc='best' , framealpha=0.5)
#		# ---- Linear Leakage terms ----------
		
		
#		for sk in [st, st1, xyst, vst]:
#			j = 1
#			hp.visufunc.mollview(sk, fig=2, coord = ['G', 'C'], norm= 'linear', 
#					     xsize = 3200, cbar=True, unit = 'mK/beam', sub = (2, 2, j),
#					     margins=(0.01, 0.01, 0.0, 0.02), notext = True) 
#		        j+=1
#		
#		plt.show()
		
		hi = np.array(pf.open('/home/narh/papercodes/SimHI_001.fits')[1].data)
		hi = [hi[ii][0]*1e-1 for ii in xrange(3145728)]
		# *****  5z ********
		i =  v1[0, 0, ...]
		mk = i.copy()
		q =  v1[0, 1, ...]
		u =  v1[0, 2, ...]
		st = np.sqrt(q**2 + u**2)*1e1
		# *****  10z ********
		i =  v[0, 0, ...]
		q =  v[0, 1, ...]
		u =  v[0, 2, ...]		
		st1 = np.sqrt(q**2 + u**2)*1e1
		# *****  5z ********
		i =  xy[0, 0, ...]
		q =  xy[0, 1, ...]
		u =  xy[0, 2, ...]		
		xyst = np.sqrt(q**2 + u**2)*1e1
		# *****  10z ********
		i =  vl[0, 0, ...]
		mv = i.copy()
		q = vl[0, 1, ...]
		u =  vl[0, 2, ...]		
		vst = np.sqrt(q**2 + u**2)*1e1
		
		
#		hp.zoomtool.mollzoom(st, coord = ['G', 'C'], norm= 'linear', 
#				     xsize = 3200, unit = 'mK/beam', 
#				     margins=(0.01, 0.01, 0.0, 0.02), 
#				    ) 
		
#		hp.visufunc.mollview(st, coord = ['G', 'C'], norm= 'linear', 
#				     xsize = 3200, cbar=True, unit = 'mK/beam', 
#				     margins=(0.01, 0.01, 0.0, 0.02), notext = True) 
#		nm = 400 #int(420)	#720
#		uk = savgol_filter(hp.sphtfunc.anafast(mk, lmax = nm), 5, 1)  
#		vk = savgol_filter(hp.sphtfunc.anafast(mv, lmax = nm), 5, 1)  
		from collections import OrderedDict
		linestyles = OrderedDict(
					    [('solid',               (0, ())),
					     ('loosely dotted',      (0, (1, 10))),
					     ('dotted',              (0, (1, 5))),
					     ('densely dotted',      (0, (1, 1))),

					     ('loosely dashed',      (0, (5, 10))),
					     ('dashed',              (0, (5, 5))),
					     ('densely dashed',      (0, (5, 1))),

					     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
					     ('dashdotted',          (0, (3, 5, 1, 5))),
					     ('densely dashdotted',  (0, (3, 1, 1, 1))),

					     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
					     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
					     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
		
		plt.figure(2)
		u1 = savgol_filter(hp.sphtfunc.anafast(st1, lmax = nm), 5, 1)  #TRUE 
		u = savgol_filter(hp.sphtfunc.anafast(st, lmax = nm), 5, 1)   # GP  
		uxy = savgol_filter(hp.sphtfunc.anafast(xyst, lmax = nm), 5, 1)  
		uv = savgol_filter(hp.sphtfunc.anafast(vst, lmax = nm), 5, 1)
		h1 = hp.sphtfunc.anafast(hi, lmax = nm) #savgol_filter(hp.sphtfunc.anafast(vst, lmax = nm), 5, 1) 
#		I = hp.sphtfunc.anafast(v[0, 0, ...], lmax = nm)
#		I = hp.sphtfunc.anafast(dat[iter,0,...], lmax = nm)
		
#		plt.semilogy(lmax,1*lmax*(lmax +1.)*I/2.0*np.pi, label = r'$ I_{T} $ ')
		fig = figur(fig_num = 2 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax = fig.add_subplot(111)
#		plt.semilogy(lmax,1*lmax*(lmax +1.)*u1/2.0*np.pi, label = r'$ |Q + iU|_{T} $ ')
#		plt.semilogy(lmax,10*lmax*(lmax +1.)*u/2.0*np.pi, linestyle=linestyles['densely dashed'], linewidth=3, label = r' $|Q + iU|_{GP} $ ' )
#	        plt.semilogy(lmax,10*lmax*(lmax +1.)*uxy/2.0*np.pi, linestyle='--',linewidth=3, label = r' $|Q + iU|_{XY}$' )
#	        plt.semilogy(lmax,10*lmax*(lmax + 1.)*uv/2.0*np.pi, linestyle = linestyles['dashdotted'],  linewidth=3, label = r'$ |Q + iU|_{JVLA} $ ')
#	        plt.semilogy(lmax,lmax*(lmax +1.)*h1/2.0*np.pi, marker='o', markersize=4, markerfacecolor='magenta', label = r' $HI \text{at} z \approx 0.67$' )
	        
	        ax.semilogy(lmax,1*lmax*(lmax +1.)*u1/2.0*np.pi, label = r'$ |Q + iU|_{T} $ ')
		ax.semilogy(lmax,100*lmax*(lmax +1.)*u/2.0*np.pi, linestyle=linestyles['densely dashed'], linewidth=3, label = r' $|Q + iU|_{GP} $ ' )
	        ax.semilogy(lmax,100*lmax*(lmax +1.)*uxy/2.0*np.pi, linestyle='--',linewidth=3, label = r' $|Q + iU|_{XY}$' )
	        ax.semilogy(lmax,100*lmax*(lmax + 1.)*uv/2.0*np.pi, linestyle = linestyles['dashdotted'],  linewidth=3, label = r'$ |Q + iU|_{JVLA} $ ')
	        ax.semilogy(lmax,10*lmax*(lmax +1.)*h1/2.0*np.pi, marker='o', markersize=4, markerfacecolor='magenta', label = r' $HI \, at \, z \approx 0.67$' )
	        # OR use use \mbox 
	        		
		ax.set_xlim(xmin = 1)		
		ax.set_xlabel(r" $l$ (multipole)",  fontsize="15")
		ax.set_ylabel(r" $\Delta T \, [K]$", fontsize="15")
		ax.legend(loc='best' , framealpha=0.5, ncol=2)
		ax.grid()
#	        plt.xlim(xmin = 1)		
#		plt.xlabel(r" $l$ (multipole)",  fontsize="15")
#		plt.ylabel(r" $\Delta T \, [K]$", fontsize="15")
#		plt.legend(loc='best' , framealpha=0.5, ncol=2)
#		plt.show()
		fig.savefig('linear_leak.png')
		fig.clear()
		raise
		'''
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
		
#		"""

		# convolved maps

		conv_healpix_plott(conv_map = v, Tmap = dat[iter], figname = '%s_Tchan_%d.png' %(strt_sav_file, iter), fig_num = 1 ,
		                   w_size = 18, h_size = 10, DPI = 500)
		conv_healpix_plott(conv_map = v1, Tmap = dat[iter], figname = '%s_Dchan_%d.png' %(strt_sav_file, iter), fig_num = 2 , 
		                   w_size = 18, h_size = 10, DPI = 500)
		f_healpix_plott(conv_map = dat[iter], figname = '%s_fchan_%d' %(strt_sav_file, iter), fig_num = 3 , 
		                   w_size = 18, h_size = 10, DPI = 500)
           	
#		                   
		                   
		                   
		                   
		print '\n >>> Plotting Absolute difference btn the conv maps'
		#
		# Relative Error
		conv_errors = [] 
		from scipy.stats import sem
		for jj in xrange(16):
			
			if jj in [3, 7, 11, 12, 13, 14, 15]:
				continue
			else:
			
				m1 = v1[0, jj, ...]  #**2 #abs(v1[0, jj, ...]) 
				m2 = v[0, jj, ...]   #**2 #abs(v[0, jj, ...])
				cc = (m1 - m2)
				cc = cc/abs(v[0, jj, ...]).mean()	#sem(cc/cc.mean(), axis=None, ddof=0)
				#print xg
				
				#cc = (m1 - m2)**2
				conv_errors.append(cc)
#				print jj
		
		np.save('gphase_error_%d' %iter, np.array(conv_errors)) 
#		m = [1e5, 1e5, 1e4, 1e5, 2000, 7000, 2e5, 4000, 8000]
#		m = np.load('mask1.npy')
#		error_conv_healpix_plott(conv_map = np.array(conv_errors), figname = '%s_Echan_%d.png' %(strt_sav_file, iter),
#		 			fig_num = 15, w_size = 18, h_size = 10, DPI = 500, mask=None)
#		print '\n >>> done !!!!!!!!!!'
##		"""
#		raise
#		
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
		'''
		nbox = 19
		#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Stokes I
		
		nm = 500 #int(400)	#720
		u1 = hp.sphtfunc.anafast(dat[iter,0,...], lmax = nm) 		# True Sky		
		m_IT = v[0, 0, ...] + savgol_filter(v[0, 1, ...] + v[0, 2, ...], 31, 5, deriv = 2) 
		m_ID = v1[0, 0, ...] + savgol_filter(v1[0, 1, ...] + v1[0, 2, ...], 31, 5, deriv = 2)
		u = hp.sphtfunc.anafast(m_IT,  lmax = nm) 				# Conv I^T
		ud = hp.sphtfunc.anafast(m_ID, lmax = nm)				# Conv. I^D
		ul = hp.sphtfunc.anafast(m_IT - m_ID, lmax = nm)
#		p = np.polyfit(x,rb[:,0], 20)
		
		#Leakage term  
		#ul = hp.sphtfunc.anafast(m_IT - m_ID, lmax = nm) 
		#
		lmax = 2.0*np.arange(len(u1))/4.8 #np.arange(len(u1))/3.5
		#np.save('lmax_%d' %iter, np.array(lmax)) 
		fig = figur(fig_num = 5 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
#		p = np.polyfit(lmax,u1, 2)
#		u1 = np.polyval(p, lmax)
#		p = np.polyfit(lmax,u, 2)
#		u = np.polyval(p, lmax)
#		p = np.polyfit(lmax,ud, 2)
#		ud = np.polyval(p, lmax)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $I$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $I$ with ${B_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$I$ with ${B_{GP}}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $I$')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $I$')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",  fontsize=fsz)
		ax1.set_ylabel(r"  $ \Delta T \, [K]$", fontsize=fsz)
				
		ax2 = ax1.twinx()		
		
		ax2.semilogy(lmax,u/u1, '-y', label = '${B_{T}}$ ')
#		print (u/u1).max()
		
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${B_{GP}}$ ') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize=fsz)		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='best', framealpha = 0.5)		

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
		
#		p = np.polyfit(lmax,u1, 2)
#		u1 = np.polyval(p, lmax)
#		p = np.polyfit(lmax,u, 2)
#		u = np.polyval(p, lmax)
#		p = np.polyfit(lmax,ud, 2)
#		ud = np.polyval(p, lmax)
#		ax1.semilogy(lmax,lmax*(lmax +1.)* u1/2.0*np.pi, label = 'True Stokes $Q$ ')
#		ax1.semilogy(lmax,lmax*(lmax +1.)* u/2.0*np.pi, label = 'Conv. $Q$ with ${B_{5z}}$ ')		 
#		ax1.semilogy(lmax,lmax*(lmax +1.)* ud/2.0*np.pi, label = 'Conv.$Q$ with ${B_{10z}}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $Q$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $Q$ with ${B_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$Q$ with ${B_{GP}}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in  $Q$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)", fontsize=fsz)
		ax1.set_ylabel(r"  $ \Delta T \, [K]$", fontsize=fsz)
				
		ax2 = ax1.twinx()		
		ax2.semilogy(lmax, u/u1, '-y', label = '${B_{T}}$')
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${B_{GP}}$') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize=fsz)		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='best', framealpha = 0.5)		

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
		
#		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $U$ ')
#		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $U$ with ${B_{5z}}$ ')		 
#		ax1.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv.$U$ with ${B_{10z}}$ ')
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True Stokes $U$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv. $U$ with ${B_{T}}$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)* ud/2.0*np.pi, label = 'Conv.$U$ with ${B_{GP}}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in  $U$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $U$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",fontsize=fsz)
		ax1.set_ylabel(r" $ \Delta T \, [K]$",fontsize=fsz)
				
		ax2 = ax1.twinx()
		ax2.semilogy(lmax, u/u1, '-y', label = '${B_{T}}$')
		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${B_{GP}}$') 		
#		ax2.semilogy(lmax,u/u1, '-y', label = '${B_{5z}}$ ')
#		ax2.semilogy(lmax,ud/u1, 'darkviolet', label = '${B_{10z}}$') 		
		ax2.set_ylabel(r'Normalized Beam Pattern', fontsize=fsz)		#	, color='b'
		
		#	Legend
		h1, l1 = ax1.get_legend_handles_labels()
		h2, l2 = ax2.get_legend_handles_labels()		
		ax1.legend(h1+h2, l1+l2, loc='best', framealpha = 0.5)		

		fig.savefig('uupwchan_%d.png' %iter)
		fig.clear()
		
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		q = hp.sphtfunc.anafast(dat[iter, 1,...], lmax = nm)
		c2 =	v[0, 5, ...]
		c3 =	savgol_filter(v[0, 1, ...], 21, 3, deriv = 2)      #+ savgol_filter(v[0, 6, ...], 51, 5, deriv = 4) 		
#		m_QD =	v1[0, 5, ...] + savgol_filter(v1[0, 6, ...], 51, 5, deriv = 4) 
		qi = hp.sphtfunc.anafast(c3, lmax = nm)
		qq = hp.sphtfunc.anafast(c2, lmax = nm)
		fig = figur(fig_num = 8 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*q/2.0*np.pi, label = ' $Q_{T}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*qi/2.0*np.pi, label = '$ Q \longrightarrow I$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*qq/2.0*np.pi, label = '$ Q \longrightarrow Q$ ')
#		ax1.semilogy(lmax,lmax*(lmax +1.)*smooth(y=ul, box_pts=nbox)/2.0*np.pi, label = 'Error in  $U$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $U$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",fontsize=fsz)
		ax1.set_ylabel(r" $ \Delta T \, [K]$",fontsize=fsz)
		ax1.legend(loc='upper right', framealpha = 0.5)
		fig.savefig('QIQchan_%d.png' %iter)
		fig.clear()
		#
		
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		u = hp.sphtfunc.anafast(dat[iter, 2,...], lmax = nm)
		c2 =	v[0, 10, ...]
		c3 =	v[0, 2, ...]    #+ savgol_filter(v[0, 6, ...], 51, 5, deriv = 4) 		
#		m_QD =	v1[0, 5, ...] + savgol_filter(v1[0, 6, ...], 51, 5, deriv = 4) 
		qi = hp.sphtfunc.anafast(c3, lmax = nm)
		qq = hp.sphtfunc.anafast(c2, lmax = nm)
		fig = figur(fig_num = 9 , w_size = 10, h_size = 6, DPI = 200) 		#plt.figure()
		ax1 = fig.add_subplot(111)
		
		ax1.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = ' $U_{T}$ ')
		ax1.semilogy(lmax,lmax*(lmax +1.)*qi/2.0*np.pi, label = '$ U \longrightarrow I$ ')		 
		ax1.semilogy(lmax,lmax*(lmax +1.)*qq/2.0*np.pi, label = '$ U \longrightarrow U$ ')
#		ax1.semilogy(lmax,lmax*(lmax +1.)*smooth(y=ul, box_pts=nbox)/2.0*np.pi, label = 'Error in  $U$ ')
		#ax1.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Error in $U$ ')
		#		
		ax1.set_xlabel(r" $l$ (multipole)",fontsize=fsz)
		ax1.set_ylabel(r" $ \Delta T \, [K]$",fontsize=fsz)
		ax1.legend(loc='upper right', framealpha = 0.5)
		fig.savefig('UIUchan_%d.png' %iter)
		fig.clear()
		"""
		if config.get('parameters', 'sav_percentage_leakage', 1) == 'yes' :		
			
			stokes_leakage.append(savgol_filter(ul/u1 * 100, 21, 3))
			d1 =  { 'l_max' : np.array(lmax), 'IQU' : np.array(stokes_leakage)} 
			np.savez('%s_fsky_percentage_leakage' %(strt_sav_file), **d1)"""  
			
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

