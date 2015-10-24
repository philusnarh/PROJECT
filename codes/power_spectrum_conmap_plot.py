import pyfits as pf
import numpy as np

from matplotlib import pylab as plt
import healpy as hp
from itertools import permutations, chain
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
	k = ['I', 'Q', 'U' ]
	#k1 =  sorted(k)
	for i in xrange(conv_map.shape[0] - 1):	
		
		if i in xrange(1):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, cbar=True, unit = 'Jy/pixel', margins=(0.01, 0.01, 0.0, 0.02), \
			notext = True, sub = (1,3 , i + 1))		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, cbar=True, unit = 'Jy/pixel',  margins=(0.01, 0.01, 0.0, 0.02),\
			notext = True, sub = (1, 3, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize=20, fontweight="bold") 
	fig.savefig(figname)
	plt.close(fig_num)
    	return None
#
def conv_healpix_plott_all(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['$Measured \, Stokes \, I^T$', '$Measured \, Stokes \, Q^T$', '$Measured \, Stokes \, U^T$', '$Measured \, Stokes \, I^D$', '$Measured \, Stokes \, Q^D$', \
	'$Measured \, Stokes \, U^D$', '$ Leakage \, in \, Stokes \, I$', '$ Leakage \, in \, Stokes \, Q$', '$ Leakage \, in \, Stokes \, U$' ]
	#k1 =  sorted(k)
	for i in xrange(conv_map.shape[0]):	
		#print conv_map.shape[0]
		if i in xrange(6):
			hp.visufunc.mollview(abs(conv_map[i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), \
			notext = True, sub = (3,3 , i + 1))		
		else:
			hp.visufunc.mollview(conv_map[i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, cbar=True, unit = 'Jy/beam',  margins=(0.01, 0.01, 0.0, 0.02),\
			notext = True, sub = (3, 3, i + 1)) 
		plt.title(r'%s'  % k[i], fontsize=20, fontweight="bold") 
	fig.savefig(figname)	
	plt.close(fig_num)
    	return None
#
def conv_healpix_plott(conv_map, figname, fig_num, w_size, h_size, DPI):
	fig = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
	k = ['I', 'Q', 'U', 'V']*4
	k1 =  sorted(k)
	j = 0
	for i in xrange(conv_map[0].shape[0]):	# [0, 1, 2, 4, 5, 6, 8, 9, 10]:	#
		if i in [0, 4, 8]:
			j = j + 1
			hp.visufunc.mollview(abs(conv_map[0, i,...]), coord = ['G', 'C'], norm= 'log', xsize = 3200, cbar=True, unit = 'Jy/beam', margins=(0.01, 0.01, 0.0, 0.02), \
			notext = True, sub = (3,3 , j))
			plt.title(r'$ %s \longleftarrow %s $' %(k1[i], k[i]), fontsize=20, fontweight="bold")
		elif i in [3, 7, 11, 12, 13, 14, 15]:
			continue
		
		else:
			j = j + 1
			hp.visufunc.mollview(conv_map[0, i,...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, cbar=True, unit = 'Jy/beam',  margins=(0.01, 0.01, 0.0, 0.02),\
			notext = True, sub = (3, 3, j )) 
			plt.title(r'$ %s \longleftarrow %s $' %(k1[i], k[i]), fontsize=20, fontweight="bold")
	fig.savefig(figname)
	plt.close(fig_num)
    	return None
#
def run_test():
	#v = np.load('IQUV_Trumap0.npy')
	#v1 = np.load('IQUV_Dismap0.npy')
#
#k =  [p for p in permutations('IQUV')]

#fig = figur(fig_num = 1, w_size = 18, h_size 10, DPI = 500 )
	#hdu_gridd = openFitsFile('560_4_786432_K2JanskypPixel.fits')
	hdu_gridd = openFitsFile('hh10IQUVR.fits')
	dat = hdu_gridd[0].data
	#conv_healpix_plott(conv_map = v, figname = 'tt.png', fig_num = 1, w_size = 18, h_size = 10, DPI = 500)
	for iter in xrange(1):
		print 'Started plotting chan:++++	', iter
		v = np.load('Trumap_%d.npy' %iter)
		v1 = np.load('Disumap_%d.npy'%iter)
		conv_healpix_plott(conv_map = v, figname = 'Tchan_%d.png' %iter, fig_num = 1 , w_size = 18, h_size = 10, DPI = 500)
		conv_healpix_plott(conv_map = v1, figname = 'Dchan_%d.png' %iter, fig_num = 2 , w_size = 18, h_size = 10, DPI = 500)
		f_healpix_plott(conv_map = dat[iter], figname = 'fchan_%d' %iter, fig_num = 3 , w_size = 14, h_size = 3, DPI = 500)
		#break
		lst = []
		ti = v[0, 0, ...] + v[0, 1, ...] + v[0, 2, ...] 
		tq = v[0, 4, ...] + v[0, 5, ...] + v[0, 6, ...] 
		tu = v[0, 8, ...] + v[0, 9, ...] + v[0, 10, ...] 
		di = v1[0, 0, ...] + v1[0, 1, ...] + v1[0, 2, ...] 
		dq = v1[0, 4, ...] + v1[0, 5, ...] + v1[0, 6, ...] 
		du = v1[0, 8, ...] + v1[0, 9, ...] + v1[0, 10, ...] 
		lst.append(ti)
		lst.append(tq)
		lst.append(tu)
		lst.append(di)
		lst.append(dq)
		lst.append(du)
		lst.append(ti - di)
		lst.append(tq - dq)
		lst.append(tu - du)
		lst = np.array(lst)
		#print lst.shape
		#raise
		#conv_healpix_plott_all(conv_map = v, figname = 'Cchan_T%d.png' %iter, fig_num = 1, w_size = 18, h_size = 10, DPI = 500)
		conv_healpix_plott_all(conv_map = lst, figname = 'All_chan_%d.png' %iter, fig_num = 4 , w_size = 18, h_size = 10, DPI = 500)
		#	
		# Stokes I
		fg = figur(fig_num = 5 , w_size = 10, h_size = 6, DPI = 200)
		nm = 800
		u1 = hp.sphtfunc.anafast(dat[iter,0,...], lmax =nm) 
		lmax = np.arange(len(u1)) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True $I$ ') 
		u = hp.sphtfunc.anafast(ti,lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv $I^T$ ')
		ud = hp.sphtfunc.anafast(di, lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv $I^D$ ')  
		ul = hp.sphtfunc.anafast(ti - di, lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Leakage $I^T - I^D$ ')
		plt.semilogy(lmax,u/u1, label = 'True Beam $I^{T}/ I$ ')
		plt.semilogy(lmax,ud/u1, label = 'Distorted Beam $I^{D}/ I$ ')
		lg = plt.legend(fancybox=True,loc = 'best') 
		lg.get_frame().set_alpha(0.5) 
		#plt.xlabel(r" Spatial Frequency $(MHz)$",fontsize="15") 
		plt.xlabel(r" $l$",fontsize="25")
		plt.ylabel(r" Power Spectrum $(dB)$",fontsize="15")
		fg.savefig('iipwchan_%d.png' %iter)
		plt.close(5)
		# Stokes Q
		fg = figur(fig_num = 6 , w_size = 10, h_size = 6, DPI = 200)
		u1 = hp.sphtfunc.anafast(dat[iter, 1,...], lmax = nm) 
		lmax = np.arange(len(u1)) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True $Q$ ') 
		u = hp.sphtfunc.anafast(v[0, 5, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv $Q^T$ ')
		ud = hp.sphtfunc.anafast(v1[0, 5, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv $Q^D$ ')   
		ul = hp.sphtfunc.anafast(v1[0, 5, ...] - v[0, 5, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Leakage $Q^T - Q^D$ ') 
		plt.semilogy(lmax,u/u1, label = 'True Beam $Q^{T}/ Q$ ')
		plt.semilogy(lmax,ud/u1, label = 'Distorted Beam $Q^{D}/ Q$ ')
		#plt.legend(loc = 'best')
		lg = plt.legend(fancybox=True,loc = 'best') 
		#lg = plt.legend(fancybox=True, loc='center left', bbox_to_anchor=(0.94, 0.5))
		lg.get_frame().set_alpha(0.5) 
		#plt.xlabel(r" Spatial Frequency $(MHz)$",fontsize="15") 
		plt.xlabel(r" $l$",fontsize="25")
		plt.ylabel(r" Power Spectrum $(dB)$",fontsize="15")
		fg.savefig('qqpwchan_%d.png' %iter)
		plt.close(6)
		#
		# Stokes Q1
		fg = figur(fig_num = 7 , w_size = 10, h_size = 6, DPI = 200)
		u1 = hp.sphtfunc.anafast(dat[iter, 1,...], lmax = nm) 
		lmax = np.arange(len(u1)) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True $Q$ ') 
		u = hp.sphtfunc.anafast(v[0, 5, ...] + v[0, 6, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv $Q^T$ ')
		ud = hp.sphtfunc.anafast( v1[0, 5, ...] + v1[0, 6, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv $Q^D$ ')   
		ul = hp.sphtfunc.anafast(v1[0, 5, ...] + v1[0, 6, ...] - v[0, 5, ...] - v[0, 6, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Leakage $Q^T - Q^D$ ') 
		plt.semilogy(lmax,u/u1, label = 'True Beam $Q^{T}/ Q$ ')
		plt.semilogy(lmax,ud/u1, label = 'Distorted Beam $Q^{D}/ Q$ ')
		lg = plt.legend(fancybox=True,loc = 'best') 
		#lg = plt.legend(fancybox=True, loc='center left', bbox_to_anchor=(0.94, 0.5))
		lg.get_frame().set_alpha(0.5)
		#plt.xlabel(r" Spatial Frequency $(MHz)$",fontsize="15") 
		plt.xlabel(r" $l$",fontsize="25")
		plt.ylabel(r" Power Spectrum $(dB)$",fontsize="15")
		fg.savefig('quchan_%d.png' %iter)
		plt.close(7)
		#
		# Stokes U
		fg = figur(fig_num = 8, w_size = 10, h_size = 6, DPI = 200)
		u1 = hp.sphtfunc.anafast(dat[iter, 2,...], lmax = nm) 
		lmax = np.arange(len(u1)) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True $U$ ') 
		u = hp.sphtfunc.anafast(v[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv $U^T$ ')
		ud = hp.sphtfunc.anafast(v1[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv $U^D$ ')   
		ul = hp.sphtfunc.anafast(v1[0, 10, ...] - v[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Leakage $U^T - U^D$ ') 
		plt.semilogy(lmax,u/u1, label = 'True Beam $U^{T}/ U$ ')
		plt.semilogy(lmax,ud/u1, label = 'Distorted Beam $U^{D}/ U$ ')
		lg = plt.legend(fancybox=True,loc = 'best') 
		#lg = plt.legend(fancybox=True, loc='center left', bbox_to_anchor=(0.94, 0.5))
		lg.get_frame().set_alpha(0.5)
		#plt.xlabel(r" Spatial Frequency $(MHz)$",fontsize="15") 
		plt.xlabel(r" $l$",fontsize="25")
		plt.ylabel(r" Power Spectrum $(dB)$",fontsize="15")
		fg.savefig('uupwchan_%d.png' %iter)
		plt.close(8)
		#
		# Stokes U1
		fg = figur(fig_num = 9, w_size = 10, h_size = 6, DPI = 200)
		u1 = hp.sphtfunc.anafast(dat[iter, 2,...],lmax = nm) 
		lmax = np.arange(len(u1)) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u1/2.0*np.pi, label = 'True $U$ ') 
		u = hp.sphtfunc.anafast(v[0, 9, ...] + v[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*u/2.0*np.pi, label = 'Conv $U^T$ ')
		ud = hp.sphtfunc.anafast(v1[0, 9, ...] + v1[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ud/2.0*np.pi, label = 'Conv $U^D$ ')   
		ul = hp.sphtfunc.anafast( v1[0, 9, ...] +v1[0, 10, ...] - v[0, 9, ...] - v[0, 10, ...], lmax = nm) 
		plt.semilogy(lmax,lmax*(lmax +1.)*ul/2.0*np.pi, label = 'Leakage $U^T - U^D$ ') 
		plt.semilogy(lmax,u/u1, label = 'True Beam $U^{T}/ U$ ')
		plt.semilogy(lmax,ud/u1, label = 'Distorted Beam $U^{D}/ U$ ')
		lg = plt.legend(fancybox=True,loc = 'best') 
		#lg = plt.legend(fancybox=True, loc='center left', bbox_to_anchor=(0.94, 0.5))
		lg.get_frame().set_alpha(0.5)
		plt.xlabel(r" $l$",fontsize="25") 
		#plt.xlabel(r" Spatial Frequency $(MHz)$",fontsize="15") 
		plt.ylabel(r" Power Spectrum $(dB)$",fontsize="15")
		fg.savefig('uqpwchan_%d.png' %iter)
		plt.close(9)
	#plt.close('all')
	print 'done wai !!!!!!!!!!!!'
if __name__ == '__main__':
	run_test()

