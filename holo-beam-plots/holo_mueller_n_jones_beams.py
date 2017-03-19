#
import glob
import pyfits 
import pywcsgrid2
import pywcs
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
from matplotlib import pylab as plt
import numpy as np
from scipy.stats import expon
import scipy.interpolate as interpolate
import warnings
warnings.filterwarnings("ignore")
import os, sys
import shutil
from scipy import signal
import time

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
	    gh = pywcsgrid2.GridHelper(wcs=header)
	    gh.locator_params(nbins=3)
	    #print "done !!!!!!!"
            #raise
	    g = axes_grid.ImageGrid(fig, 111,
		                    nrows_ncols=(nrows, ncols),
		                    ngrids=None,
		                    direction='row',
		                    axes_pad=0.02, add_all=True,
		                    share_all=True, aspect=True,
		                    label_mode='L', cbar_mode=None,
		                    cbar_location='right',
		                    axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh))) ##
	    # make colorbar
	    cax = []
	    for iter in range(i):		
		ax = inset_axes(g[iter],
		             width="45%", # width = 10% of parent_bbox width ---- "5%"
		             height="5%", # height : 50%  ------- "70%"
		             loc=9)
		cax.append(ax)
	    return g, cax
	    
#
def beam_model(beam_amplitude, beam_phase):
	    #	
	    #	Generating Jones Terms    
	    xx = beam_amplitude[0, :, :] +  1j*beam_phase[0, :, :]
	    xy = beam_amplitude[1, :, :] +  1j*beam_phase[1, :, :]
	    yx = beam_amplitude[2, :, :] +  1j*beam_phase[2, :, :]
	    yy = beam_amplitude[3, :, :] +  1j*beam_phase[3, :, :]
	    #
	    
	    '''cI = np.complex(0.,1.)
	    xx = np.transpose(np.array(0.5*(xx1 + xy1 + yx1 + yy1)))
	    xy = np.transpose(np.array(0.5*cI*(xx1 - xy1 + yx1 - yy1)))
	    yx = np.transpose(np.array(0.5*cI*(-xx1 - xy1 + yx1 + yy1)))
	    yy = np.transpose(np.array(0.5*(xx1 - xy1 - yx1 + yy1)))'''
	    #	Generating Mueller Terms
	    #
	    print "\n >> Mueller Matrix Conversion:"
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
	    
	    return M
	    
	    
#
def Mueller_image(fig, header, nrows, ncols, i , Mueller_matrix, figname):
	    
	    g2, cax2 = setup_axes(fig, header, nrows, ncols, i)
	    #
	    # Plot M Images
	    #
	    
	    cmap = plt.get_cmap('jet')	
	    images = []
	    for i, ax in enumerate(g2):
		channel = Mueller_matrix[i]			
		im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos")
		if i in [0, 5, 10, 15]:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel),  np.max(channel)], orientation = "horizontal") 
		    # ticks = [0,  np.max(channel)]
		   cax2[i].xaxis.set_ticks_position("bottom")  # ------ just added
		else:
		   plt.colorbar(im, cax=cax2[i], format = '%.0e', ticks = [np.min(channel), np.max(channel)], orientation = "horizontal")
		   cax2[i].xaxis.set_ticks_position("bottom")  # ------ just added
	        ax.axis["bottom"].major_ticklabels.set(fontsize=15)
		ax.axis["left"].major_ticklabels.set(fontsize=15)
		images.append(im)
	    plt.savefig(figname)
	    return None 
	    
	    
#  
def Jones_image(fig, header, nrows, ncols, i , Jones_matrix, figname):
    g1, cax1 = setup_axes(fig, header, nrows, ncols, i)		
    #
    # Plot Jones Images
    #
    print "\n >> Jones Matrix Conversion:"
    cmap = plt.get_cmap('jet')		
    images = []
    for i, ax in enumerate(g1):
        channel =  Jones_matrix[i]		
        im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos")
        plt.colorbar(im, cax=cax1[i], format = '%.0e', ticks = [np.min(channel),  np.max(channel)], orientation = "horizontal")
        cax1[i].xaxis.set_ticks_position("bottom")
        images.append(im)
    plt.savefig(figname)
    return None
    
    
# 
def figur(n1, n2, n3, n4):
	    fx = plt.figure(n1, figsize=(n2, n3), dpi= n4)
	    return fx
	    
#    
def headers(Fits_filename):
	    hd = pyfits.getheader(Fits_filename)
	    return hd
	    
#
def run_test():
	#chan_num = 0
	ant_num = 5
	start = time.time()
        startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Start at %s" % startime
        #
        #jns = ['XX', 'XY', 'YX', 'YY'] 
        
        jns = ['LL', 'LR', 'RL', 'RR']
        ty = ['real', 'imag']
        k = ['Eampl.fits', 'Ephase.fits']
        
        filename = 'holo_a%d' %ant_num
	if os.path.isdir("%s" %filename):		
       		os.popen('rm -rf %s' %filename)
	os.mkdir('%s' %filename)
       		#os.system('mv -f *png %s '  %(filename))
#	else:		        
#		#os.makedirs('%s' %filename)
#		os.mkdir('%s' %filename)
	#
	cwd = os.getcwd()                            # Current directory path
        for iter in [0, 100, 230, 340, 450]:		#[0, 9, 19, 29, 39, 49]:	[700, 750, 800, 850, 900, 950]	# channel frequencies
        	
		for ii in xrange(len(ty)):
			m = []			         	    
			for jj in xrange(len(jns)):			
				beam_file = 'ant%d%s%s.fits' %(ant_num, jns[jj], ty[ii]) 
				hdu = openFitsFile('%s' %(beam_file))
				beam_data = hdu[0].data[iter]
				m.append(beam_data)			#
			
			beam_header = hdu[0].header
			beam_header.update('CTYPE1','RA', after='CRVAL1')
			beam_header.update('CTYPE2','DEC', after='CRVAL2')			 
			pyfits.writeto('%s' %(k[ii]), np.array(m), header=beam_header, clobber=True)
			#os.chdir('%s' %cwd)
			
		#os.chdir('%s' %filename)		
		beam_amp = 'Eampl.fits' 
		hdu = openFitsFile('%s' %(beam_amp))
		beam_data_amp = hdu[0].data
		beam_ph = 'Ephase.fits'		# phase 
		hdu = openFitsFile('%s' %(beam_ph))
		beam_data_ph = hdu[0].data
		beam_header = hdu[0].header
		mueller_beam = beam_model(beam_data_amp, beam_data_ph)
		
		os.chdir('%s' %filename)          	     # change current directory
		fv1 = figur(1, 12, 12, 80)
		
	    	Mueller_image(fv1, beam_header, 4, 4, 16 , mueller_beam, 
	    		     'vla_ant%d_chan%s_mueller_im.png' %(ant_num, str(iter).zfill(4)))
	        
	    	pyfits.writeto('vla_ant%d_chan%s_mueller_im.fits' %(ant_num, str(iter).zfill(4)),
	    	               mueller_beam, header = beam_header, clobber=True)
		#        
		fv1 = figur(2, 12, 12, 70)
	    	Jones_image(fv1, beam_header, 2, 2, 4 ,  beam_data_amp, 
	    	            'vla_ant%d_chan%s_jones_re_im.png' %(ant_num, str(iter).zfill(4)))		    	
	    	fv2 = figur(3, 12, 12, 70)
	    	Jones_image(fv2, beam_header, 2, 2, 4 , beam_data_ph, 
	    		    'vla_ant%d_chan%s_jones_im_im.png' %(ant_num, str(iter).zfill(4)))
	        os.chdir('%s' %cwd)
	        #chan_num += 1
    	#
        print 'Done !!!!!!!!!!!!!!!!!!!!'   	
     
        stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0		
        print "Total run time: %7.2f hours" % elasped_time


if __name__ == '__main__':
	run_test()
	
	

