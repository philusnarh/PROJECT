#!/usr/bin/env python
"""
primary_beam_model.py
====================

This program generates random dipoles corresponding to Aperture Illumination to model the Primary Beam of 
a Circular Aperture Array.

python setup.py install --user
"""
#
#from __future__ import division
import glob
import pyfits 
import pywcsgrid2
import pywcs
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
from matplotlib import pylab as plt
#from matplotlib import pylab as plt2
#from matplotlib import pylab as plt3
import numpy as np
from scipy.stats import expon
import scipy.interpolate as interpolate
import warnings
warnings.filterwarnings("ignore")
import os, sys
import shutil
#from ConfigParser import SafeConfigParser
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
def inverse_transform_sampling(hist, bin_edges, num_dipoles):
    cum_values = np.cumsum(hist[1:]*np.diff(bin_edges)) 	# Cumulative Distribution Function for each bin
    #inv_cdf1 = interpolate.UnivariateSpline(cum_values, bin_edges[1:])	# Inverse Transform using scipy interpolate function
    inv_cdf1 = interpolate.InterpolatedUnivariateSpline(cum_values, bin_edges[1:], k=5)
    pdf = np.arange(0,1, 1./num_dipoles)
    return inv_cdf1(pdf)
#
def setup_axes(fig, header, nrows, ncols, i):
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
                            axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh))) ##
    # make colorbar
    cax = []
    for iter in range(i):
        #ax = g[iter]
        ax = inset_axes(g[iter],
                     width="45%", # width = 10% of parent_bbox width ---- "5%"
                     height="5%", # height : 50%  ------- "70%"
                     loc=9,
                     #bbox_to_anchor=(1.05, 0., 1, 1), # (1.01, 0, 1, 1), --------(0.025, 0, 1.05, 0.95)
                     #bbox_transform=g[iter].transAxes,
                     #borderpad=0.
                     )
        cax.append(ax)
    return g, cax
    #
def beam_model(beam_voltage, beam_phase, header):
    #
    #	Checking for any nan & inf
    beam_voltage = np.where((beam_voltage == 0.) | (np.isnan(beam_voltage)) |(np.isinf(beam_voltage)), 0.000012, beam_voltage) 
    beam_phase = np.where((beam_phase == 0.) | (np.isnan(beam_phase)) | (np.isinf(beam_phase)), 0.000012, beam_phase) 
    #
    #	Extracting the Amplitude
    #
    X_Amp_theta = beam_voltage[0, 0, ...] 
    X_Amp_phi = beam_voltage[0, 1, ...]  
    Y_Amp_theta = beam_voltage[0, 2, ...]  
    Y_Amp_phi = beam_voltage[0, 3, ...]
    #
    #	Extracting the Phase
    #
    X_Phase_theta = beam_phase[0, 0, ...]    
    X_Phase_phi = beam_phase[0, 1, ...]      
    Y_Phase_theta = beam_phase[0, 2, ...]    
    Y_Phase_phi = beam_phase[0, 3, ...]
    #
    #	Generating Jones Terms
    #
    Re_X_theta = X_Amp_theta*np.cos(X_Phase_theta)
    Im_X_theta = X_Amp_theta*np.sin(X_Phase_theta)
    Re_X_phi = X_Amp_phi*np.cos(X_Phase_phi) 
    Im_X_phi = X_Amp_phi*np.sin(X_Phase_phi) 
    Re_Y_theta = Y_Amp_theta*np.cos(Y_Phase_theta)
    Im_Y_theta = Y_Amp_theta*np.sin(Y_Phase_theta)
    Re_Y_phi = Y_Amp_phi*np.cos(Y_Phase_phi)
    Im_Y_phi = Y_Amp_phi*np.sin(Y_Phase_phi)
    Xc_theta = Re_X_theta + 1j*Im_X_theta 
    Xc_phi = Re_X_phi + 1j*Im_X_phi 
    Yc_theta = Re_Y_theta + 1j*Im_Y_theta 
    Yc_phi = Re_Y_phi + 1j*Im_Y_phi 
    xx = Xc_theta*np.conjugate(Xc_theta) +  Xc_phi*np.conjugate(Xc_phi) 
    xy = Xc_theta*np.conjugate(Yc_theta) +  Xc_phi*np.conjugate(Yc_phi) 
    yx = Yc_theta*np.conjugate(Xc_theta) +  Yc_phi*np.conjugate(Xc_phi) 
    yy = Yc_theta*np.conjugate(Yc_theta) +  Yc_phi*np.conjugate(Yc_phi)
    JJr = []
    JJi = []
    JJr.append(xx.real)
    JJr.append(xy.real)
    JJr.append(yx.real)		#	yx.real
    JJr.append(yy.real)
    JJi.append(xx.imag)
    JJi.append(xy.imag)
    JJi.append(xy.imag)		#	yx.imag
    JJi.append(yy.imag)
    JJr = np.array(JJr)
    JJi = np.array(JJi)
    #
    pyfits.writeto('Jones_2_x_2_Real_images.fits',JJr, header[0:38], clobber=True)
    pyfits.writeto('Jones_2_x_2_Imag_images.fits',JJi, header[0:38], clobber=True)
    #
    #	Generating Mueller Terms
    #
    print "\n >> Mueller Matrices Conversion:"
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
    return JJr, JJi, M
#
def Jones_image(fig, header, nrows, ncols, i , Jones_matrix, figname):
    g1, cax1 = setup_axes(fig, header, nrows, ncols, i)		
    #
    # Plot Jones Images
    #
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
        images.append(im)
    plt.savefig(figname)
    return None 
#    
def figur(n1, n2, n3, n4):
    fx = plt.figure(n1, figsize=(n2, n3), dpi= n4)
    return fx
    
def headers(Fits_filename):
    hd = pyfits.getheader(Fits_filename)
    return hd
#    
def run_test():
    #
    option = OptionParser()
    option.set_description('Adds Kernel size for convolution & FITS filename')
    option.set_usage('%prog [options] arg1 arg2')
    option.add_option('-n','--nd',dest='ndipoles',default=1000,type = int,help='Enter the number of dipoles for the model')
    option.add_option('-d','--ap',dest='dish_size',default=1,type = int,help='Enter the value of Aperture size')
    option.add_option('-s','--subf',dest='subref_size',default=1,type = int,help='Enter the value of Subreflector size')
    option.add_option('-a','--na',dest='num_ant',default=1,type = int,help='Enter the number of antennas')
    option.add_option('-t','--typ',dest='rad_type',default=1,type = int,help='Enter the value 0 for Gaussian OR 1 for Exponential,\
    optional')
    #option.add_option('-f','--path', dest='file_path', default=None, action="store", type="string", help='Enter setup path')
    options,args = option.parse_args(sys.argv[1:])
    #path = args[0]
    #print path
    #
    #
    bin_edges = np.arange(0.0,1,0.01) 
    #
    #	Determines the Type of Radial Distribution for the Aperture Illumination
    # 
    Type = options.rad_type
    if Type == 0:
       Jpdf = (1/(2*np.pi))*np.exp(-(bin_edges**2)/2)  # Gaussian distribution of random deviates
       plt.figure(1)
       rad = inverse_transform_sampling(Jpdf, bin_edges, options.ndipoles)
       plt.hist(rad,80,normed=True)
       plt.xlabel('Radius')
       plt.ylabel('Aperture Illumination Function')
       plt.title('Radial Distribution of Aperture Illumination')
       plt.savefig('radial_distribution.png')		# Saves radial distribution for Gaussian
    elif Type == 1:
       #Jrad = np.exp(-bin_edges) 			# Exponential distribution of random deviates
       Jrad = expon.pdf(bin_edges, loc=0,scale=1)
       plt.figure(1)
       rad = inverse_transform_sampling(Jrad, bin_edges, options.ndipoles)
       plt.hist(rad, 80, normed=True)
       plt.xlabel('Raduis')
       plt.xlim(0, 4)
       plt.ylabel('Aperture Illumination Function')
       plt.title('Radial Distribution of Aperture Illumination')
       plt.savefig('radial_distribution.png')		# Saves radial distribution for Exponential
    else:
       print "ENTER THE VALUE 0 OR 1 TO CHOOSE RADIAL TYPE"
       raise
    #
    #	Circular Aperture Array Illumination layout
    #
    r1 = options.subref_size/6.0
    a = rad[np.where((rad<=1) & (rad > r1))]
    U = np.random.uniform(0, 2*np.pi, len(a))
    plt.figure(2)
    x = options.dish_size/2*a*np.cos(U)			# X Orientation of dipoles
    y = options.dish_size/2*a*np.sin(U)			# Y Orientation of dipoles
    s1 = np.where(np.logical_and(x >= - r1, x <= r1)) 
    s1 = np.take(s1, range(len(s1[0])- len(s1[0])/8))
    x1 = np.delete(x, s1)
    y1 = np.delete(y, s1)
    s2 = np.where(np.logical_and(y1 >= - r1, y1 <= r1))
    s2 = np.take(s2, range(len(s2[0])- len(s2[0])/8))
    x1 = np.delete(x1, s2)
    y1 = np.delete(y1,s2)
    plt.plot(x1, y1, 'b+') 
    plt.xlabel('X/ m', fontsize="15")
    plt.ylabel('Y/ m', fontsize="15")
    plt.title('Station Setup')
    plt.savefig('statin_setup.png') 			#Saves Station Setup
    plt.show()
    #
    #	Save Dipoles
    #
    for num in range(0,options.num_ant,1):        
        os.makedirs('folder/station00%d'%num)
    os.chdir('folder/')
    np.savetxt('layout.txt' ,np.array([x1,y1, np.zeros(len(x1))]).T, delimiter = ',') 
    #shutil.copyfile('layout.txt', 'station???')
    os.system('echo station00* | xargs -n 1 cp layout.txt -f')
    #os.remove('layout.txt')
    os.chdir('../')
    if os.path.isdir("station000"):		# OR use os.path.exists("/home/el/myfile.txt")
       os.popen('rm -rf station*')
    #for path in glob.glob("station*/layout.txt"):
        #shutil.rmtree('path')
        #os.remove('path')
    os.system('mv -f folder/station* ./')
    os.system('rm -rf folder')
    #
    # Oskar Generates Non-distorted Beam Patterns  
    #
    #
    print " ------------------------------------------------------"
    print "\n >> Oskar Beampattern Simulation Begins:"
    print " ------------------------------------------------------"
    if os.path.exists("true_model_oskar_congifuration_setup1.ini"):
       os.system('oskar_sim_beam_pattern true_model_oskar_congifuration_setup1.ini')
    else:
       print "Error: config filename must be: 'true_model_oskar_congifuration_setup1.ini' "
       #print "Error: cannot open %s"%filename
       raise
    #
    config = ConfigParser.ConfigParser()
    config.read('true_model_oskar_congifuration_setup1.ini')
    filename = config.get('beam_pattern', 'root_path', 1)
    beam_filename = filename
    phase_filename = filename
    beam_filename+='_VOLTAGE.fits'
    phase_filename+='_PHASE.fits'
    #
    print "\n >> Opening Voltage & Phase Beampattern:"
    hduAmp = pyfits.open(beam_filename) 
    hduPh = pyfits.open(phase_filename) 
    datAmp = hduAmp[0].data[0]
    datPh = hduPh[0].data[0]
    #
    #	header
    #
    h = hduAmp[0].header 
    #
    jones_real_data, jones_imag_data, meuller_data = beam_model(datAmp, datPh, h)
    #fig = plt.figure(1, figsize=(12, 12), dpi=70)
    fv1 = figur(1, 12, 12, 70)
    Jones_image(fv1, h, 2, 2, 4 , jones_real_data, 'True_Jones_Real_Image.png')
    #fig1 = plt.figure(2, figsize=(12, 12), dpi=70)
    fv2 = figur(2, 12, 12, 70)
    Jones_image(fv2, h, 2, 2, 4 , jones_imag_data, 'True_Jones_Imag_Image.png')
    #fig2 = plt3.figure(3, figsize=(12, 12), dpi=80)
    fv3 = figur(3, 12, 12, 80)
    Mueller_image(fv3, h, 4, 4, 16 , meuller_data, 'True_Meuller_4_X_4_Images.png')
    pyfits.writeto('mueller_4_x_4_beam.fits',meuller_data,h[0:38], clobber=True)
    #
    #
    # Oskar Generates Distorted Beam Patterns  
    #
    #
    if os.path.exists("distort_model_oskar_congifuration_setup1.ini"):
       os.system('oskar_sim_beam_pattern distort_model_oskar_congifuration_setup1.ini')
    else:
       print "Error: config filename must be: 'distort_model_oskar_congifuration_setup1.ini' "
       raise
    #
    config = ConfigParser.ConfigParser()
    config.read('distort_model_oskar_congifuration_setup1.ini')
    filename = config.get('beam_pattern', 'root_path', 1)
    beam_filename = filename
    phase_filename = filename
    beam_filename+='_VOLTAGE.fits'
    phase_filename+='_PHASE.fits'
    #
    # Opening Distorted Voltage & Phase Beampattern
    #
    hduAmp = pyfits.open(beam_filename) 
    hduPh = pyfits.open(phase_filename) 
    datAmp = hduAmp[0].data[0]
    datPh = hduPh[0].data[0]
    #
    #	header
    #
    #h = hduAmp[0].header
    h = headers(beam_filename) 
    #
    dis_jones_real_data, dis_jones_imag_data, dis_meuller_data = beam_model(datAmp, datPh, h)
    fu1 = figur(4, 12, 12, 70)  #plt.figure(1, figsize=(12, 12), dpi=70)
    #
    Jones_image(fu1, h, 2, 2, 4 , dis_jones_real_data, figname = 'Distorted_Jones_Real_Image.png')
    fu2 = figur(5, 12, 12, 70)
    Jones_image(fu2, h, 2, 2, 4 , dis_jones_imag_data, figname = 'Distorted_Jones_Imag_Image.png')
    fu3 = figur(6, 12, 12, 80)
    #fig = plt.figure(1, figsize=(12, 12), dpi=70)
    Mueller_image(fu3, h, 4, 4, 16 , dis_meuller_data, figname = 'Distorted_Meuller_4_X_4_Images.png')
    pyfits.writeto('dis_mueller_4_x_4_beam.fits', dis_meuller_data,h[0:38], clobber=True)
    #
    # Compare Distorted & Non-distorted Beams
    #
    compare_beam = []
    for num in xrange(16):
        dif = dis_meuller_data[num] - meuller_data[num]
        compare_beam.append(dif)
    compare_beam = np.array(compare_beam)
    pyfits.writeto('mueller_difference_beam.fits',compare_beam,h[0:38], clobber=True)
    fw1 = figur(7, 12, 12, 80)
    #fig = plt.figure(1, figsize=(12, 12), dpi=70)
    Mueller_image(fw1, h, 4, 4, 16 , compare_beam, figname = 'Difference_Meuller_4_X_4_Images.png')
    ###############################################################################################
    '''hdu_new = openFitsFile('mueller_4_x_4_beam.fits')
    mueller_beam = hdu_new[0].data
    x_new, y_new = mueller_beam[0].shape
    measured_map= np.zeros(shape = (x_new, y_new), dtype=np.float64)
    Healpix_grid_I = np.load('mapII.npy')
    Healpix_grid_Q = np.load('mapQQ.npy')
    Healpix_grid_U = np.load('mapUU.npy')
    Healpix_grid_V = np.load('mapVV.npy')
    import copy
    temp22 = []
    fgrid_x_beam = []
    for iter1 in xrange(16):
        if iter1 in range(16)[0::4]:
           for k1 in xrange(x_new):
               for k2 in xrange(y_new):
                   measured_map[k1, k2] = np.sum(Healpix_grid_I[k1:x_new + k1, k2:y_new + k2]*mueller_beam[iter1])
	   #fgrid_x_beam.append(measured_map)
	   #print 'I',  fgrid_x_beam
	elif iter1 in range(16)[1::4]:
           for k1 in xrange(x_new):
               for k2 in xrange(y_new):
                   measured_map[k1, k2] = np.sum(Healpix_grid_Q[k1:x_new + k1, k2:y_new + k2]*mueller_beam[iter1])
	   #fgrid_x_beam.append(measured_map) 
	elif iter1 in range(16)[2::4]:
           for k1 in xrange(x_new):
               for k2 in xrange(y_new):
                   measured_map[k1, k2] = np.sum(Healpix_grid_U[k1:x_new + k1, k2:y_new + k2]*mueller_beam[iter1])
	   #fgrid_x_beam.append(measured_map) 
	else:
           for k1 in xrange(x_new):
               for k2 in xrange(y_new):
	           measured_map[k1, k2] = np.sum(Healpix_grid_V[k1:x_new + k1, k2:y_new + k2]*mueller_beam[iter1])
	fgrid_x_beam.append(measured_map)
	temp22.append(copy.deepcopy(measured_map))
	#print 'I',  np.array(fgrid_x_beam).shape
	AA = fgrid_x_beam
	print AA
	print np.array(AA).shape
    #fgrid_x_beam = np.array(fgrid_x_beam)
    #print 'I',  fgrid_x_beam
    print AA	
    print temp22
    pyfits.writeto('for_grid_x_trubeam.fits', temp22, h[0:38], clobber=True)
    fu4 = figur(8, 12, 12, 80)
    #fig = plt.figure(1, figsize=(12, 12), dpi=70)
    Mueller_image(fu4, h, 4, 4, 16 , temp22, figname = 'fgrid_x_trubeam_4_X_4_Images.png')'''
       
    ###############################################################################################
    print 'Done !!!!!!!!!!!'
    '''#
    # Oskar Interferometer Generation
    #
    print " ------------------------------------------------------"
    print "\n >> Checking for healpix map to run Oskar Interferometer:"
    print " ------------------------------------------------------"
    config = ConfigParser.ConfigParser()
    config.read('true_model_oskar_congifuration_setup1.ini')
    if config.get('sky', 'healpix_fits\map_units', 1)!= 'None':
       print "\n >> Oskar Interferometer Begins:"
       #os.system('oskar_sim_interferometer true_model_oskar_congifuration_setup1.ini')
    else:
       print "\n >> Oskar Configuration has no healpix map:"
    #shutil.move('folder/station???', './')'''
        #
    if os.path.exists("ap_array_model1.py~"):
       os.remove('ap_array_model1.py~')
    if os.path.exists("foreground_gridding.py~"):           # OR use os.path.exists("/home/el/myfile.txt")
           os.popen('rm -rf foreground_gridding.py~')

if __name__ == '__main__':
    run_test()
# ++++++++++++++++++++++++++++++++++++++++++++
#	run with the command below
# python ap_array_modelclass1.py -n 80000 -d 12 -a 7 -t 1 

