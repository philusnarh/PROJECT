#!/usr/bin/env python
"""
primary_beam_model.py
====================

This program generates random dipoles corresponding to Aperture Illumination to model the Primary Beam of 
a Circular Aperture Array.

"""
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
import scipy.interpolate as interpolate
import os, sys
import shutil
#from ConfigParser import SafeConfigParser
import ConfigParser
#
def inverse_transform_sampling(hist, bin_edges, num_dipoles):
    cum_values = np.cumsum(hist[1:]*np.diff(bin_edges)) 	# Cumulative Distribution Function for each bin
    inv_cdf1 = interpolate.UnivariateSpline(cum_values, bin_edges[1:])	# Inverse Transform using scipy interpolate function
    pdf = np.arange(0,1, 1./num_dipoles)
    return inv_cdf1(pdf)
#
def setup_axes(fig, header, nrows, ncols):
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
    ax = g[-1]
    cax = inset_axes(ax,
                     width="8%", # width = 10% of parent_bbox width
                     height="100%", # height : 50%
                     loc=3,
                     bbox_to_anchor=(1.01, 0, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0.
                     )
    return g, cax
    #
def run_test():
    #
    option = OptionParser()
    option.set_description('Adds Kernel size for convolution & FITS filename')
    option.set_usage('%prog [options] arg1 arg2')
    option.add_option('-n','--nd',dest='ndipoles',default=1000,type = int,help='Enter the number of dipoles for the model')
    option.add_option('-d','--ap',dest='dish_size',default=1,type = int,help='Enter the value of Aperture size')
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
       Jrad = np.exp(-bin_edges) 			# Exponential distribution of random deviates
       plt.figure(1)
       rad = inverse_transform_sampling(Jrad, bin_edges, options.ndipoles)
       plt.hist(rad, 80, normed=True)
       plt.xlabel('Raduis')
       plt.xlim(0, 3)
       plt.ylabel('Aperture Illumination Function')
       plt.title('Radial Distribution of Aperture Illumination')
       plt.savefig('radial_distribution.png')		# Saves radial distribution for Exponential
    else:
       print "ENTER THE VALUE 0 OR 1 TO CHOOSE RADIAL TYPE"
       raise
    #
    #	Circular Aperture Array Illumination layout
    #
    a = rad[np.where(rad<=1)]
    U = np.random.uniform(0, 2*np.pi, len(a))
    theta = np.random.uniform(0, np.pi, len(a))
    plt.figure(2)
    x = options.dish_size/2*a*np.cos(U)			# X Orientation of dipoles
    y = options.dish_size/2*a*np.sin(U)			# Y Orientation of dipoles		
    plt.plot(x, y, 'b+') 
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
    np.savetxt('layout.txt' ,np.array([x,y, np.zeros(len(x))]).T, delimiter = ',') 
    #print "Dipoles are saved as ....: dipoless.txt"
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
    #os.remove('aperture_array_model.py~')
    #
    # Oskar Beam Pattern Generation
    #
    print " ------------------------------------------------------"
    print "\n >> Oskar Beampattern Simulation Begins:"
    print " ------------------------------------------------------"
    os.system('oskar_sim_beam_pattern true_model_oskar_congifuration_setup1.ini')
    #
    print "\n >> Opening Voltage & Phase Beampattern:"
    hduA = pyfits.open('beam_pattern_VOLTAGE.fits') 
    hduPh = pyfits.open('beam_pattern_PHASE.fits') 
    datA = hduA[0].data
    #print datA.shape 
    datPh = hduPh[0].data
    #
    #	converting nan & inf
    #
    datA = np.where((datA == 0.) | (np.isnan(datA)) |(np.isinf(datA)), 0.000012, datA) 
    datPh = np.where((datPh == 0.) | (np.isnan(datPh)) | (np.isinf(datA)), 0.000012, datPh) 
    #
    #	Extracting the Amplitude
    #
    X_Amp_theta = datA[0, :, 0, ...] 
    X_Amp_phi = datA[0, :, 1, ...] 
    Y_Amp_theta = datA[0, :, 2, ...] 
    Y_Amp_phi = datA[0, :, 3, ...] 
    #
    #	Extracting the Phase
    #
    X_Phase_theta = datPh[0, :, 0, ...] 
    X_Phase_phi = datPh[0, :, 1, ...] 
    Y_Phase_theta = datPh[0, :, 2, ...] 
    Y_Phase_phi = datPh[0, :, 3, ...]
    #
    #	header
    #
    h = hduA[0].header 
    #
    #	Generating Complex Terms
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
    #
    #	Jones Matrix
    '''xx = Xc_theta
    xy = Xc_phi
    yx = Yc_theta
    yy = Yc_phi
    JJ = []
    JJ.append(xx.real)
    JJ.append(xy.real)
    JJ.append(yx.real)
    JJ.append(yy.real)
    np.save('JJ', JJ)'''
    ##############
    xx = Xc_theta*np.conjugate(Xc_theta) +  Xc_phi*np.conjugate(Xc_phi) 
    xy = Xc_theta*np.conjugate(Yc_theta) +  Xc_phi*np.conjugate(Yc_phi) 
    yx = Yc_theta*np.conjugate(Xc_theta) +  Yc_phi*np.conjugate(Xc_phi) 
    yy = Yc_theta*np.conjugate(Yc_theta) +  Yc_phi*np.conjugate(Yc_phi)
    JJr = []
    JJi = []
    JJr.append(xx.real)
    JJr.append(xy.real)
    JJr.append(yx.real)
    JJr.append(yy.real)
    JJi.append(xx.imag)
    JJi.append(xy.imag)
    JJi.append(yx.imag)
    JJi.append(yy.imag)
    JJr = np.array(JJr)
    JJi = np.array(JJi)
    '''pyfits.writeto('Jones_Real.fits',JJr[:,33,...],h[0:38], clobber=True)
    pyfits.writeto('Jones_Imag.fits',JJi[:,33,...],h[0:38], clobber=True)
    np.save('JJ_R', JJr[:,33,...]) 
    np.save('JJ_Im', JJi[:,33,...])'''
    np.save('JJ_Rr', JJr) 
    np.save('JJ_Img', JJi)
    ##############
    
    #
    #
    #fits_cube = pyfits.open("kat7_my_16_List.fits")-----	#
    #header = fits_cube[0].header --------------
    #
    fig = plt.figure(1, figsize=(12, 12), dpi=60)
    g, cax = setup_axes(fig, h, 2, 2)		#
    #
    # Plot Jones Real Images
    #
    cmap = plt.get_cmap('jet')	
    images = []
    for i, ax in enumerate(g):
        channel =  JJr[:,33,...][i]			#fits_cube[0].data[i]	#channel_number
        im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos") 	
        images.append(im)

    # make colorbar
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("Jy / Beam")
    plt.savefig('J_Real.png') 
    plt.show()
    #
    fig = plt.figure(1, figsize=(12, 12), dpi=60)
    g, cax = setup_axes(fig, h, 2, 2)
    # Plot Jones Imaginary Images
    #
    cmap = plt.get_cmap('jet')	
    images = []
    for i, ax in enumerate(g):
        channel = JJi[:,33,...][i]			#fits_cube[0].data[i]	#channel_number
        im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos") 	
        images.append(im)

    # make colorbar
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("Jy / Beam")
    plt.savefig('J_imag.png') 
    plt.show()
    #
    #Mueller_Matrices Conversion
    #
    print "\n >> Mueller Matrices Conversion:"
    #
    m_ii = 0.5*(xx*np.conjugate(xx) + xy*np.conjugate(xy) + yx*np.conjugate(yx)+ yy*np.conjugate(yy)) 
    m_iq = 0.5*(xx*np.conjugate(xx) - xy*np.conjugate(xy) + yx*np.conjugate(yx) - yy*np.conjugate(yy)) 
    m_iu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) + yx*np.conjugate(yy)+ yy*np.conjugate(yx)) 
    m_iv = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
    m_qi = 0.5*(xx*np.conjugate(xx) + xy*np.conjugate(xy) - yx*np.conjugate(yx) - yy*np.conjugate(yy)) 
    m_qq = 0.5*(xx*np.conjugate(xx) - xy*np.conjugate(xy) - yx*np.conjugate(yx)+ yy*np.conjugate(yy)) 
    m_qu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx)) 
    m_qv = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx)) 
    #m_ui = (xx*np.conjugate(yy) + yx*np.conjugate(xx) + xy*np.conjugate(yy) + yy*np.conjugate(xy))/2.
    m_ui = 0.5*(xx*np.conjugate(yx) + yx*np.conjugate(xx) + xy*np.conjugate(yy)+ yx*np.conjugate(xy)) 
    #m_uq = 0.5*(xy*np.conjugate(yx) + yx*np.conjugate(xx) - xy*np.conjugate(yy) - yy*np.conjugate(xy))
    #m_qu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx))
    m_uq = 0.5*(xx*np.conjugate(xy) - xy*np.conjugate(xx) + yx*np.conjugate(yy) + yy*np.conjugate(yx))  
    m_uu = 0.5*(xx*np.conjugate(yy) + yy*np.conjugate(xx) + xy*np.conjugate(yx)+ yx*np.conjugate(xy))
    m_uv = 0.5*1j*(xx*np.conjugate(yy) + yx*np.conjugate(xy) - yy*np.conjugate(xx) - xy*np.conjugate(yx))  
    #m_uv = 1j*(-xx*np.conjugate(yx) + xy*np.conjugate(yy) - yx*np.conjugate(xx) + yy*np.conjugate(xy))/2. 
    m_vi = 0.5*1j*(-xx*np.conjugate(yx) + yx*np.conjugate(xx) - xy*np.conjugate(yy) + yy*np.conjugate(xy)) 
    m_vq = 0.5*1j*(-xx*np.conjugate(yx) + yx*np.conjugate(xx) + xy*np.conjugate(yy) - yy*np.conjugate(xy)) 
    m_vu = 0.5*1j*(-xx*np.conjugate(yy) + yy*np.conjugate(xx) - xy*np.conjugate(yx) + yx*np.conjugate(xy)) 
    m_vv = 0.5*(xx*np.conjugate(yy) - yx*np.conjugate(xy) + yy*np.conjugate(xx) - xy*np.conjugate(yx)) 
    #
    M = []
    M.append(m_ii.real[33,...])
    M.append(m_iq.real[16,...])
    M.append(m_iu.real[33,...])
    M.append(m_iv.real[33,...])
    M.append(m_qi.real[16,...])
    M.append(m_qq.real[33,...])
    M.append(m_qu.real[33,...])
    M.append(m_qv.real[33,...])
    M.append(m_ui.real[33,...])
    M.append(m_uq.real[33,...])
    M.append(m_uu.real[33,...])
    M.append(m_uv.real[33,...])
    M.append(m_vi.real[33,...])
    M.append(m_vq.real[33,...])
    M.append(m_vu.real[33,...])
    M.append(m_vv.real[33,...])
    M = np.array(M)
    #pyfits.writeto('mueller.fits',M,h[0:38], clobber=True)
    pyfits.writeto('muellerfull.fits',M,h, clobber=True)
    #np.save('mueller',M[:,33,...])
    np.save('muellerF',M)
    #
    # Plot Mueller Images
    #
    fig = plt.figure(1, figsize=(12, 12), dpi=70)
    g, cax = setup_axes(fig, h, 4, 4)
    cmap = plt.get_cmap('jet')	
    images = []
    for i, ax in enumerate(g):
        channel = M[i]			#fits_cube[0].data[i]	#channel_number
        im = ax.imshow(channel, origin="lower", cmap=cmap, interpolation="lanczos") 	
        images.append(im)

    # make colorbar
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("Jy / Beam")
    plt.savefig('mueller.png') 
    plt.show()
    print 'Done !!!!!!!!!!!'
    # Oskar Interferometer Generation
    #
    print " ------------------------------------------------------"
    print "\n >> Checking for healpix map to run Oskar Interferometer:"
    print " ------------------------------------------------------"
    config = ConfigParser.ConfigParser()
    config.read('true_model_oskar_congifuration_setup1.ini')
    if config.get('sky', 'healpix_fits\map_units', 1)!= 'None':
       print "\n >> Oskar Interferometer Begins:"
       os.system('oskar_sim_interferometer true_model_oskar_congifuration_setup1.ini')
    else:
       print "\n >> Oskar Configuration has no healpix map:"
    #shutil.move('folder/station???', './')'''
        #
if __name__ == '__main__':
    run_test()
#	python aperture_array_model.py -n 10000 -d 12 -a 7 -t 1 
#../theo/comp_sc/emma_n_theo/plott/hlpx300.fits
