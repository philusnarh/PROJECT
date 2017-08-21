#!/usr/bin/python
import cPickle,numpy,sys,os,pyfits
import numpy as np
from numpy.ma import masked_array
import pylab as plt
import matplotlib.pylab as plt
import gaussfitsutil
from scipy.optimize import curve_fit

from scipy.ndimage.interpolation import map_coordinates
from copy import deepcopy as deepcp
from scipy.ndimage.interpolation import shift
import warnings
import traceback
warnings.filterwarnings('ignore')

        
# with open('/home/iheanetu/images/beamfits-masked-norot-circ-clip0.20.cp', 'rb') as input:
#         rr, ll = cPickle.load(input)
with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-clip-0.2.cp', 'rb') as input:
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-holobm2-clip-0.6.cp', 'rb') as input:
        rr, ll = cPickle.load(input)
        
# ANTS = set(rr.keys()) - set([None])
ANTS = rr.keys()

print 'The antennas are : %s'%ANTS
#sys.exit()
# sr,sl = rr[None], ll[None]  # sim LL/RR  parameters

maskout_nans = lambda a : np.nan_to_num(a) # np.ma.masked_array(a, np.isnan(a)) # np.ma.masked_where(np.isnan(a), a, copy=True)

#Flag on high variance
# this is the threshold at which high variance is flagged

VAR_THRESHOLD = 0.01
for iant,ant in enumerate(ANTS):
    rr[ant].set_mask((rr[ant].var>VAR_THRESHOLD)|(rr[ant].var<=0))
    ll[ant].set_mask((ll[ant].var>VAR_THRESHOLD)|(ll[ant].var<=0))
    rr[ant].var[rr[ant].var<=0] = 0
    ll[ant].var[ll[ant].var<=0] = 0
    
        
# compute mean values of everything across antennas

path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm1/'
#path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm2/'
path = './'

terget_chan = -1  # 

for ant in ANTS[:] :
    
    hdu = {}
    newdata = {}
    
    print 'Antenna m%s'%ant

    for pol in ["XX","YX","XY","YY"][:1]: # 'JVLA-L_interp_by_weight_freq-ll-im.fits'
        for reim in "re","im": 
            #hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim))[0]
            hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_%s_%s.fits'%(ant,pol,reim))[0]
            newdata[pol,reim] = []#; newdata[pol,reim].append(hdu[pol,reim].data[0,:,:])

    header = hdu["XX","re"].header
 
    crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
    
    nchan, ypl, xpl = hdu["XX","re"].data.shape


    # Holography beam centres in pixel coordinates    
    xcr1 = crpix1-1 + (rr[ant].x0.data - crval1)/cdelt1
    ycr1 = crpix2-1 + (rr[ant].y0.data - crval2)/cdelt2
    xcl1 = crpix1-1 + (ll[ant].x0.data - crval1)/cdelt1
    ycl1 = crpix2-1 + (ll[ant].y0.data - crval2)/cdelt2

    # Simulated cassbeam centres in pixel coordinates calc for exponential fit   
#     xcr2 = crpix1-1 + (fit_sr_x0 - crval1)/cdelt1
#     ycr2 = crpix2-1 + (fit_sr_y0 - crval2)/cdelt2
#     xcl2 = crpix1-1 + (fit_sl_x0 - crval1)/cdelt1
#     ycl2 = crpix2-1 + (fit_sl_y0 - crval2)/cdelt2
    
    
    for ichan in range(nchan)[:]:


        # exclude frequency channals with high variance
        #if ~rr[ant].mask[ichan] and rr[ant].xw[ichan] != 0 and rr[ant].yw[ichan] != 0:
        if rr[ant].xw.data[ichan] != 0 and rr[ant].yw.data[ichan] != 0:

		#x = xcr[ichan] + sr.xw[ichan]*(numpy.arange(xpl) - xcr1[iantch])/rr[ant].xw[iantch]
		#y = ycr[ichan] + sr.yw[ichan]*(numpy.arange(ypl) - ycr1[iantch])/rr[ant].yw[iantch]
	    x = xcr1[ichan] + rr[ant].xw.data[ichan]*(numpy.arange(xpl) - xcr1[terget_chan])/(rr[ant].xw.data[terget_chan]*1.0)
	    y = ycr1[ichan] + rr[ant].yw.data[ichan]*(numpy.arange(ypl) - ycr1[terget_chan])/(rr[ant].yw.data[terget_chan]*1.0)
	    #print '\t\t\t\t\t\t\t--------====== Okay for Ant %s Chan%s=====---------'%(ant,ichan)
	#                 x = xcr[isimchan+1] + sr.xw[isimchan+1]*(numpy.arange(xpl) - xcr1[iantch])/rr[ant].xw[iantch]
	#                 y = ycr[isimchan+1] + sr.yw[isimchan+1]*(numpy.arange(ypl) - ycr1[iantch])/rr[ant].yw[iantch]

            for pol in ["XX","YX"][:1]:
                for reim in ["re","im"]:
                    data = hdu[pol,reim].data 
		    
#                             dat = map_coordinates(data[isimchan+1,:,:],numpy.array(numpy.meshgrid(y,x))).T newdata[pol,reim]
                    hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))#;print hdu[pol,reim].data[ichan,150:256,150:256];sys.exit()
                    #hdu[pol,reim].data[ichan,:,:] = map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x)))
                    #newdata[pol,reim].append(map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x))))

        else: # rr[ant].xw[ichan] != 0 and rr[ant].yw[ichan] != 0:
	    x = xcr1[ichan] + rr[ant].xw.data[ichan]*(numpy.arange(xpl) - xcr1[terget_chan])/(rr[ant].xw.data[terget_chan]*1.0)
	    y = ycr1[ichan] + rr[ant].yw.data[ichan]*(numpy.arange(ypl) - ycr1[terget_chan])/(rr[ant].yw.data[terget_chan]*1.0)
  	    #print 'No for Ant %s Chan%s'%(ant,ichan)
            for pol in ["XX","YX"][:1]:
                for reim in ["re","im"]:
                    #newdata[pol,reim].append(hdu[pol,reim].data[ichan,:,:])
                    data = hdu[pol,reim].data ;#print hdu[pol,reim].data[ichan,150:256,150:256]
                    #hdu[pol,reim].data[ichan,:,:] = map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x)))
                    hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))#;print '\n', rr[ant].xw.data[ichan],rr[ant].xw.data[-1], rr[ant].xw[ichan]/rr[ant].xw[-1],'\n' ;print hdu[pol,reim].data[ichan,150:256,150:256];sys.exit()
                    #newdata[pol,reim].append(map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x))))

        #if ~ll[ant].mask[ichan] and ll[ant].xw[ichan] != 0 and ll[ant].yw[ichan] != 0:
        if ll[ant].xw.data[ichan] != 0 and ll[ant].yw.data[ichan] != 0:

            # convert pixels on high-freq grid into fractional coordinates on low-freq grid
                 #x = xcl[isimchan+1] + sl.xw[isimchan+1]*(numpy.arange(xpl) - xcl1[iantch])/ll[ant].xw[iantch]
                # y = ycl[isimchan+1] + sl.yw[isimchan+1]*(numpy.arange(ypl) - ycl1[iantch])/ll[ant].yw[iantch]
            x = xcl1[ichan] + ll[ant].xw.data[ichan]*(numpy.arange(xpl) - xcl1[terget_chan])/(ll[ant].xw.data[terget_chan]*1.0)
            y = ycl1[ichan] + ll[ant].yw.data[ichan]*(numpy.arange(ypl) - ycl1[terget_chan])/(ll[ant].yw.data[terget_chan]*1.0)


            for pol in ["YY","XY"][:0]:
                for reim in ["re","im"]:
                    data = hdu[pol,reim].data ;#print hdu[pol,reim].data[ichan,150:256,150:256]
                    #hdu[pol,reim].data[ichan,:,:] = map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x))) 
		    hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))
                    #newdata[pol,reim].append(map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x))))

        else: # rr[ant].xw[ichan] != 0 and rr[ant].yw[ichan] != 0:
            x = xcl1[ichan] + ll[ant].xw.data[ichan]*(numpy.arange(xpl) - xcl1[terget_chan])/(ll[ant].xw.data[terget_chan]*1.0)
            y = ycl1[ichan] + ll[ant].yw.data[ichan]*(numpy.arange(ypl) - ycl1[terget_chan])/(ll[ant].yw.data[terget_chan]*1.0)
            for pol in ["YY","XY"][:0]:
                for reim in ["re","im"]:
                    data = hdu[pol,reim].data 
                    #newdata[pol,reim].append(hdu[pol,reim].data[ichan,:,:])
                    #hdu[pol,reim].data[ichan,:,:] = map_coordinates(data[ichan,:,:],numpy.array(numpy.meshgrid(y,x))) 
		    hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))



    for pol in ["XX","YX","XY","YY"][:1]:
#     for pol in ["rr","ll"][:]:
        for reim in ['re','im']:
            path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_%s_%s.fits'%(ant,pol,reim)
            #path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_down_%s_%s.fits'%(ant,pol,reim)
            #path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_holobm2_scaled_%s_%s.fits'%(ant,pol,reim)
            hdu[pol,reim].writeto(path_fname,clobber=True)
#             pyfits.writeto("scaled_JVLA-L_interp_direct_from_sim_ant%s-%s_%s.fits"%(ant,pol,reim), newdata[pol,reim],
#                            header=header,clobber=True)
            #pyfits.writeto(path_fname, newdata[pol,reim], header=hdu[pol,reim].header,clobber=True)
           
            print "\t--- %s --- %s Done!"%(path_fname,hdu[pol,reim].data.shape) 
