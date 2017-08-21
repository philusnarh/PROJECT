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
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-clip-0.2.cp', 'rb') as input:
with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-holobm2-clip-0.6.cp', 'rb') as input:
        rr, ll = cPickle.load(input)
        
# ANTS = set(rr.keys()) - set([None])
ANTS = rr.keys()
print ANTS

print 'The antennas are : %s'%ANTS
#sys.exit()

srr,sll = cPickle.load(file('beamfits-mearkat-em-clip-0.2.cp'))  # sim LL/RR  parameters

sr,sl = srr['em'],sll['em']

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

holonumb = 2
path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm%s/'%holonumb
#path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm2/'
#path = './'


header2 = pyfits.open('/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/em-meerkatbeam/emodel_6deg_900_1100_1285_1350MHz_XX_re.fits')[0].header   # EM file header

terget_chan = -1  # 

for ant in ANTS[:] :
    
    hdu = {}
    newdata = {}
    
    print 'Antenna m%s'%ant

    for pol in ["XX","YX","XY","YY"][:]: # 'JVLA-L_interp_by_weight_freq-ll-im.fits'
        for reim in "re","im": 
            hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim))[0]
            newdata[pol,reim] = []#; newdata[pol,reim].append(hdu[pol,reim].data[0,:,:])

    header = hdu["XX","re"].header  # Holography file header
 
    crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]


    crpix11, crval11, cdelt11 = [ header2.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix22, crval22, cdelt22 = [ header2.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
    
    nchan, ypl, xpl = hdu["XX","re"].data.shape


    # calc. Holography beam centres in pixel coordinates    
    xcr1 = crpix1-1 + (rr[ant].x0.data - crval1)/cdelt1
    ycr1 = crpix2-1 + (rr[ant].y0.data - crval2)/cdelt2
    xcl1 = crpix1-1 + (ll[ant].x0.data - crval1)/cdelt1
    ycl1 = crpix2-1 + (ll[ant].y0.data - crval2)/cdelt2
    print sr, sr.xw
    # calc. EM beams centres in pixel coordinates    
    xcr2 = crpix11-1 + (sr.x0 - crval11)/cdelt11
    ycr2 = crpix22-1 + (sr.y0 - crval22)/cdelt22
    xcl2 = crpix11-1 + (sl.x0 - crval11)/cdelt11
    ycl2 = crpix22-1 + (sl.y0 - crval22)/cdelt22
    
    
    for ichan in range(nchan)[:]:


        # exclude frequency channals with high variance
        #if ~rr[ant].mask[ichan] and rr[ant].xw[ichan] != 0 and rr[ant].yw[ichan] != 0:

		#x = xcr[ichan] + sr.xw[ichan]*(numpy.arange(xpl) - xcr1[iantch])/rr[ant].xw[iantch]
		#y = ycr[ichan] + sr.yw[ichan]*(numpy.arange(ypl) - ycr1[iantch])/rr[ant].yw[iantch]
	x = xcr1[ichan] + rr[ant].xw.data[ichan]*(numpy.arange(xpl) - xcr2[terget_chan])/(sr.xw[terget_chan])
	y = ycr1[ichan] + rr[ant].yw.data[ichan]*(numpy.arange(ypl) - ycr2[terget_chan])/(sr.yw[terget_chan])

        for pol in ["XX","YX"][:]:
            for reim in ["re","im"]:
                data = hdu[pol,reim].data 
	        
                hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))


	x = xcl1[ichan] + ll[ant].xw.data[ichan]*(numpy.arange(xpl) - xcl2[terget_chan])/(sl.xw[terget_chan])
	y = ycl1[ichan] + ll[ant].yw.data[ichan]*(numpy.arange(ypl) - ycl2[terget_chan])/(sl.yw[terget_chan])


        for pol in ["YY","XY"][:]:
            for reim in ["re","im"]:
                data = hdu[pol,reim].data 
	        
                hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:]),numpy.array(numpy.meshgrid(y,x)))

    # Saving files to disk
    for pol in ["XX","YX","XY","YY"][:]:
        for reim in ['re','im']:
            path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb%s_scaled_down_to_em_grid_%s_%s.fits'%(ant,holonumb,pol,reim)
            hdu[pol,reim].writeto(path_fname,clobber=True)
#             pyfits.writeto("scaled_JVLA-L_interp_direct_from_sim_ant%s-%s_%s.fits"%(ant,pol,reim), newdata[pol,reim],
#                            header=header,clobber=True)
            #pyfits.writeto(path_fname, newdata[pol,reim], header=hdu[pol,reim].header,clobber=True)
           
            print "\t--- %s --- %s Done!"%(path_fname,hdu[pol,reim].data.shape) 
