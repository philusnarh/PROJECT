#!/usr/bin/python
import pylab as plt
import gaussfitsutil as ut
import cPickle,numpy,sys,os,pyfits,scipy
import numpy as np

maskout_nans = lambda a : np.nan_to_num(a)

path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm1/'
#path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm2/'

#sys.exit()
ANTS = ['00', '012', '017']
#ANTS = ['000', '012', '017']

for ant in ANTS[:] :
    
    hdu = {}
    newdata = {}
    
    print 'Antenna m%s'%ant

    for pol in ["XX","YX","XY","YY"][:]:
        for reim in "re","im": 
            hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim))[0]

    nchan, _,_ = hdu[pol,reim].data.shape
    for ichan in range(nchan)[:]:


        for pol in ["XX","YX","XY","YY"][:]:
            for reim in ["re","im"]:
                data = hdu[pol,reim].data
                hdu[pol,reim].data[ichan,:,:] = maskout_nans(data[ichan,:,:])


    for pol in ["XX","YX","XY","YY"][:]:
        for reim in ['re','im']:
            path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_%s_%s.fits'%(ant,pol,reim)
            hdu[pol,reim].writeto(path_fname,clobber=True)
           
            print "\t--- %s --- %s Done!"%(path_fname,hdu[pol,reim].data.shape) 

