#!/usr/bin/python
import pylab as plt
import gaussfitsutil as ut
import cPickle,numpy,sys,os,pyfits,scipy
import numpy as np

plotchan = [500]
clip = 0.2
#clip = 0.6

#cachefile = 'beamfits-mearkat-clip-%s.cp'%clip
cachefile = 'beamfits-mearkat-holobm2-clip-%s.cp'%clip
#cachefile = 'beamfits-mearkat_scaled-clip-%s.cp'%clip
#cachefile = 'beamfits-mearkat_scaled_down-clip-%s.cp'%clip
cachefile = 'beamfits-mearkat_mb1_scaled_down_to_em_grid-clip-%s.cp'%clip

xx = {}
yy = {}

#for ifname in ['00', '012', '017'][:]:
for ifname in ['000', '012', '017']:
    
    #path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm1/'
    path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm2/' # 1487813282_m012_900MHz_1MHz_500channels_Jones_scaled_down_YX_im.fits
    #path = './'
    #1487813282_m012_900MHz_1MHz_500channels_Jones_scaled_XX_im.fits
    fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_XX_re.fits'%ifname
    fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_XX_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_XX_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_XX_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_down_XX_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_down_XX_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb1_scaled_down_to_em_grid_XX_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb1_scaled_down_to_em_grid_XX_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_XX_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_XX_im.fits'%ifname
    xx[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    #cPickle.dump(xx[ifname],file('1487813282_m0%s_XX.cp'%ifname,'w'),2)
    #cPickle.dump(xx[ifname],file('1487813282_holobm2_m0%s_XX.cp'%ifname,'w'),2)
    
    fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_YY_re.fits'%ifname
    fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_YY_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_YY_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_YY_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_down_YY_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_scaled_down_YY_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb1_scaled_down_to_em_grid_YY_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb1_scaled_down_to_em_grid_YY_im.fits'%ifname
    #fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_YY_re.fits'%ifname
    #fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_nonans_YY_im.fits'%ifname
    yy[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    #cPickle.dump(yy[ifname],file('1487813282_m%s_YY.cp'%ifname,'w'),2)
    #cPickle.dump(yy[ifname],file('1487813282_holobm2_m%s_YY.cp'%ifname,'w'),2)
    
cPickle.dump((xx,yy),file(cachefile,'w'))
print '\nDone !!!\n'

sys.exit()

for ifname in ['000', '012', '017']:
#for ifname in ['00', '012', '017']:
    
    #path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm1/'
    path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/holobm2/'

    fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_XX_re.fits'%ifname
    fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_XX_im.fits'%ifname
    xx[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    #cPickle.dump(xx[ifname],file('1487813282_m0%s_XX.cp'%ifname,'w'),2)
    cPickle.dump(xx[ifname],file('1487813282_holobm2_m0%s_XX.cp'%ifname,'w'),2)
    
    fileReal = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_YY_re.fits'%ifname
    fileImag = path + '1487813282_m%s_900MHz_1MHz_500channels_Jones_YY_im.fits'%ifname
    yy[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    cPickle.dump(yy[ifname],file('1487813282_m%s_YY.cp'%ifname,'w'),2)
    #cPickle.dump(yy[ifname],file('1487813282_holobm2_m%s_YY.cp'%ifname,'w'),2)
    
cPickle.dump((xx,yy),file(cachefile,'w'))
print '\nDone !!!\n'

