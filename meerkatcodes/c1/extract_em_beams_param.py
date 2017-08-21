#!/usr/bin/python
import pylab as plt
import gaussfitsutil as ut
import cPickle,numpy,sys,os,pyfits,scipy
import numpy as np

plotchan = [3]
clip = 0.2
#clip = 0.6

cachefile = 'beamfits-mearkat-em-clip-%s.cp'%clip

xx = {}
yy = {}


for ifname in ['em']:#['00', '012', '017'][:]: # 
#for ifname in ['000', '012', '017']:
    
    path = '/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/em-meerkatbeam/'

    fileReal = path + '%sodel_6deg_900_1100_1285_1350MHz_XX_re.fits'%ifname
    fileImag = path + '%sodel_6deg_900_1100_1285_1350MHz_XX_im.fits'%ifname
    xx[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    #cPickle.dump(xx[ifname],file('1487813282_m0%s_XX.cp'%ifname,'w'),2)
    #cPickle.dump(xx[ifname],file('1487813282_holobm2_m0%s_XX.cp'%ifname,'w'),2)
    
    fileReal = path + '%sodel_6deg_900_1100_1285_1350MHz_YY_re.fits'%ifname
    fileImag = path + '%sodel_6deg_900_1100_1285_1350MHz_YY_im.fits'%ifname
    yy[ifname]  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    #cPickle.dump(yy[ifname],file('1487813282_m%s_YY.cp'%ifname,'w'),2)
    #cPickle.dump(yy[ifname],file('1487813282_holobm2_m%s_YY.cp'%ifname,'w'),2)
    
cPickle.dump((xx,yy),file(cachefile,'w'))
print '\nDone !!!\n'
