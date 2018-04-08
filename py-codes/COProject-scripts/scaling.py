#!/usr/bin/python
import cPickle,numpy,sys,os,pyfits
import numpy as np
from numpy.ma import masked_array
import pylab as plt
import matplotlib.pylab as plt
import gaussfitsutil as ut
from numpy.linalg import svd
import datetime

from scipy.ndimage.interpolation import map_coordinates
import warnings
import traceback
#warnings.filterwarnings('ignore')

from utilities import *
from gaussfitsutil import *

def extract_params(pref, _clip=0.2):
    s = time.time()
    clip = _clip # Used to maskout very high pixel values 

    cachefile = '%sfbp_clip%s.cp'%(pref,clip)

    data = fits.getdata(pref+'XX_R.fits')
    nchan = data.shape[0]
    
    xx = get2DGaussFitParams3(pref+'XX_R.fits', pref+'XX_I.fits', clip=clip)
    yy = get2DGaussFitParams3(pref+'YY_R.fits', pref+'YY_I.fits', clip=clip)

    cPickle.dump((xx,yy),file(cachefile,'w'))

    e = time.time()
    print "Time taken %.2f min"%((e-s)/60.)

def scale_size(pref, fbp, ref_chan=494, var_thresh=0.3):
    start = time.time()
    rr,ll = cPickle.load(file(fbp))
    header = fits.getheader(pref+'XX_R.fits')
    data = fits.getdata(pref+'XX_R.fits')
    data_scaled = np.zeros(data.shape)

    # Flag on high variance
    rr.set_mask((rr.var>var_thresh)|(rr.var<=0))
    ll.set_mask((ll.var>var_thresh)|(ll.var<=0))
    rr.var[rr.var<=0] = 0
    ll.var[ll.var<=0] = 0
 
    crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
    cdelt3, crval3 = [ header.get(x) for x in "CDELT3", "CRVAL3" ]

    nchan, ypl, xpl = data.shape

    # calculate beam centres in pixel coordinates
    xcr = crpix1-1 + (rr.x0.data - crval1)/cdelt1
    ycr = crpix2-1 + (rr.y0.data - crval2)/cdelt2
    xcl = crpix1-1 + (ll.x0.data - crval1)/cdelt1
    ycl = crpix2-1 + (ll.y0.data - crval2)/cdelt2

    print xcr, ycr, xcl, ycl
    return

    P, C = ['R', 'I'], ['X', 'Y']

    for p in range(2):
        for i in range(2):
            for j in range(2):
                filename = pref+'%s%s_%s'%(C[i],C[j],P[p])
                print '%s%s_%s'%(C[i],C[j],P[p])
                data = fits.getdata(filename+'.fits')

                for ichan in range(nchan):
                    if j==0:
                        x = xcr[ichan] + rr.xw.data[ichan]*(np.arange(xpl) - xcr[ref_chan])/(rr.xw[ref_chan])
                        y = ycr[ichan] + rr.yw.data[ichan]*(np.arange(ypl) - ycr[ref_chan])/(rr.yw[ref_chan])
                    elif j==1:
                        x = xcl[ichan] + ll.xw.data[ichan]*(np.arange(xpl) - xcl[ref_chan])/(ll.xw[ref_chan])
                        y = ycl[ichan] + ll.yw.data[ichan]*(np.arange(ypl) - ycl[ref_chan])/(ll.yw[ref_chan])

                    data_scaled[ichan,:,:] = map_coordinates(data[ichan,:,:], np.array(np.meshgrid(y,x)))

                fits.writeto(filename+'_scaled.fits', data_scaled, header, clobber=True)

    end = time.time()
    print "... Time taken: %.2f minutes"%((end-start)/60.)

#path = '/home/asad/data/meerkat/beam/holography/'
#pref = path+'1487813282_m017_256px_856MHz_1MHz_857channels_Jones_'
#scale_size(pref, pref+'fbp_clip0.2.cp', 820, 0.3)

def genBeamsParameters(filenamePrefix, filePol=['XX','YX', 'XY', 'YY'], fileReim=['re','im'], _clip=0.2, fileInputPath='./', outPath='./'):   


    clip = _clip # Used to maskout very high pixel values 

    cachefile = 'beamfits-%s-clip-%s.cp'%(filenamePrefix,clip) 


    fileReal = '%s%s_%s_%s.fits'%(fileInputPath,filenamePrefix,filePol[0],fileReim[0])
    fileImag = '%s%s_%s_%s.fits'%(fileInputPath,filenamePrefix,filePol[0],fileReim[1])

    data = pyfits.open(fileReal)[0].data
    plotchan = [data.shape[0]] # [nchans] -  number of frequencies
    print 'generating beams parameter for %s beams...\n'%filePol[0]
    xx  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    print '....Done!\n'
 
    fileReal = '%s%s_%s_%s.fits'%(fileInputPath,filenamePrefix,filePol[-1],fileReim[0])
    fileImag = '%s%s_%s_%s.fits'%(fileInputPath,filenamePrefix,filePol[-1],fileReim[1])
    print 'generating beams parameter for %s beams...\n'%filePol[-1]
    yy  = ut.get2DGaussFitParams3(fileReal,fileImag,plotchan=plotchan,clip=clip)
    print '....Done!\n'

        
    print 'Saving beams parameters to disk....\n'
    cPickle.dump((xx,yy),file('%s%s'%(outPath,cachefile),'w'))
    print 'Beams\' parameters saved in \'%s\'\nDone !!!\n'%cachefile

    return cachefile,(xx,yy)


def scaleBeansToRefBeam(filenamePrefix, refFilePrefix, picklefilename, refpicklefilename, filePol=['XX','YX', 'XY', 'YY'], 
                        refFilePol=['XX','YX', 'XY', 'YY'], fileReim=['re','im'], refFileReim=['re','im'], tergetFreq_chan=0,
                        var_threshold=0.301, fileInputPath='./', refFileInputPath='./', outPath='./', picklefileInputPath='./', relfPickleInputPath='./' ):   


    print 'Loading parameters from \'%s\' .....'%(picklefileInputPath + picklefilename)
    rr,ll = cPickle.load(file(picklefileInputPath + picklefilename)) #'beamfits-mearkat-clip-0.2.cp'))         
    print '.....Done\n'

    print 'Loading parameters from \'%s\' .....'%(relfPickleInputPath + refpicklefilename)
    sr,sl = cPickle.load(file(relfPickleInputPath + refpicklefilename)) #'beamfits-mearkat-em-clip-0.2.cp'))  # sim LL/RR  parameters
    print '.....Done\n'


    fname = '%s_%s_%s.fits'%(refFileInputPath + refFilePrefix, refFilePol[0], refFileReim[0])
    print 'Loading reference beam file\'s \'%s\' header .....'%fname
    #header2 = pyfits.open('/net/elwood/home/narh/MEERKAT_HOLOBEAMS/meerkatholo/em-meerkatbeam/emodel_6deg_900_1100_1285_1350MHz_XX_re.fits')[0].header   # EM file header
    header2 = pyfits.open(fname)[0].header   # EM file header
    print '.....Done\n'

    print 'The terget frequency is \'%s all beams will be scaled to the size of this beam'%tergetFreq_chan
    ref_chan = tergetFreq_chan  # 



    #sr,sl = srr[sll.keys()[0]],sll[sll.keys()[0]]

    #ANTS = rr.keys()
    #print ANTS
    #print 'The antennas are : %s'%ANTS


    maskout_nans = lambda a : np.nan_to_num(a) 

    #Flag on high variance
    # this is the threshold at which high variance is flagged

    print 'Variance threshold = %s\n'%var_threshold

    VAR_THRESHOLD = var_threshold #0.301

    rr.set_mask((rr.var>VAR_THRESHOLD)|(rr.var<=0))
    ll.set_mask((ll.var>VAR_THRESHOLD)|(ll.var<=0))
    rr.var[rr.var<=0] = 0
    ll.var[ll.var<=0] = 0
        


        
    hdu = {}
    newdata = {}
    
    #print 'Antenna m%s'%ant

    for pol in filePol: #["XX","YX","XY","YY"][:]: # 'JVLA-L_interp_by_weight_freq-ll-im.fits'
        for reim in fileReim: #"re","im": 
            #hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim))[0]
            fname = '%s_%s_%s.fits'%(fileInputPath + filenamePrefix, pol, reim)
            print 'Loading beam file \'%s\' .....'%fname
            #hdu[pol,reim] = pyfits.open(path+'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim))[0]
            hdu[pol,reim] = pyfits.open(fname)[0]
            newdata[pol,reim] = []#; newdata[pol,reim].append(hdu[pol,reim].data[0,:,:])
            print '.....Done\n'

   
    header = hdu[filePol[0], fileReim[0]].header  # Holography file header
 
    crpix1, crval1, cdelt1 = [ header.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
    cdelt3, crval3 = [ header.get(x) for x in "CDELT3", "CRVAL3" ] 


    crpix11, crval11, cdelt11 = [ header2.get(x) for x in "CRPIX1", "CRVAL1", "CDELT1" ]
    crpix22, crval22, cdelt22 = [ header2.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]

    #emfreq = [header2.get("CRVAL3")+i*header2.get("CDELT3") for i in range(header2.get("NAXIS3"))]
    #freq = [header.get("CRVAL3")+i*header.get("CDELT3") for i in range(header.get("NAXIS3"))]
    #freq = [header.get("CRVAL3")+i*1e+6 for i in range(header.get("NAXIS3"))]

     
    nchan, ypl, xpl = hdu[filePol[0], fileReim[0]].data.shape

    # calc. Holography beam centres in pixel coordinates    
    xcr1 = crpix1-1 + (rr.x0.data - crval1)/cdelt1
    ycr1 = crpix2-1 + (rr.y0.data - crval2)/cdelt2
    xcl1 = crpix1-1 + (ll.x0.data - crval1)/cdelt1
    ycl1 = crpix2-1 + (ll.y0.data - crval2)/cdelt2
    #print sr, sr.xw
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
    	#x = xcr1[ichan] + rr[ant].xw.data[ichan]*(numpy.arange(xpl) - xcr2[ref_chan])/(sr.xw[ref_chan])
        #print xcr1.shape, rr.xw.data[ichan].shape,numpy.arange(xpl).shape, xcr2[ref_chan].shape,sr.xw[ref_chan].shape
    	x = xcr1[ichan] + rr.xw.data[ichan]*(numpy.arange(xpl) - xcr2[ref_chan])/(sr.xw[ref_chan])
    	y = ycr1[ichan] + rr.yw.data[ichan]*(numpy.arange(ypl) - ycr2[ref_chan])/(sr.yw[ref_chan])

    	for pol in filePol[:2] : #["XX","YX"][:]: 
    	    for reim in fileReim: # ["re","im"]:
    	        data = hdu[pol,reim].data 
    		
    	        hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:].T),numpy.array(numpy.meshgrid(y,x)))
        # else:
        #     for pol in ["XX","YX"]:
        #         for reim in ["re","im"]:
        #             hdu[pol,reim].data[ichan,:,:] = 0
                    


        #if ~ll[ant].mask[ichan] and ll[ant].xw[ichan] != 0 and ll[ant].yw[ichan] != 0:
    	x = xcl1[ichan] + ll.xw.data[ichan]*(numpy.arange(xpl) - xcl2[ref_chan])/(sl.xw[ref_chan])
    	y = ycl1[ichan] + ll.yw.data[ichan]*(numpy.arange(ypl) - ycl2[ref_chan])/(sl.yw[ref_chan])
    	

    	for pol in filePol[2:] : #["YY","XY"][:]:
    	    for reim in fileReim: # ["re","im"]:
    	        data = hdu[pol,reim].data 
    		
    	        hdu[pol,reim].data[ichan,:,:] = map_coordinates(maskout_nans(data[ichan,:,:].T),
                                                                numpy.array(numpy.meshgrid(y,x)))
        # else:
        #     for pol in ["YY","XY"]:
        #         for reim in ["re","im"]:
        #             hdu[pol,reim].data[ichan,:,:] = 0
                    

    # Saving files to disk
    for pol in filePol[:] : #["XX","YX","XY","YY"][:]:
        for reim in fileReim: # ['re','im']: #'1487813282_m%s_900MHz_1MHz_500channels_Jones_%s_%s.fits'%(ant,pol,reim)
            #path_fname = '1487813282_m%s_900MHz_1MHz_500channels_Jones_mb%s_scaled_up_to_em_grid_%s_%s.fits'%(ant,holonumb,pol,reim)  
            fname = 'scaled_to_chan%s_%s_%s_%s.fits'%(ref_chan,filenamePrefix, pol, reim)

            hdu[pol,reim].writeto(outPath + fname,clobber=True)

            #pyfits.writeto(path_fname, newdata[pol,reim], header=hdu[pol,reim].header,clobber=True)
           
            print "\t--- %s --- %s Done!"%(fname,hdu[pol,reim].data.shape) 
    print '\t\t---------------===== Done scaling all beams! =====-----------------'





def rescale_data_size(data, newsizex, newsizey):
# ++++++++++ interpolating to have same size ++++++++
    dshape = data.shape
    # define new size
    outKSize_x = newsizex
    outKSize_y = newsizey
    
    # Rescale Data Size
    x_old = np.linspace(-dshape[0]/2., dshape[0]/2., dshape[0])      
    y_old = np.linspace(-dshape[-1]/2., dshape[-1]/2., dshape[-1])
    xnew = np.linspace(x_old.min(), x_old.max(), outKSize_x)
    ynew =  np.linspace(y_old.min(), y_old.max(), outKSize_y)
    
    # Perform Interpolation
    interp_Fxn = interpolate.RectBivariateSpline(np.sort(x_old),
                                 np.sort(y_old),data, kx=3,ky=3)           
    return interp_Fxn(xnew,ynew)


outPutString = lambda s: '%s %s %s'%('-'*9+'='*4,s,'='*4+'-'*9)

def calc_coverienc_matrix(_z):

    # ==========================================================
    print '\n\t%s Computing the covarience matrix... %s'%('-'*5+'='*3,'='*3+'-'*5)
    _z = masked_array(_z)
    m,n = _z.shape
    # rZ = np.array([(Z[i]/Z.sum(axis=1)[i]) for i in range(m)])
    rZ = np.array([_z[i]/_z[i].sum() for i in range(m)])
    # rZ = Z/np.array([Z.sum(axis=1)]*n).T
    mean = masked_array(rZ).mean(axis=0)
    rZ = rZ - masked_array(rZ).mean(axis=0)
    rZTrZ = np.dot(rZ.T,rZ)/(1.0*m)

    rZ = 0  # Added to free up memory
    print outPutString('done!!')
    return rZTrZ

# ========================================================== """
