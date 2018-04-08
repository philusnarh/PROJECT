""" This module provides some utility functions for beam investigation """
# ---------========= Accommodating frequency in 'exclude pixel variable =======----------

# gaussutil.py
from numpy import *
from scipy import optimize
import numpy as np
import numpy.ma
import traceback
from numpy.ma import masked_array

#import pylab as plt
import sys
import pyfits
import math
import os
import os.path
import pylab as plt


from numpy import *
from scipy import optimize
from scipy.ndimage.interpolation import map_coordinates 
from scipy import stats

def moments (data,circle,rotate,vheight,amplitude):
    """Returns (height, amplitude, x, y, width_x, width_y, rotation angle)
    the gaussian parameters of a 2D distribution by calculating its
    moments.  Depending on the input parameters, will only output 
    a subset of the above"""
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    width = ( width_x + width_y ) / 2.
    height = stats.mode(data.ravel())[0][0] if vheight else 0;
    if amplitude is None:
        amplitude = data.max()-height
        mylist = [amplitude,x,y]
    else:
        mylist = [x,y]
    if vheight==1:
        mylist = [height] + mylist
    if circle==0:
        mylist = mylist + [width_x,width_y]
    else:
        mylist = mylist + [width]
    if rotate==1:
        mylist = mylist + [0.] #rotation "moment" is just zero...
    return tuple(mylist)

def twodgaussian(inpars, circle, rotate, vheight, amplitude):
    """Returns a 2d gaussian function of the form:
        x' = cos(rota) * x - sin(rota) * y
        y' = sin(rota) * x + cos(rota) * y
        (rota should be in degrees)
        g = b + a exp ( - ( ((x-center_x)/width_x)**2 +
        ((y-center_y)/width_y)**2 ) / 2 )
        where x and y are the input parameters of the returned function,
        and all other parameters are specified by this function
        However, the above values are passed by list.  The list should be:
        inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)
        You can choose to ignore / neglect some of the above input parameters using the following options:
            circle=0 - default is an elliptical gaussian (different x, y widths), but can reduce
                the input by one parameter if it's a circular gaussian
            rotate=1 - default allows rotation of the gaussian ellipse.  Can remove last parameter
                by setting rotate=0
            vheight=1 - default allows a variable height-above-zero, i.e. an additive constant
                for the Gaussian function.  Can remove first parameter by setting this to 0
            amplitude=None: fit an amplitude, else set amplitude to a fixed value
        """
    inpars_old = inpars
    inpars = list(inpars)
    if vheight == 1:
        height = inpars.pop(0)
        height = float(height)
    else:
        height = float(0)
    if amplitude is None:
        amplitude = float(inpars.pop(0))
    amplitude = float(amplitude)
    center_x, center_y = inpars.pop(0),inpars.pop(0)
    center_x = float(center_x)
    center_y = float(center_y)
    if circle == 1:
        width = inpars.pop(0)
        width_x = float(width)
        width_y = float(width)
    else:
        width_x, width_y = inpars.pop(0),inpars.pop(0)
        width_x = float(width_x)
        width_y = float(width_y)
    if rotate == 1:
        rota = inpars.pop(0)
        rota = pi/180. * float(rota)
        rcen_x = center_x * cos(rota) - center_y * sin(rota)
        rcen_y = center_x * sin(rota) + center_y * cos(rota)
    else:
        rcen_x = center_x
        rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + " circle=%d, rotate=%d, vheight=%d" % (circle,rotate,vheight) )
            
    def rotgauss(x,y):
        if rotate==1:
            xp = x * cos(rota) - y * sin(rota)
            yp = x * sin(rota) + y * cos(rota)
        else:
            xp = x
            yp = y
        g = height+amplitude*exp(
            -(((rcen_x-xp)/width_x)**2+
            ((rcen_y-yp)/width_y)**2)/2.)
        return g
    return rotgauss

def gaussfit(data,err=None,params=[],autoderiv=1,return_all=0,circle=0,rotate=1,vheight=0,amplitude=None):
    """
    Gaussian fitter with the ability to fit a variety of different forms of 2-dimensional gaussian.
    
    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)
            if not input, these will be determined from the moments of the system, 
            assuming no rotation
        autoderiv=1 - use the autoderiv provided in the lmder.f function (the alternative
            is to us an analytic derivative with lmdif.f: this method is less robust)
        return_all=0 - Default is to return only the Gaussian parameters.  See below for
            detail on output
        circle=0 - default is an elliptical gaussian (different x, y widths), but can reduce
            the input by one parameter if it's a circular gaussian
        rotate=1 - default allows rotation of the gaussian ellipse.  Can remove last parameter
            by setting rotate=0
        vheight=1 - default allows a variable height-above-zero, i.e. an additive constant
            for the Gaussian function.  Can remove first parameter by setting this to 0
        amplitude=None - default fits peak amplitude, else set to a fixed value
    Output:
        Default output is a set of Gaussian parameters with the same shape as the input parameters
        Can also output the covariance matrix, 'infodict' that contains a lot more detail about
            the fit (see scipy.optimize.leastsq), and a message from leastsq telling what the exit
            status of the fitting routine was
        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    if params == []:
        params = (moments(data,circle,rotate,vheight,amplitude))
    if err == None:
        errfunc = lambda p: (twodgaussian(p,circle,rotate,vheight,amplitude)(*indices(data.shape)) - data)
    else:
        errfunc = lambda p: (twodgaussian(p,circle,rotate,vheight,amplitude)(*indices(data.shape)) - data)/err
    # if data is masked, modify the error function to ignore the masked pixels
    if numpy.ma.isMA(data):
        def errorfunction(p):
            res = errfunc(p)
            res[data.mask] = 0
            return np.ravel(res)
    else:
        def errorfunction(p):
            return np.ravel(errfunc(p))
    if autoderiv == 0:
        # the analytic derivative, while not terribly difficult, is less efficient and useful.  I only bothered
        # putting it here because I was instructed to do so for a class project - please ask if you would like 
        # this feature implemented
        raise ValueError("I'm sorry, I haven't implemented this feature yet.")
    else:
        p, cov, infodict, errmsg, success = optimize.leastsq(errorfunction, params, full_output=1)
    if  return_all == 0:
        return p
    elif return_all == 1:
        return p,cov,infodict,errmsg
# ============================================

def rotfit (data,x0,y0,radius=None):
    """
    Fits a rotation axis to the data.
    """
    nx, ny = data.shape
    def total_moment (angle):
        x, y = np.meshgrid(np.arange(nx)-x0, np.arange(ny)-y0)
        dist = x*math.cos(angle) - y*math.sin(angle)
        return abs((dist*data).sum())

    return optimize.minimize(total_moment,[0],options=dict(disp=True))



def getMaskedBeam(data, mask, data_re, threshold=0.25):
    ## OMS: Kela, note you can do this much more efficiently with boolean arrays
    ## your original function can be rewritten as one line:
    ##    data[data<=threshold] = 0
    ## it will also perform a lost faster!

    # #::===============::====================::

    # tsize = len(data.diagonal())
    # maskedsheet = np.zeros((tsize,tsize),dtype=int)
    # for i in range(tsize):
    #     for j in range(tsize):
    #         if data[i,j] > threshold:
    #             maskedsheet[i,j] = 1
    # maskedsheet = maskedsheet*data   
    
    # return maskedsheet
    # #::===============::====================::
    from scipy.ndimage.measurements import label
    nx,ny = data.shape
    # real part is either positive or negative in main lobe, set initial mask based on that
    region = data_re > 0 if data_re[nx/2,ny/2] > 0 else data_re < 0
    # add mask of points above threshold
    if threshold:
        region &= data > threshold
    # find connected regions, we want the center region only
    reglab,_ = label(region)
    # mask out all regions whose label is not equal to the centre pixel. Also mask high values
    mask1 = (reglab != reglab[nx/2,ny/2])
    # now find corresponding surrounding box
    x0 = sum(mask1[:nx/2,ny/2])   # number of 1 pixels in first half of centre row
    x1 = sum(mask1[nx/2:,ny/2])   # number of 1 pixels in second half of centre row
    y0 = sum(mask1[nx/2,:ny/2])   # number of 1 pixels in first half of centre column
    y1 = sum(mask1[nx/2,ny/2:])   # number of 1 pixels in second half of centre column
    # if box is too small, this beam is broken -- don't fit it
    if (nx-x0-x1)*(ny-y0-y1)<100 or (data.size - mask1.sum())<100:
        return None,x0,y0,reglab,mask1
    # return center box, along with x0,y0 coordinates
    mask |= mask1
    return masked_array(data[x0:-x1,y0:-y1],mask[x0:-x1,y0:-y1]),x0,y0,reglab,mask1


class FittedBeamParm (object):
    """This object represents the fitted parameters for one beam pattern.
    The attributes "peak0","h0","x00","y00","xw0","yw0","rot0","freq0"
    are arrays representing the raw unmasked parameters at all frequencies.
    The attribute "var" is the variance of the fit at every frequency.
    The attribute "mask" can be set via set_mask() to mask out bad channels 
    (based on variance, presumably)
    The attrbutes without a trailing 0 are then the masked versions of all the
    fit parameters.
    """
    def __init__ (self,nchan):
        for attr in "peak","h","x0","y0","xw","yw","rot","freq":
            z = np.zeros(nchan)
            setattr(self,attr+"0",z)
            setattr(self,attr,z)
        self.var = np.zeros(nchan)
        self.perr = [None]*nchan
        self.freq = None
        self.mask = None

    def set_mask (self,mask):
        self.mask = mask|(self.var<=0)
        self.var[self.var<=0] = 0
        for attr in "peak","h","x0","y0","xw","yw","rot","freq":
            setattr(self,attr,masked_array(getattr(self,attr+"0"),mask),)



def get2DGaussFitParams3(fitsfileReal='', fitsfileImag='',startchan=0,nchan=None,plotchan=set(),
                         clip=0,
                         rotate=1,circle=0,vheight=0,amplitude=None,fit_xform=False):
    """  This function returns the parameters of a 2-dimenssional Gaussian fit in 4 numpy array elements.
         It uses the Real beam to locat the center lobe then fits a 2D gaussian to it.

    Input Parameters:
        fitsfileReal - File name a fits file - antenna real file name eg. ''nt5RRreal.fits'
        fitsfileImag - File name a fits file - antenna imaginary file name eg. ''ant5RRimag.fits'
        startchan,nchan  - range of channels to process. Default is all
        plotchan    makes diagnostic plots for a set of channels
        clip        threshold to clip the main lobe at. If 0, clips main lobe using the real part>0

    Output:
        returns params[], err_in_args[]
    
        params is an array of peak[], x0[], y0[], sigma_x[], sigma_y[] in that order
        err_in_args[] is the 'params' return argurements 
    """
    
    FWHM = math.sqrt (math.log (256) )
    
    myfitfilefReal = pyfits.open(fitsfileReal) #creat a variable of the fits file
    myfitfilefImag = pyfits.open(fitsfileImag) #creat a variable of the fits file
    dataReal = myfitfilefReal[0].data.T # creatint a variable of the last image data
    dataImag = myfitfilefImag[0].data.T # creatint a variable of the last image data
    
    hdr = myfitfilefImag[0].header
    npix = hdr['NAXIS1']
    crvalx, crpixx, cdeltx = hdr['CRVAL1'], hdr['CRPIX1']-1, hdr['CDELT1']   # note ref pixel is 1-based in FITS so we subtract 1 
    crvaly, crpixy, cdelty = hdr['CRVAL2'], hdr['CRPIX2']-1, hdr['CDELT2']
    if abs(cdeltx) != abs(cdelty):
        raise RuntimeError,"sorry, we can't cope with non-square FITS file pixels yet"
    cell = abs(cdeltx)

    startFreqValue = hdr['CRVAL3'] #  =  1008000000.0 

    totchan = hdr['NAXIS3'] 
    if nchan is None:
        nchan = totchan - startchan
    else:
        nchan = min(nchan,totchan - startchan)

    ##OMS: commenting this out for now, let's use all channels in the fit and look at the errors
    # if fitsfile1[:3] == 'ant' :
    #     excludepixel = 32
    #     selectedChanls = getChannels() 
    # else:
    #     selectedChanls = range(selectedChanls)
    #     excludepixel = 192
        
    endchan = startchan + nchan           #myfitfilefImag[0].header['NAXIS3'] #  =  1024
    
    fbp = FittedBeamParm(totchan)

    ## compute an amplitude cube in one go
    ## transpose it so that axes are now x,y,frequency
    beamampl = np.sqrt( dataReal**2 + dataImag**2 )

    fbp.freq0 = fbp.freq = hdr['CRVAL3'] + (np.arange(startchan,endchan)-hdr['CRPIX3']-1)*hdr['CDELT3']

    for channel in range(endchan-1,startchan-1,-1):
        ##OMS: removing this, as getMaskedBeam should pull out the main lobe anyway
        ## making the exclude pixel frequency dependent - to accommodate the decrease in main as frequency increases

        beam0 = beamampl[:,:,channel] 
        # mask out high pixels (they're bad data)
        mask = beam0>2
        beam0[mask] = 0
    
        try:
            beam,xoff,yoff,region_labels,region_mask = getMaskedBeam(beam0,mask,dataReal[:,:,channel],clip)
            if beam is None:
                raise RuntimeError,"failed to find beam mask"
            params, cov, infodict, errmsg = gaussfit(beam,return_all=1,
                rotate=rotate,circle=circle,vheight=vheight,amplitude=amplitude)
            print channel,
        # if there are any errors in the fit, plot the offending beam,mark it with variance=-1,
        # and go on to next channel
        except:
            print "Failed to fit beam in channel",channel
            # traceback.print_exc()
            # plt.figure(figsize=(30,10))
            # plt.subplot(1,3,1)
            # plt.imshow(beam0)
            # plt.colorbar()
            # plt.title("Failed to fit beam in channel %d"%channel)  
            # plt.subplot(1,3,2)
            # plt.title("Region labels")  
            # plt.imshow(region_labels)
            # plt.colorbar()
            # plt.subplot(1,3,3)
            # plt.imshow(region_mask)
            # plt.title("Region mask")  
            # plt.colorbar()
            fbp.var[channel] = 0
            continue

        residual = infodict['fvec'].reshape(beam.shape)
        fbp.var[channel] = masked_array(residual,beam.mask).std()
        
        # fit parameters returned in this order, but some are optional depending on
        # arguments
        h = peak_0 = x0 = y0 = sigma_x = sigma_y = rot = 0
        params = list(params)
        # so we pop them from the head of list one by one
        h = params.pop(0) if vheight else 0
        peak_0 = params.pop(0) if amplitude is None else amplitude
        x0, y0 = params.pop(0), params.pop(0)
        if circle:
            sigma_x = sigma_y = params.pop(0)
        else:
            sigma_x,sigma_y = params.pop(0), params.pop(0)
        rot = params.pop(0) if rotate else 0
        # sanity check -- all params should be accounted for now
        if params:
            raise RuntimeError,"oops, some parameters left over. This is a bug"

        # since getMaskedBeam cut out a box, add coordinates of corner of box to fitted x0, y0
        # to get the "true" centre coordinate in pixels
        x0 += xoff
        y0 += yoff

        # experimental: try to fit transform to map to highest channel
        if fit_xform:
            if channel == endchan-1:
                beam0ref,x0ref,y0ref,xwref,ywref = beam0, x0, y0, sigma_x, sigma_y
                xgridref, ygridref = numpy.arange(npix)-x0ref, numpy.arange(npix)-y0ref
            else:
                def xform (beam,params,verbose=False):
                    rot, xw1, xw2 = params[:3]
                    yw1, yw2 = (xw1, xw2) if circle else params[4:5]
                    r = numpy.sqrt(xgridref**2 + ygridref**2)
                    x1 = xw1*(xw2**(r/100))*xgridref
                    y1 = yw1*(yw2**(r/100))*ygridref
                    c, s = math.cos(rot), math.sin(rot)
                    x2 =  x0 + x1*c + y1*s
                    y2 =  y0 - x1*s + y1*c
                    if verbose:
                        print [int(x) for x in x2]
                        print [int(y) for y in y2]
                    return map_coordinates(beam,numpy.array(numpy.meshgrid(x2,y2,indexing='ij')),cval=-.01)
                def rotform (beam,params,verbose=False):
                    rot = params[0]
                    x1 = sigma_x/xwref*xgridref
                    y1 = sigma_y/ywref*ygridref
                    c, s = math.cos(rot), math.sin(rot)
                    x2 =  x0 + x1*c + y1*s
                    y2 =  y0 - x1*s + y1*c
                    if verbose:
                        print [int(x) for x in x2]
                        print [int(y) for y in y2]
                    return map_coordinates(beam,numpy.array(numpy.meshgrid(x2,y2,indexing='ij')),cval=-.01)

                # errfunc translates beam into beam0 frame and returns the difference
                def errfunc (params):
                    beam1 = xform(beam,params)
                    diff = beam1 - beam0ref
                    diff[beam1<0] = 0
                    return np.ravel(diff)
                xpars0 = [0, sigma_x/xwref, 2]
                if not circle:
                    xpars0 += [sigma_y/ywref, 1]
                xpars, xcov, xinfodict, xerrmsg, xsuccess = optimize.leastsq(errfunc, xpars0, full_output=1)
                print "Fitted transform, channel",channel,xpars,xpars[0]*180/math.pi
                beam1a = xform(beam0,xpars0)
                beam1b = xform(beam0,xpars,verbose=True)
                # plt.figure(figsize=(20,5))
                # plt.subplot(141)
                # plt.imshow(beam0ref)
                # plt.subplot(142)
                # plt.imshow(beam1a)
                # plt.subplot(143)
                # plt.imshow(beam1b)
                # plt.subplot(144)
                # diff = beam1b - beam0ref
                # diff[beam1b<0] = 0
                # plt.imshow(diff)
                # plt.show()

        fbp.peak0[channel] = peak_0 
        fbp.h0[channel] = h
        fbp.rot0[channel] = rot

        # now convert to world coordinates in fits
        fbp.x00[channel] = x0 = crvalx + (x0 - crpixx)*cdeltx
        fbp.y00[channel] = y0 = crvaly + (y0 - crpixy)*cdelty
        fbp.xw0[channel] = width_x = abs(FWHM*sigma_x*cdeltx)
        fbp.yw0[channel]  = width_y = abs(FWHM*sigma_y*cdeltx)

        if cov is None:
            print "channel",channel,"of",fitsfileReal,fitsfileImag,": no covariance"
        else:
            fbp.perr[channel] = np.sqrt(cov.diagonal()*cell)/(1.0*beam.size)
        # if cov is None or (plotchan and channel in plotchan):
        #     print fitsfileReal,fitsfileImag
        #     ## beam = np.transpose(beam) no need for this since we already transpose above
        #     plt.figure(figsize=(15,5))
        #     plt.subplot(131)
        #     plt.imshow(beam)
        #     plt.colorbar()
        #     plt.title("Masked beam for channel %d"%channel)  
        #     print "channel",channel,beam.shape,xoff,yoff
        #     print h,peak_0,x0,y0,sigma_x,sigma_y,rot
        #     print fbp.perr[channel],"variance",fbp.var[channel]
        #     plt.subplot(132)
        #     plt.imshow(residual)
        #     plt.title("Residual of fit")  
        #     plt.colorbar()
        #     plt.subplot(133)
        #     plt.imshow(region_labels)
        #     plt.title("Regions for masking")  
        #     plt.colorbar()
        #     plt.show()

    #return peak, centerx, centery, sizecenterlobehorzontal, sizecenterlobevertical   #, covs, infodicts
    return fbp   #, covs, infodicts

