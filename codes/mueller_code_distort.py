#!/usr/bin/env python
#
import pyfits #as pf
import numpy as np
from matplotlib import pylab as plt
#
hduA = pyfits.open('beam_pattern_distortedGP_VOLTAGE.fits') 
#
hduPh = pyfits.open('beam_pattern_distortedGP_PHASE.fits') 
#
datA = hduA[0].data
print datA.shape 
#
datPh = hduPh[0].data 
#
#converting nan & inf
#
datA = np.where((datA == 0.) | (np.isnan(datA)) |(np.isinf(datA)), 0.000012, datA) 
datPh = np.where((datPh == 0.) | (np.isnan(datPh)) | (np.isinf(datA)), 0.000012, datPh) 
#
#Extracting the Amplitude
#
X_Amp_theta = datA[0, :, 0, ...] 
X_Amp_phi = datA[0, :, 1, ...] 
Y_Amp_theta = datA[0, :, 2, ...] 
Y_Amp_phi = datA[0, :, 3, ...] 
#
#Extracting the Phase
#
X_Phase_theta = datPh[0, :, 0, ...] 
X_Phase_phi = datPh[0, :, 1, ...] 
Y_Phase_theta = datPh[0, :, 2, ...] 
Y_Phase_phi = datPh[0, :, 3, ...] 
#
#header
#
h = hduA[0].header 
#
#Complex Values
#
Re_X_theta = X_Amp_theta*np.cos(X_Phase_theta)
Im_X_theta = X_Amp_theta*np.sin(X_Phase_theta)
Re_X_phi = X_Amp_phi*np.cos(X_Phase_phi) 
Im_X_phi = X_Amp_phi*np.sin(X_Phase_phi) 
Re_Y_theta = Y_Amp_theta*np.cos(Y_Phase_theta)
Im_Y_theta = Y_Amp_theta*np.sin(Y_Phase_theta)
Re_Y_phi = Y_Amp_phi*np.cos(Y_Phase_phi)
Im_Y_phi = Y_Amp_phi*np.sin(Y_Phase_phi)
#
#then
#
Xc_theta = Re_X_theta + 1j*Im_X_theta 
Xc_phi = Re_X_phi + 1j*Im_X_phi 
Yc_theta = Re_Y_theta + 1j*Im_Y_theta 
Yc_phi = Re_Y_phi + 1j*Im_Y_phi 
#
#compute the complex beams
#
print 'Compute Complex Beams'
#
xx = Xc_theta*np.conjugate(Xc_theta) +  Xc_phi*np.conjugate(Xc_phi) 
xy = Xc_theta*np.conjugate(Yc_theta) +  Xc_phi*np.conjugate(Yc_phi) 
yx = Yc_theta*np.conjugate(Xc_theta) +  Yc_phi*np.conjugate(Xc_phi) 
yy = Yc_theta*np.conjugate(Yc_theta) +  Yc_phi*np.conjugate(Yc_phi) 
#
'''xx = Xc_theta
xy = Xc_phi
yx = Yc_theta
yy = Yc_phi'''
#
#kat7__Jones_YY_R.fits
pyfits.writeto('kat7__Jones_XX_R.fits',xx.real,h[0:38], clobber=True)
pyfits.writeto('kat7__Jones_XX_Im.fits',xx.imag,h[0:38], clobber=True)  
pyfits.writeto('kat7__Jones_XY_R.fits',xy.real,h[0:38], clobber=True)
pyfits.writeto('kat7__Jones_XY_Im.fits',xy.imag,h[0:38], clobber=True)
pyfits.writeto('kat7__Jones_YX_R.fits',yx.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Jones_YX_Im.fits',yx.imag,h[0:38], clobber=True)   
pyfits.writeto('kat7__Jones_YY_R.fits',yy.real,h[0:38], clobber=True)
pyfits.writeto('kat7__Jones_YY_Im.fits',yy.imag,h[0:38], clobber=True)  
#Mueller_Matrices Conversion
#
print 'Mueller_Matrices Conversion'
m_ii = (xx*np.conjugate(xx) + xy*np.conjugate(xy) + yx*np.conjugate(yx)+ yy*np.conjugate(yy))/2. 
m_iq = (xx*np.conjugate(xx) - xy*np.conjugate(xy) + yx*np.conjugate(yx) - yy*np.conjugate(yy))/2. 
m_iu = (xy*np.conjugate(xx) + yy*np.conjugate(yx) + xx*np.conjugate(xy)+ yx*np.conjugate(yy))/2. 
m_iv = 1j*(xy*np.conjugate(xx) + yy*np.conjugate(yx) - xx*np.conjugate(xy) - yx*np.conjugate(yy))/2. 
m_qi = (xx*np.conjugate(xx) + xy*np.conjugate(xy) - yx*np.conjugate(yx) - yy*np.conjugate(yy))/2. 
m_qq = (xx*np.conjugate(xx) - xy*np.conjugate(xy) - yx*np.conjugate(yx)+ yy*np.conjugate(yy))/2. 
m_qu = (xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(xx) - yy*np.conjugate(yx))/2. 
m_qv = 1j*(xy*np.conjugate(xx) + yx*np.conjugate(yy) - yy*np.conjugate(yx) - xx*np.conjugate(xy))/2. 
m_ui = (xx*np.conjugate(yy) + yx*np.conjugate(xx) + xy*np.conjugate(yy) + yy*np.conjugate(xy))/2. 
m_uq = (xx*np.conjugate(yx) + yx*np.conjugate(xx) - xy*np.conjugate(yy) - yy*np.conjugate(xy))/2. 
m_uu = (xx*np.conjugate(yy) + xy*np.conjugate(yx) + yx*np.conjugate(xy)+ yy*np.conjugate(xx))/2. 
m_uv = 1j*(-xx*np.conjugate(yx) + xy*np.conjugate(yy) - yx*np.conjugate(xx) + yy*np.conjugate(xy))/2. 
m_vi = 1j*(xx*np.conjugate(yx) + xy*np.conjugate(yy) - yx*np.conjugate(xx) - yy*np.conjugate(xy))/2. 
m_vq = 1j*(xx*np.conjugate(yx) - xy*np.conjugate(yy) - yx*np.conjugate(xx) + yy*np.conjugate(xy))/2. 
m_vu = 1j*(xx*np.conjugate(yy) + xy*np.conjugate(yx) - yx*np.conjugate(xy) - yy*np.conjugate(xx))/2. 
m_vv = (xx*np.conjugate(yy) - xy*np.conjugate(yx) - yx*np.conjugate(xy) + yy*np.conjugate(xx))/2. 
#
# pyfits
#
print 'Mueller_Matrices into FITS'
#
pyfits.writeto('kat7__Mueller_GP_II.fits',m_ii.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_IQ.fits',m_iq.real,h[0:38], clobber=True)
pyfits.writeto('kat7__Mueller_GP_IU.fits',m_iu.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_IV.fits',m_iv.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_QI.fits',m_qi.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_QQ.fits',m_qq.real,h[0:38], clobber=True)
pyfits.writeto('kat7__Mueller_GP_QU.fits',m_qu.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_QV.fits',m_qv.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_UI.fits',m_ui.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_UQ.fits',m_uq.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_UU.fits',m_uu.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_UV.fits',m_uv.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_VI.fits',m_vi.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_VQ.fits',m_vq.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_VU.fits',m_vu.real,h[0:38], clobber=True) 
pyfits.writeto('kat7__Mueller_GP_VV.fits',m_vv.real,h[0:38], clobber=True) 

#
print 'DONE !!!!!!!!!!'
#
