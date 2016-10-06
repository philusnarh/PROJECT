from scipy import signal
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
#
hdulist = pf.open('mueller_4_x_4_beam.fits')
beam_data = hdulist[0].data
header = hdulist[0].header
crpix2, crval2, cdelt2 = [ header.get(x) for x in "CRPIX2", "CRVAL2", "CDELT2" ]
bm_data = beam_data[0]
#print plt.imshow(bm_data)
#	masked_beam= beam_data[iter]
shape = bm_data.shape
r = 1.0
rad = np.linspace(-shape[0]/2,shape[-1]/2,shape[0])
rad2d =  np.sqrt(rad[np.newaxis,:]**2+rad[:,np.newaxis]**2)
#mask = rad2d <= radius/abs(cdelt2)
mask = rad2d <= r/abs(cdelt2)
bm_2_bx = np.ones(shape = (shape[0],shape[-1]), dtype=np.float64)
#masked_beam = bm_data*mask
masked_beam = bm_2_bx*mask
# write to fits
pf.writeto('boxcar_beam.fits', masked_beam, header = header, clobber=True)
#pf.writeto('beam.fits', masked_beam, header = header, clobber=True)  
#print masked_beam
plt.imshow(masked_beam)
plt.colorbar()
plt.axis('off')
#masked_beam = np.where(masked_beam > 0,  1, masked_beam)

#x = np.linspace(-2.5, 2.5, 256)
plt.figure(2)
freq = rad*(5./shape[0])
radius = 0.25
while radius <= 2.50:
	mask = rad2d <= radius/cdelt2
#print mask
	masked_beam0 = masked_beam*mask
	#masked_beam = np.ma.masked_where(masked_beam == 0., masked_beam)
#u = response.copy()
	plt.plot(freq, masked_beam0.diagonal(), label = 'radius $%s^\circ$'%radius)
	radius = radius + 0.25
#plt.plot(freq, response)
#plt.axis([-0.5, 0.5, -120, 0])
plt.title("Angular response of the boxcar window" ,fontsize="15")
plt.ylabel("Normalized magnitude [dB]" ,fontsize="15")
plt.xlabel("Angular coordinate / deg" ,fontsize="15")
plt.legend(loc = 'best') 

plt.show()
