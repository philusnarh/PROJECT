import numpy as np
import pyfits as pf
from matplotlib import pylab as plt
from scipy.fftpack import fftfreq
#
#
def openFitsFile(filename):
	""" Opens a FITS Files using PyFITS
        Parameters
        ----------
        filename:str
        path to file  """
        try:
        	hdu = pf.open(filename)        	
        	return hdu
        except:
        	print "Error: cannot open %s"%filename
		raise
#
def shift_fft(input_array,shift):	
	    shift_rows,shift_cols = shift
	    nr,nc = input_array.shape
	    Nr, Nc = fftfreq(nr), fftfreq(nc)
	    Nc,Nr = np.meshgrid(Nc,Nr)
	    fft_inputarray = np.fft.fft2(input_array)
	    fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
	    output_array = np.fft.ifft2(fft_inputarray*fourier_shift)
   	    return np.real(output_array)
#
def run_shiftfun():
#
	hdu_gridd = openFitsFile('./oskarr_beam_XY__dis_mueller_4_x_4_beam.fits')
	datt = hdu_gridd[0].data[0]
	print plt.figure(1);	plt.imshow(datt); plt.colorbar()
	#
	newdat = shift_fft(input_array = datt, shift = [2,1])
	#header = hduAmp[0].header 
	print datt == newdat
	print plt.figure(2);	plt.imshow(newdat); plt.colorbar()
	print plt.figure(3); plt.plot(datt.diagonal()); plt.plot(newdat.diagonal())
	print datt - newdat
	plt.grid()
	plt.show()
	#
#
if __name__ == '__main__':
	run_shiftfun()

