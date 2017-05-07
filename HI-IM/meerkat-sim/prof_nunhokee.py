#
import numpy as np
import sys, os

from astropy.io import fits

from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import cygrid
import healpy as hp
import warnings
warnings.filterwarnings("ignore")
#
import argparse

def galactic_to_equicoord(lon, lat):
	"""
	Convert Galactic to Equatorial coordinates (J2000.0)
	(use at own risk)

	Input: [l,b] in decimal degrees
	Returns: [ra,dec] in decimal degrees   

	"""
	l = np.radians(lon)
	b = np.radians(lat)

	# North galactic pole (J2000) -- according to Wikipedia
	pole_ra = np.radians(192.859508)
	pole_dec = np.radians(27.128336)
	posangle = np.radians(122.932-90.0)    

	ra = np.arctan2( (np.cos(b)*np.cos(l-posangle)), (np.sin(b)*np.cos(pole_dec) - np.cos(b)*np.sin(pole_dec)*np.sin(l-posangle)) ) + pole_ra
	dec = np.arcsin( np.cos(b)*np.cos(pole_dec)*np.sin(l-posangle) + np.sin(b)*np.sin(pole_dec) )

	return np.rad2deg(ra), np.rad2deg(dec)
    
    #

#
def run_test():
    
	parser = argparse.ArgumentParser(description='Grid Specific position in Healpix into 2D FITS')

	parser.add_argument('-m', '--hlpx', type=str, dest='healpix', required=True, help='healpix filename')
	parser.add_argument('-x', '--xsize', type=int, dest='xaxis', default=256, help='x-axis for 2D FITS, Optional')
	parser.add_argument('-y', '--ysize', type=int, dest='yaxis', default=256, help='y-axis for 2D FITS, Optional')
	parser.add_argument('-p', '--psize', type=float, dest='pixsize', default=3.5/60., help='pixsize for 2D FITS, Optional')
	parser.add_argument('-c', '--dtype', type=int, dest='datatype', default=32, help='BITPIX for 2D FITS, Optional')
	parser.add_argument('-r', '--ra', type=float, dest='ra', default=0., help='Right Ascension, Optional')
	parser.add_argument('-d', '--dec', type=float, dest='dec', default=-30., help='Declination, Optional')
	parser.add_argument('-k', '--fwhm', type=float, dest='fwhm', default=60.0, help='image_fov in arcminutes, Optional')
	parser.add_argument('-e', '--ctype1', type=str, dest='ctype1', default='RA--SIN', help='Projection for CTYPE1, Optional')
	parser.add_argument('-f', '--ctype2', type=str, dest='ctype2', default='DEC--SIN', help='Projection for CTYPE2, Optional')

	args = parser.parse_args()    
	    
	    
	wmap_map_I = hp.read_map('%s' %args.healpix)

	
	header = fits.Header()

	header['SIMPLE'] = 'T'
	header['BITPIX'] = -args.datatype
	header['NAXIS'] = 3

	header['NAXIS1'] = args.xaxis
	header['NAXIS2'] = args.yaxis
	header['NAXIS3'] = 1

	header['CDELT1'] = -args.pixsize
	header['CDELT2'] = args.pixsize
	header['CDELT3'] = 1.

	header['CRPIX1'] = args.xaxis/2.
	header['CRPIX2'] = args.yaxis/2.
	header['CRPIX3'] = 0

	header['CRVAL1'] = args.ra
	header['CRVAL2'] = args.dec
	header['CRVAL3'] = 0.
	header['LATPOLE'] = 90.

	header['CTYPE1'] = 'GLON'  #args.ctype1  #
	header['CTYPE2'] =  'GLAT' # args.ctype2
	header['CTYPE3'] = 'VRAD'

	wcs = WCS(header)

	print '\n uploading header must in wcs representation into cygrid'
	gridder = cygrid.WcsGrid(header)


	NSIDE = hp.get_nside(wmap_map_I)

	# get_map_size

	NPIX = hp.nside2npix(NSIDE) 			

	#   phi - longitude in radians
	#   theta - latitude in radians

	theta, phi = hp.pix2ang(NSIDE, np.arange(NPIX))


	# convert phi, theta to degrees

	lons = np.rad2deg(phi).astype(np.float64)
	lats = (90. - np.rad2deg(theta)).astype(np.float64)
	

	kernelsize_fwhm = args.fwhm/60. 			#	 converts arcminutes to degrees
	kernelsize_sigma = kernelsize_fwhm / 2.355 		#     https://en.wikipedia.org/wiki/Full_width_at_half_maximum
	sphere_radius = 3. * kernelsize_sigma

	gridder.set_kernel(
	    'gauss1d',
	    (kernelsize_sigma,),
	    sphere_radius,
	    kernelsize_sigma / 2.
	    )
	#
	if (args.ctype1 != 'GLON') & (args.ctype1 != 'GLAT'):
		ra, dec = galactic_to_equicoord(lon = lons, lat = lats)
		gridder.grid(ra, dec, wmap_map_I[:, None])
		
	gridder.grid(lons, lats, wmap_map_I[:, None])

	datacube = gridder.get_datacube().squeeze()


#	fig = plt.figure()


#	cel_header = wcs.celestial.to_header()
#        import pywcsgrid2
#	
#	pywcsgrid2.subplot(111, header=cel_header)
#	plt.imshow(datacube, origin="lower")
#	
#	plt.show()

	fits.writeto('fullsky_regridded.fits', datacube, cel_header, clobber=True)
	
if __name__ == '__main__':
	run_test()

