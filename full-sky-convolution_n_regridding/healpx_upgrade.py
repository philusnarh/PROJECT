#!/usr/bin/env python
"""
	healpx.py
====================

	* This program write healpix map into healpix file. 
	* It also perfomrs regridding


"""
#

import healpy as hp
import pyfits
from argparse import ArgumentParser
import sys
import numpy as np
#

def openFitsFile(filename):
	    """ Opens a FITS Filen using PyFITS
	    Parameters
	    ----------
	    filename:str
		path to file
	    """
	    try:
		hdu = pyfits.open(filename)
		return hdu
	    except:
		print "Error: cannot open %s"%filename
		raise
		
#
def run_sim_test():
	#	
	
	for i, arg in enumerate(sys.argv):		
        	if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
 
	parser = ArgumentParser(description="Writes a healpix map into a healpix file.")

	add = parser.add_argument

	add("healpixname",
	    help="HEALpix Map")

	add("-n", "--chan", type=int,  default=0,
	    help="Channel Number")
	
	add("-s", "--stokes", type=str,  default="IQU", required=True,
	    help="Choose Stokes I or IQU")
	    
	add("-f", "--filename", type=str, required=True, 
	    help="HEALpix filename")
	    
	add("-u", "--upgrade", type=int, default=0,
	    help="Choose 1 to upgrade OR 0 otherwise")
	    
    	add("-o", "--newnside", type=int, default=128,
	    help="Choose new NSIDE to upgrade")

	args = parser.parse_args()	

	hdu = openFitsFile('%s' %args.healpixname)
	hp_data = hdu[0].data[args.chan]	
	
	
	if (args.stokes == 'I') | (args.stokes == 'i'):
		mp = hp_data[0,...]
		
		
		if args.upgrade != 0:
			mp = hp.pixelfunc.ud_grade(mp, args.newnside)
			print np.array(mp).shape
		nside = hp.pixelfunc.get_nside(m=mp)			
		hp.fitsfunc.write_map(filename ='%s_%s_NSIDE_%d.fits' %(args.filename, args.stokes, nside),m = mp)
		
	if (args.stokes == 'IQU') | (args.stokes =='iqu'):
		mp = hp_data[0:3,...]
		
		if args.upgrade != 0:
			mp = hp.pixelfunc.ud_grade(mp, args.newnside)
			print np.array(mp).shape
		nside = hp.pixelfunc.get_nside(m=mp)			
		hp.fitsfunc.write_map(filename ='%s_%s_NSIDE_%d.fits' %(args.filename, args.stokes, nside),m = mp)
		
	
if __name__ == '__main__':
	run_sim_test()
	
