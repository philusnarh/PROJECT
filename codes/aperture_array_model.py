#!/usr/bin/env python
"""
dipol_generation.py
================

This program generates random dipoles corresponding to Aperture Illumination to model the beampattern of an antenna.

Requirements
-------------
ds9
numpy, scipy, matplotlib

"""

#
from matplotlib import pylab as plt
import numpy as np
import scipy.interpolate as interpolate
import os
import optparse,sys
#
def inverse_transform_sampling(hist, bin_edges):
    cum_values = np.cumsum(hist[1:]*np.diff(bin_edges)) 	# Cumulative Distribution Function for each bin
    #inv_cdf1 = interpolate.InterpolatedUnivariateSpline(cum_values, bin_edges[1:])# Inverse Transform using scipy interpolate function
    inv_cdf1 = interpolate.UnivariateSpline(cum_values, bin_edges[1:])
    pdf = np.arange(0,1, 0.0001)
    #pdf = np.random.uniform(0,1, 1000)
    return inv_cdf1(pdf)
#
def run_test():
    #
    #	CREATING ANTENNA POSITIONS & OSKAR CONFIGURATION SETUP
    #
    print " ------------------------------------------------------"
    print "\n >> Enter the Telescope Directory: "
    print " ------------------------------------------------------"
    #telescope_dir =raw_input()
    #print telescope_dir
    #
    print " ------------------------------------------------------"
    print "\n >> Opening Text File to Enter Antenna Position: "
    print " ------------------------------------------------------"
    #os.system(" vi telescope_dir/layout.txt")
    #
    print " ------------------------------------------------------"
    print "\n >> Opening Oskar Configuration Setup: "
    print " ------------------------------------------------------"
    #os.system(" vi telescope_dir/config_setup.ini")
    #
    #	CREATING  DIRECTORIES FOR ANTENNA STATIONS 
    #
    print " ------------------------------------------------------"
    print "\n >> Creating Stations to Store Dipoles: "
    print " ------------------------------------------------------"
    #posnSiz = input('Number of Aperture array: ')
    #tel_dir = sys.argv[1]
    #posnSiz = int(sys.argv[2])
    #p2 = sys.argv[3]
    #p3 = sys.argv[4]
    os.system('mkdir $(seq -f "$HOME/telescopemodel/station%03g" 7)')
    #os.system('mkdir $(seq -fr "$Home/telescope_dir/station%03g" 7)')
    #for k in range (posnSiz):
    #	return os.system("mkdir telescope_dir/station%d%d"%(0,k))
    #	os.system("for dir in *; do [ -d "$dir" ] && cp -vfr hh.txt "$dir"; done ")
    #
    #	CREATING  RADIAL DISTRIBUTION OF DIRECTORIES FOR ANTENNA STATIONS 
    #
    bin_edges = np.arange(0.0,1,0.01) #100 values: 0, 0.01, 0.02, etc.
    #bin_edges = np.random.uniform(0.0,1,500)
    #Determines the Distributions of Aperture Illumination 
    Type = input("Enter the value 0 for Gaussian OR 1 for Exponential ...: ")
    if Type == 0:
       Jpdf = (1/(2*np.pi))*np.exp(-(bin_edges **2)/2)  # Gaussian distribution of random deviates(U)
       plt.figure(1)
       rad = inverse_transform_sampling(Jpdf, bin_edges)
       plt.hist(rad,80,normed=True)
       plt.xlabel('Radius')
       plt.ylabel('Illumination Function')
       plt.title('Distribution of Aperture Illumination')
    elif Type == 1:
       Jrad = np.exp(-bin_edges) 	# Exponential distribution of random deviates(U)
       fig = plt.figure(1)
       rad = inverse_transform_sampling(Jrad, bin_edges)
       plt.hist(rad, 80, normed=True)
       plt.xlabel('Raduis')
       plt.xlim(0,3)
       plt.ylabel('Illumination Function')
       plt.title('Distribution of Aperture Illumination')
    else:
       print "ENTER THE VALUE 0 OR 1"
       raise
    #
    #	CREATING Circular Aperture Array Illumination layout
    #
    a = rad[np.where(rad<=1)]
    U = np.random.uniform(0, 2*np.pi, len(a))
    theta = np.random.uniform(0, np.pi, len(a))
    fig1 = plt.figure(2)
    x = 6*a*np.cos(U)		
    y = 6*a*np.sin(U)		
    plt.plot(x, y, 'b+') 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Station Setup')
    plt.show()
    # Save Dipoles
    np.savetxt('layout.txt' ,np.array([x,y, np.zeros(len(x))]).T, delimiter = ',') 
    #print "Dipoles are saved as ....: dipoless.txt" 
    #os.system("cp -vrf layout.txt /home/narh/sw/KAT7_ANTENNA1/telescopeS1/station???")
    #os.system("rm -vrf /home/narh/sw/KAT7_ANTENNA1/telescopeS1/station006/station00*")
    #os.system('cd')
    # OSKAR SIMULATION
    os.system('for dir in $HOME/telescopemodel/station*; do [ -d "$dir" ] && cp -vfr layout.txt "$dir"; done ')
    os.system('oskar_sim_beam_pattern $HOME/telescopemodel/oskar_congifuration_setup.ini')
    os.system('ds9 -tile beam_pattern_VOLTAGE.fits -cmap Rainbow -zscale -zoom to fit -cube play')
    #os.system('tigger beam_pattern_VOLTAGE.fits')
if __name__ == '__main__':
    run_test()

