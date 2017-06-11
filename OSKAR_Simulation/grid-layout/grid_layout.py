#!/usr/bin/env python
import pylab
import numpy
import sys
from argparse import ArgumentParser


def grid_layout_enu(ants_x, ants_y, step_x, 
                    step_y, lonlat=None, filename=None, 
                    dish_diameter=13.5, mount="ALT-AZ", 
                    savefig=None):
    """
    ants_x : Antennas along x-axis (East)
    ants_y : Antennas along y-axis (North)
    step_x : Step size along x-axis (metres)
    step_y : Step size along y-axis (metres)
    lonlat : Telescope location (degrees)
    """
    
    nants = ants_x * ants_y
    enu = []
    
    if lonlat in [None, [], ""]:
        lonlat = 0,0
    
    zero_x = ants_x/2 - 0.5 if ants_x%2==0 else ants_x/2
    zero_y = ants_y/2 - 0.5 if ants_y%2==0 else ants_y/2
    
    for i in xrange(ants_x):
        for j in xrange(ants_y):
            enu.append( ( (i-zero_x)*step_x + lonlat[0], 
                          (j-zero_y)*step_y + lonlat[1]))
    
    if filename:
        with open(filename, "w") as wstd:
            wstd.write("# East North Up Dish_Diameter Station Mount\n")
            for i,(x,y) in enumerate(enu):
                
                wstd.write("%.6f %.6f 0.0 %.3f ST-%d %s\n"%(x, y, 
                           dish_diameter, i, mount))
                
    enu = numpy.array(enu).T                

    if savefig:
        pylab.figure(figsize=(12,12))
        pylab.scatter(enu[0], enu[1])
        pylab.xlabel("East [m]")
        pylab.ylabel("North [m]")
        pylab.grid()
        pylab.savefig(savefig)
        
    return enu


if __name__=="__main__":

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
 
    parser = ArgumentParser(description="Create an antenna layout in ENU coordinates")

    add = parser.add_argument

    add("prefix", nargs="?",
            help="Prefix for output products")

    add("-nx", "--ants-x", type=int, default=8,
            help="Number of antennas along the x-axis")

    add("-ny", "--ants-y", type=int, default=1,
            help="Number of antennas along y-axis")

    add("-dx", "--delta-x", type=float, default=100,
            help="Step size along x-axis")

    add("-dy", "--delta-y", type=float, default=100,
            help="Step size along y-axis")

    add("-dd", "--dish-diameter", type=float, default=2,
            help="Dish diameter [metres]")

    add("-ll", "--lon-lat",
            help="Telescope location. Comma separated longitude and lattitude values in degrees")

    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = "grid_layout_ENU-%dx%d_%dx%d"%(args.ants_x, args.ants_y,
                                                    args.delta_x, args.delta_y)

    if args.lon_lat:
        lonlat = map(float, args.lon_lat.split(","))
    else:
        lonlat = 0,0

    grid_layout_enu(args.ants_x, args.ants_y, 
                    args.delta_x, args.delta_y,
                    lonlat=lonlat,
                    dish_diameter=args.dish_diameter,
                    filename=prefix+".txt",
                    savefig=prefix+".png")

