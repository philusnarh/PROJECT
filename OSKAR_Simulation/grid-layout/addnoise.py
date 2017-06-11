#!/usr/bin/env python

from pyrap.tables import table
import numpy
import sys
from argparse import ArgumentParser


def vis_noise(msname, sefd):

    tab = table(msname)
    spwtab = table("%s/SPECTRAL_WINDOW"%msname)

    freq0 = spwtab.getcol("CHAN_FREQ")[0, 0]
    wavelength = 300e+6/freq0
    bw = spwtab.getcol("CHAN_WIDTH")[0, 0]
    dt = tab.getcol("EXPOSURE", 0, 1)[0]
    dtf = (tab.getcol("TIME", tab.nrows()-1, 1)-tab.getcol("TIME", 0, 1))[0]

    tab.close()
    spwtab.close()

    noise = sefd/numpy.sqrt(abs(2*bw*dt))

    return noise


def addnoise(msname, sigma=None, sefd=None, column="DATA"):
    
    tab = table(msname, readonly=False)

    if sefd:
        noise = vis_noise(msname, sefd)

    nrows = tab.nrows()
    data = tab.getcol(column)
    dshape = list(data.shape)

    rowchunk = nrows/10 if nrows > 100 else nrows

    for row0 in range(0, nrows, rowchunk):
        nr = min(rowchunk, nrows-row0)
        dshape[0] = nr

        data = noise*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape))

        print("Adding noise to %s (rows %d to %d)"%(column ,row0, row0+nr-1))

        tab.putcol(column, data, row0, nr) 
    tab.close()
    
    print("Done!")

if __name__=="__main__":
    
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
 
    parser = ArgumentParser(description="Add noise to an MS")

    add = parser.add_argument

    add("msname",
            help="MS Name")

    add("-n", "--noise", type=float,
            help="Noise per visibility")

    add("-s", "--sefd", type=float,
            help="System Equivalent Flux Density")

    add("-c", "--column", default="DATA",
            help="Column")

    args = parser.parse_args()

    if args.sefd:
        addnoise(args.msname, sefd=args.sefd, column=args.column)
        
    else:
        addnoise(args.msname, noise=args.noise, column=args.column)

