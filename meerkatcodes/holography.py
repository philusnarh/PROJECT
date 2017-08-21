''' Created by Khan M. B. Asad (khmbasad@gmail.com), 17 May 2017 '''

import numpy as np, time
import katholog
import astropy.io.fits as fits
import datetime

def load_data(input_file, diameter=6, freqs=range(1400,1421), ant=0, Npix=512):
    """
    Read MeerKAT holographic data and return katholog beamcube as npy
    Input: .h5 file, diameter [deg], freqs [list]
    """
    # download data from kat archive
    #os.system('wget http://kat-archive.kat.ac.za/archive/data/RTS/telescope_products/2017/02/23/1487813282.h5')
    start = time.time()

    # Read the dataset
    dataset = katholog.Dataset(input_file,'meerkat',method='direct',dobandpass=True,onaxissampling=0.1)
    
    # The variables
    ants = dataset.radialscan_allantenna
    scanners = dataset.scanantennas
    trackers = dataset.trackantennas
    extent = diameter
    clipextent = diameter
    bandwidth = freqs[1] - freqs[0]
    nchan = len(freqs)

    scanner = ants[scanners[ant]]
    tracker = ants[trackers[0]]

    # Flag hours
    flags_hrs=dataset.findworstscanflags(freqMHz=freqs,dMHz=bandwidth, scanantennaname=scanner, trackantennaname=tracker, doplot=True)
    dataset.flagdata(flags_hrs=flags_hrs)

    # Create the beamcube
    b = katholog.BeamCube(dataset,freqMHz=freqs,dMHz=bandwidth, scanantennaname=scanner,\
            interpmethod='scipy', applypointing=None,extent=extent, gridsize=Npix)
    
    # Write as Jones matrices
    beams = np.zeros((nchan,2,2,Npix,Npix), dtype='c16')
    filename = "%s_%s_%ipx_%iMHz_%iMHz_%ichannels_Jones"%(input_file[:-3],ants[scanners[ant]],Npix,freqs[0],bandwidth,nchan)
    # open text file to write beam widths and offsets of all channels
    t = open(filename+'.txt', 'w')
    t.write("widthGxH, widthGxV, widthGyH, widthGyV, offsetGxH, offsetGxV, offsetGyH, offsetGyV\n")
    for i,f in enumerate(freqs):
        print "-> %i"%i
        beams[i,0,0,:,:] = b.Gx[i,...]
        beams[i,0,1,:,:] = b.Dx[i,...]
        beams[i,1,0,:,:] = b.Dy[i,...]
        beams[i,1,1,:,:] = b.Gy[i,...]
        t.write("%f, %f, %f, %f, %f, %f, %f, %f\n"%\
                          (b.beamwidthGx[i][0],b.beamwidthGx[i][1],b.beamwidthGy[i][0],b.beamwidthGy[i][1],\
                           b.beamoffsetGx[i][0],b.beamoffsetGx[i][1],b.beamoffsetGy[i][0],b.beamoffsetGy[i][1]) )
    t.close()

    # Save as fits
    write_fits(dataset, beams, filename, freqs, diameter, scanner, tracker)

    end = time.time()
    print "-> Time taken: %.2f minutes"%((end-start)/60.)
    return beams

def write_fits(dataset, beamcube, filename, freqs, diameter, scanner, tracker):
    # Calculate parameters from the dataset
    el = dataset.env_el[0]
    az = np.mean(dataset.scanaz*(180./np.pi))
    ID = dataset.filename.split('/')[-1][:-3]
    
    # Create header
    hdr = fits.Header()
    hdr['ID'] = ID
    ctypes = ['AZIMUTH', 'ELEVATION', 'JONES0', 'JONES1', 'FREQ']
    crvals = [az, el, 0,0, freqs[0]/1e3]
    cdelts = [diameter/beamcube.shape[-2], diameter/beamcube.shape[-1], 1,1, freqs[1]/1e3-freqs[0]/1e3]
    cunits = ['deg', 'deg', '', '', 'GHz']
    for i in range(len(beamcube.shape)):
        ii = str(i+1)
        hdr['CTYPE'+ii] = ctypes[i]
        hdr['CRVAL'+ii] = crvals[i]
        hdr['CDELT'+ii] = cdelts[i]
        hdr['CUNIT'+ii] = cunits[i]
    hdr['TELESCOP'] = 'MeerKAT'
    hdr['OBSSTART'] = dataset.rawtime[0]
    hdr['OBSEND'] = dataset.rawtime[-1]
    hdr['DURATION'] = (dataset.rawtime[-1]-dataset.rawtime[0])/3600.
    hdr['SCANNER'] = scanner
    hdr['TRACKER'] = tracker
    hdr['TARGET'] = dataset.target.name
    hdr['TARGRA'] = dataset.target.radec()[0]*(180./np.pi)
    hdr['TARGDEC'] = dataset.target.radec()[1]*(180./np.pi)
    hdr['DATE'] = str(datetime.datetime.now())
    
    # Write real and imag parts of data
    hdr['PART'] = 'REAL'
    hdu_r = fits.PrimaryHDU(beamcube.real, header=hdr)
    hdr['PART'] = 'IMAG'
    hdu_i = fits.PrimaryHDU(beamcube.imag, header=hdr)
    hdu_r.writeto(filename+'_real.fits', clobber=True)
    hdu_i.writeto(filename+'_imag.fits', clobber=True)

if __name__=='__main__':
    path = "/home/asad/data/meerkat/beam/holography/"
    load_data(path+"1487813282.h5", freqs=range(856,856+857), ant=2, Npix=128)

"""
d = katholog.Dataset('1487813282.h5', 'meerkat')
b = katholog.BeamCube(d,freqMHz=range(900,1421,1),dMHz=1,scanantennaname='m017',interpmethod='scipy',applypointing='Gx',extent=6)
One antenna takes 25.72 minutes for (521, 2, 2, 256, 256) shape
"""
