import numpy as np
import healpy as hp
import fmt, ephem
from math import *
import pyfits as pf
 
def ga2equ(ga):
    """
    Convert Galactic to Equatorial coordinates (J2000.0)
    The coordinates should be given in degrees
    """

    l = radians(ga[0])
    b = radians(ga[1])
    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)

    ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
    return np.array([degrees(ra), degrees(dec)])


def interBeam(xpol,freq,jd,ra,dec):
  """
  Interpolate beam at given frequency and sky position (RA and DEC)
  """

  # Reading FEKO BEAMS ..."
  fekoX=fmt.FEKO(xpol)
  fields = fekoX.fields
  feko_xpol=fekoX.fields[0]
  phi=feko_xpol.phi*np.pi/180. # azimuth
  theta=feko_xpol.theta*np.pi/180. # zenith
  theta = np.pi/2 - theta # pyephem wants the elevation rather than the zenith angle
  gxx=feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi)
  gxx = gxx/np.max(gxx) # normalizing to one
  ind1 = np.where(gxx<=1e-6) 
  ind2 = np.where(gxx>=1e-6)
  gxx[ind1] = np.nan
  #gxx[ind2] = 1

  # Create OBSERVER
  paper = ephem.Observer()
  paper.lat, paper.long, paper.elevation = '-30:43:17', '21:25:40.08', 0.0
  j0 = ephem.julian_date(0)
  paper.date = float(jd) - j0 + 5./60/24

  beamFreqs = np.zeros(11)
  beamFreqs[0] = 100e6
  for i in range(1,11):
     beamFreqs[i] =  beamFreqs[i-1] + 10e6
  
  beamRA = np.ndarray(shape=(phi.shape[0]),dtype=float)
  beamDEC = np.ndarray(shape=(phi.shape[0]),dtype=float)
  for k in range(beamRA.shape[0]):
    ra0,dec0 = paper.radec_of(phi[k],theta[k])
    beamRA[k] = ra0  # RA in radians
    beamDEC[k] = dec0    # DEC in radians
    
  print "FREQUENCY INTERPOLATION"
  InterBeamF = np.ndarray(shape=(gxx.shape[0]),dtype=complex)
  dist = np.abs(freq - np.array(beamFreqs))
  ind = np.argsort(dist)
  ind.flatten()
  print gxx.shape
  if dist[ind[0]] ==0:
     InterBeamF[:] = gxx[ind[0]]
  else:
     weights = dist[ind[0]]**(-2) + dist[ind[1]]**(-2)
     InterBeamF[:] = (gxx[ind[0]]*dist[ind[0]]**(-2) + gxx[ind[1]]*dist[ind[1]]**(-2))/weights

  InterHealBeam = np.ndarray(shape=(ra.shape[0]),dtype=complex)
  print "SPATIAL INTERPOLATION"

  for r in range(ra.shape[0]):
    ra[r] = ra[r] - 2*np.pi if ra[r]>(2*np.pi) else ra[r]
    dist = np.sqrt((ra[r] - beamRA)**2 + (dec[r] - beamDEC)**2)
    ind = np.argsort(dist)
    ind = ind.flatten()

    InterHealBeam[r] = gxx[ind[0]]

  np.savez(open('Beam-f%.4g_j%.5f.npz'%(freq*1e-6,np.float(jd)),'wb'),beam=InterHealBeam)
  return InterHealBeam
     
if '__name__==__main__':
   xpol = 'PAPER_FF_X.ffe'
   #nside = 64
   #interBeam(xpol,freq,jd,ra,dec)
   nside = 32
   npix = hp.nside2npix(nside)   
   freqs = np.linspace(100e6,200e6,100)
   #freq = 155e6
   #si = -2
   si=pf.open('The_spectral_index_using_landecker_and_haslam1.fits')
   si=si[0].data
   si= hp.ud_grade(si,nside)

   
   rad2deg = lambda val: val * 180./np.pi
   deg2rad = lambda val: val * np.pi/180
   coord = hp.pix2ang(nside,np.arange(npix))
   theta = 90-(rad2deg(coord[0])) # galactic latitute
   phi = rad2deg(coord[1])   # galactic longitude
   ra = np.ndarray(shape=(npix),dtype=float)
   dec = np.ndarray(shape=(npix),dtype=float)

    # converting galactic coordinate to equatorial coordinates
   for ii in range(len(theta)):
       eq_coord = ga2equ([phi[ii],theta[ii]])
       ra[ii] = deg2rad(eq_coord[0])
       dec[ii] = deg2rad(eq_coord[1])
   for freq in freqs:
      _data = np.loadtxt('haslam.dat')
      _data = hp.ud_grade(_data, nside)
      _data_f150 = _data*(freq/408e6)**si
      np.save('mmap_%s' % (freq), _data_f150)

      jd = 2455819.46147 #in LST= 6 Hours,48mins 41.20 seconds

   #jds = np.load('julian_dates.npz')['jd'] 
   #for jd in jds:
   #for freq in freqs:
      #_data_f150 = _data*(i/408e6)**si
      #for jd in jds:
      jones = interBeam(xpol,freq,jd,ra,dec)
      #jones = interBeam(xpol,i,jd,ra,dec)
      out_data = _data_f150 * jones  # integrating the  map with the beam
      outfile = 'map_f%.4g_j%.5f.npz'%(freq*1e-6,np.float(jd))
      #outfile = 'map_f%.4g_j%.5f.npz'%(i*1e-6,np.float(jd))
      np.savez(open(outfile,'wb'),beam=jones,data=out_data)


#print "out_data",out_data

#integral = sum(out_data)
#print integral

#primary_hdu=pf.PrimaryHDU(out_data) # generate a primary header list input the beta array
#main_hdu=pf.HDUList([primary_hdu]) # create the main header lists and put in the primary header 
#main_hdu.writeto("beta.fits") # write the entire contents of the main hdulist into the new fits file,thus creating a new fits file




