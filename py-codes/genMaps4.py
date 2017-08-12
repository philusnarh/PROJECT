# This script generates interpolated 
# beams at 100 to 200 MHz in steps of 1 MHz frequencies 

import numpy as np
import healpy as hp
import fmt, ephem
from math import *
import pyfits as pf
import datetime, time
import pylab as pl


def interplotate1DMissingFreqsbeams(data,dataStartFreq,dataEndFreq,newbeamsFreqStep):
    
    nchan, xpl = np.array(data).shape             # Get shape of the data

    #sfreq = np.array([dataStartFreq + i*dataStepFreq for i in range(nchan)])  # Calc. source data freqs.
    sfreq = np.linspace(dataStartFreq, dataEndFreq, nchan)  # Calc. source data freqs.  
	
    new_chan_numb = (dataEndFreq - dataStartFreq)/newbeamsFreqStep
    nfreq = np.arange(sfreq[0],sfreq[-1],newbeamsFreqStep)  # Calc. new beams freqs. # np.linspace(900e+6,1300e+6,7)
    #print len(nfreq)
    #print sfreq[0]
    #print sfreq[1]
    #print newbeamsFreqStep
    #raise 
    new_beams = []                                    # variable for storing new beams
    
    for ifreq, freq in zip(range(len(nfreq)), nfreq):      # for each new beam ...
	#print 'ifreq:', ifreq

        index = {}                                 # variable for storing freq. diff. indexes
        
        freq_diff = abs(sfreq-freq) 		    # Calc. the frequency difference
        #print "freq_diff=",len(freq_diff)

        if 0 in freq_diff:                           # Checking if current frequency already exist in source beam then add it

            tindex = np.where(freq_diff == 0)[0][0]   # Get the index of the beams
            new_beams.append(data[tindex])
            
        else:
            
            it = 0
            for ivalue in freq_diff:                        # Store the indexes of the frequency differencies before sorting
                
                index[ivalue] = it
                it += 1

            freq_diff.sort()                                    # Sort difference in increasing order

            beam1_index = index[freq_diff[0]]                   # Get the index/channel number of 1st freq.
            #print  "beam1_index=", beam1_index 
            beam2_index = index[freq_diff[1]]                   # Get the index/channel number of 1st freq. 

            dv = abs(sfreq[beam2_index]- sfreq[beam1_index])

            weight1 = freq_diff[1]*1.0/dv                       # Compute the 1st freq. weight
            weight2 = freq_diff[0]*1.0/dv                       # Compute the 2nd freq. weight
            beam1 = data[beam1_index]
            beam2 = data[beam2_index]

            new_beams.append((beam1*weight1 + beam2*weight2))   # Compute the target beams using the freq. weight

    return np.array(new_beams)



 
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
  fekoX=fmt.FEKO(xpol)    #reading the "x polarization data i.e the antenna gain" stored in tnkhe FEKO format. eg <fmt.FEKO instance at 0x7f97136b8368>
  #fields = fekoX.fields   # there are about 11  different fields in  the FEKO format
  feko_xpol=fekoX.fields[0]  # out of the 11 the first field  is selected because that is where our data is.
  phi=feko_xpol.phi*np.pi/180. # the  azimuth angle in radians 
  theta=feko_xpol.theta*np.pi/180. #  the zenith angle in radians
  theta = np.pi/2 - theta # pyephem wants the elevation rather than the zenith angle
  gxx=feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi) 
  #the antenna gain in complex array is multiplied by it's complex conjugate to give the amplitude of the signal  
  #print "before_gxx=",gxx
  gxx = gxx/np.max(gxx) # normalizing to one
  #print "after_gxx=",gxx
  
  ind1 = np.where(gxx<=1e-6) 
  ind2 = np.where(gxx>=1e-6)
  gxx[ind1] = np.nan


  bm = np.zeros(shape=(11,gxx.shape[0]),dtype=complex)
  print "bm.shape=",bm.shape 
  for jj in xrange(11):
	
	  feko_xpol=fekoX.fields[jj]  # out of the 11 the first field  is selected because that is where our data is.
	  phi=feko_xpol.phi*np.pi/180. # the  azimuth angle in radians 
	  theta=feko_xpol.theta*np.pi/180. #  the zenith angle in radians
	  theta = np.pi/2 - theta # pyephem wants the elevation rather than the zenith angle
	  bm[jj,:]=feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi) 
          #the antenna gain in complex array is multiplied by it's complex conjugate to give the amplitude of the signal  
	  #print "before_gxx=",gxx
	  bm[jj,:] = bm[jj]/np.max(bm[jj]) # normalizing to one
	  #print "after_gxx=",gxx
	  

	  # masking the northern hemisphere with nans

	  ind1 = np.where(bm[jj]<=1e-6)  
          #print "ind1=", ind1
	  #ind2 = np.where(bm[jj]>=1e-6)
	  bm[jj,:][ind1] = np.nan        
  
  print bm.shape
  print bm[0] == bm[1]
  #hp.mollview(hp.ud_grade(bm[0].real, 128),cmap=pl.get_cmap('jet'))
  print 'Interbeam'
  InterBeamF = interplotate1DMissingFreqsbeams(data=np.abs(np.nan_to_num(bm)),
						dataStartFreq=float(freq[0]),
						dataEndFreq=float(freq[1]),
						newbeamsFreqStep=float(freq[2]))
  #pl.show()
  InterBeamF.shape
  #raise

  
  
  #gxx[ind2] = 1

  # Create OBSERVER

  """
   PyEphem provides an ephem Python package for performing high-precision astronomy computations.
   
   "ephem" is a short word for ephemeris, which is the traditional term for a table giving the position of a planet, asteroid, or comet for a series of dates. It can:
   
   i.compute where in the sky an object appears for a particular observer on Earth when given the latitude,longitude and elevation
   ii.Return the Julian Date corresponding to any calendar date
    
   """
  paper = ephem.Observer() # it creates an observation , it tells the current date and time, position etc eg.date='2017/4/16 19:12:12' epoch='2000/1/1 12:00:00' lon='0:00:00.0' lat='0:00:00.0' elevation=0.0m horizon=0:00:00.0 temp=15.0C pressure=1010.0mBar nb it does not take any argument
  paper.lat, paper.long, paper.elevation = '-30:43:17', '21:25:40.08', 0.0  # giving ephem the desired latitude,longitude and elevation 
  j0 = ephem.julian_date(0) # gives the julian date of the present day
  paper.date = float(jd) - j0 + 60./60/24    #   derive the normal calendar date  from the julian date

  '''
	 #   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 """
	  generate the beam frequency from 100 to 200 MHz in steps of 10MHz
	  
	  """
	  beamFreqs = np.zeros(11)   # create an array of eleven zeros 
	  beamFreqs[0] = 100e6       # replace the first element in the arrays of zeros with "100000000" i.e 100e6 
	  for i in range(1,11):
	     beamFreqs[i] =  beamFreqs[i-1] + 10e6  # increase each element in the array above in steps of 10e6
	  
	  beamRA = np.ndarray(shape=(phi.shape[0]),dtype=float)   # generate a random array of equal size or length as in phi eg. phi has a length of 6643 hence we generate an array of 6643 elements  
	  beamDEC = np.ndarray(shape=(phi.shape[0]),dtype=float)  
	  for k in range(beamRA.shape[0]):
	    ra0,dec0 = paper.radec_of(phi[k],theta[k])  # convert from Alt Al system to ra and dec 
	    beamRA[k] = ra0  # RA in radians
	    beamDEC[k] = dec0    # DEC in radians
	    
	  print "FREQUENCY INTERPOLATION"
	  InterBeamF = np.ndarray(shape=(gxx.shape[0]),dtype=complex)  # generate a random array of elements equal to the length of the antenna gain elements i.e 6643 and in complex number
	  dist = np.abs(freq - np.array(beamFreqs))   # take a frequency(100MHz) and subtract the array of beam frequencies (i.e 100 to 200MHz in steps of 10MHz) to generate an array of different frequencies, then find the absolute value
	  #print "dist=", dist  
	  #print "dist=",len(dist)
	  #print "np.array(beamFreqs)=",len(np.array(beamFreqs))
	  #print "freq=",freq
	  #print  " array-freq =" ,np.array(beamFreqs)
	  
	  ind = np.argsort(dist)      # gives the array of indices that will do the sorting eg. [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10] 
	  ind.flatten()
	  if dist[ind[0]] ==0:         # If the first element from dist is 0, the first element of the antenna gain is the same as the interpolated beam  
	     InterBeamF[:] = gxx[ind[0]]
	  else:                                                   # else, interpolate by multiplying the gain with the weighted average distance as will be seen below
	     weights = dist[ind[0]]**(-2) + dist[ind[1]]**(-2)    # find two values within the array of "dist" find the inverse square and sum them to generate the weight between the two values    
	     InterBeamF[:] = (gxx[ind[0]]*dist[ind[0]]**(-2) + gxx[ind[1]]*dist[ind[1]]**(-2))/weights  # to generate the interpolated beam frequencies  multiply the gain of the antenna by the inverse square of the two values and divide by the weights 
	  #print "gxx[ind[0]]=",gxx[ind[0]]
	  #print "InterBeamF[:]=",InterBeamF[:]
	  #print "gxx[ind[0]]*dist[ind[0]]**(-2)=",gxx[ind[0]]*dist[ind[0]]**(-2)
	  #print "gxx[ind[0]]=",gxx[ind[0]]
	  #print "dist[ind[0]]**(-2)=",dist[ind[0]]**(-2)
	  #print  "gxx[ind[1]]*dist[ind[1]]**(-2)=",gxx[ind[1]]*dist[ind[1]]**(-2)
	     #print   "weights=", weights 
	     #print "dist[ind[0]]**(-2)=", dist[ind[0]]**(-2)
	     #raise
  '''
  #@  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  beamRA = np.ndarray(shape=(phi.shape[0]),dtype=float)   # generate a random array of equal size or length as in phi eg. phi has a length of 6643 hence we generate an array of 6643 elements  
  beamDEC = np.ndarray(shape=(phi.shape[0]),dtype=float)  
  for k in range(beamRA.shape[0]):
    ra0,dec0 = paper.radec_of(phi[k],theta[k])  # convert from Alt Al system to ra and dec 
    beamRA[k] = ra0  # RA in radians
    beamDEC[k] = dec0    # DEC in radians
	    

  InterHealBeam = np.ndarray(shape=(InterBeamF.shape[0], ra.shape[0]),dtype=np.float_) #,dtype=complex)
  print "InterHealBeam.shape=", InterHealBeam.shape
  print "SPATIAL INTERPOLATION"
  c = 0
  for r in range(ra.shape[0]):
    ra[r] = ra[r] - 2*np.pi if ra[r]>(2*np.pi) else ra[r]     # when the ra is greater than 2 pi , subtract 2 pi from that ra value else maintain it   
    dist = np.sqrt((ra[r] - beamRA)**2 + (dec[r] - beamDEC)**2)  # difference between the ra and the beam ra , the difference between the dec and the beamDEC find the distance between them 
    ind = np.argsort(dist)  # gives the array of indices that will do the sorting eg. [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10]
    ind = ind.flatten()

    for p in xrange(InterBeamF.shape[0]):

       InterHealBeam[p, r] = InterBeamF[p, ind[0]]
       ind2 = np.where(InterHealBeam<=0)
       InterHealBeam[ind2] = np.nan
    
    print 'p', p, c
    c +=1
  # np.savez(open('Beam_2-f%.4g_j%.5f.npz'%(freq*1e-6,np.float(jd)),'wb'),beam=InterHealBeam)  # save the spatially interpolated beam in the .npz format
  return InterHealBeam
     
     
if '__name__==__main__':
   
   #
   #jd = np.load('./dee/julian_dates.npz')['jd']
   #print jd  #[60:100]
   
   	
   #	start-time
   start = time.time()
   startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
   print "Start at %s" % startime


   xpol = 'PAPER_FF_X.ffe'
   #nside = 64
   nside = 64
   npix = hp.nside2npix(nside)   
   freqs = np.linspace(100e6,200e6,100)      # generate frequencies from 100 to 200 MHz in steps of 1MHz i.e (start,stop,range of 100)
   #jds=[2455820.455625,2455819.707673611,2455819.9569907407,2455820.1647569444]
   si=pf.open('The_spectral_index_using_landecker_and_haslam1.fits')  # open the spectral index map (the spectral index map derived from the ratio of 408 and 150MHz template map) saved in the fits format
   si=si[0].data
   print si.shape                           
   si= hp.ud_grade(si,nside)                  # upgrade the nside to the desired one

   
   rad2deg = lambda val: val * 180./np.pi    # defining the conversion "scheme" from radians to degrees
   deg2rad = lambda val: val * np.pi/180     # defining the conversion "scheme" from degree to radians
   coord = hp.pix2ang(nside,np.arange(npix))
   theta = 90-(rad2deg(coord[0])) # galactic latitute
   phi = rad2deg(coord[1])   # galactic longitude
   ra = np.ndarray(shape=(npix),dtype=float)  # creating an array containing the right ascension  
   dec = np.ndarray(shape=(npix),dtype=float)  # creating an array containing the declination
   
    # converting galactic coordinate to equatorial coordinates
   iter = 0
   for ii in range(len(theta)):
       eq_coord = ga2equ([phi[ii],theta[ii]])
       ra[ii] = deg2rad(eq_coord[0])  #convert from galactic longitude to right ascension 
       dec[ii] = deg2rad(eq_coord[1]) #convert from galactic latitude to declination
   #jds=[2455820.164814815,2455820.455625,2455819.707673611,2455819.9569907407,2455820.1647569444]
   jds=[2455819.46147,2455819.71131,2455819.96115,2455820.16935]
   
   for jd in jds:
	   jones = interBeam(xpol,[100,200,1],jd,ra,dec)
	   np.save("healpxbeam_j%s" %jd,jones)
   print 'beam', iter
   iter +=1
#   a=np.load("jones.npy")
#b=a.shape
#   hp.mollview(a[100],unit="K",cbar=True,title="",min=0,max=1, notext=True, coord=["G","C"])
#hp.graticule()
#   hp.graticule(linestyle='-')
#   pl.show()

  



# print 'yes !!!'
# hp.mollview(hp.ud_grade(jones[0].real, 128),cmap=pl.get_cmap('jet'), coord=('G', 'C'))
# hp.mollview(hp.ud_grade(jones[1].real, 128),cmap=pl.get_cmap('jet'), coord=('G', 'C'))
# hp.mollview(hp.ud_grade(jones[0].real - jones[1].real, 128),cmap=pl.get_cmap('jet'), coord=('G', 'C'))
# pl.show()


#jds = np.load('./dee/julian_dates.npz')['jd']  # follow the path and take all the Julian dates in thad directory
#counter = 0
#for jd in jds:                                 # for each Julian date 
#     print 'jd %d out of %d' %(jd,len(jds))
#     jones = interBeam(xpol,[100,200,1],jd,ra,dec)   # use the function interBeam and generate beams at freq 100 to 200MHz in steps of 1MHz
#     np.save('interHealBeam_j%s' %jd, jones)         # save all the beams i.e 100 beams corresponding to each julian date


#    print 'jd:', jd


   import os
   directory = 'newbeams4jds'      # create a new directory  with name newbeams
   if not os.path.exists(directory):
            os.mkdir(directory)
   os.system('mv healpxbeam_j* %s' %directory)
   
   print 'ALL INTER-BEAMS SAVED IN ---->> %S DIRECTORY' %directory
  
   #stop-time
   stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
   print "Stop at %s" % stoptime 
   end = time.time()
   elasped_time = (end - start)/3600.0
   print "Total run time: %7.2f hours" % elasped_time


"""
   for freq in freqs:
      _data = np.loadtxt('haslam.dat')  #load the Haslam map       
      _data = hp.ud_grade(_data, nside) # upgrade to the same nside as above
      _data_f150 = _data*(freq/408e6)**si   # generating a simulated map with the given spectra index map(si) at a given frequency 
      #_data_f150 = _data_f150 /_data_f150 
      np.save('mmap_%s' % (freq), _data_f150)   # saving the simlated maps

      #jd = 2455819.46147                        # a Julian date corresponding to LST= 6 Hours,48mins 41.20 seconds
      
     

#The Julian Dates corresponds to 14 to 15 September,2011

    #jds = np.load('julian_dates.npz')['jd'] 
   #for jd in jds:
   #for freq in freqs:
      #_data_f150 = _data*(i/408e6)**si
      for jd in jds:
          jones = interBeam(xpol,freq,jd,ra,dec)                    # equating the interpolated beam function to variable jones 
      #jones = interBeam(xpol,i,jd,ra,dec)
	  #print jones/jones.sum()
          #print max(jones)
          out_data = _data_f150 * jones # multiplying the simulated maps with the antenna gain (i.e the jones to give the simulated sky observation in complex array) as out_data
          outfile_2 = 'map_f%.4g_j%.5f.npz'%(freq*1e-6,np.float(jd))
      #outfile = 'map_f%.4g_j%.5f.npz'%(i*1e-6,np.float(jd))
          np.savez(open(outfile_2,'wb'),beam=jones,data=out_data)     # saving the outfile as a .npz file      
      
    
"""

