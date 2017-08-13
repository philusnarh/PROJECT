#
import os
#%matplotlib inline
# from zernike2 import *
from astropy.io import fits
import numpy as np
#


def interplotateMissingFreqsbeams(data,dataStartFreq,dataEndFreq,newbeamsFreqStep):
    
#     nchan, xpl = np.array(data).shape             # Get shape of the data
    nchan = np.array(data).shape[0] 
    
    #sfreq = np.array([dataStartFreq + i*dataStepFreq for i in range(nchan)])  # Calc. source data freqs.
    sfreq = np.linspace(dataStartFreq, dataEndFreq, nchan)  # Calc. source data freqs.  
    
    new_chan_numb = (dataEndFreq - dataStartFreq)/newbeamsFreqStep
    nfreq = np.arange(sfreq[0],sfreq[-1],newbeamsFreqStep)  # Calc. new beams freqs. # np.linspace(900e+6,1300e+6,7)
    print len(nfreq)
    #print sfreq[0]
    #print sfreq[1]
    #print newbeamsFreqStep
    #raise 
    
#    raise
    new_beams = []                                    # variable for storing new beams
    
    for ifreq, freq in zip(range(len(nfreq)), nfreq):      # for each new beam ...
	#print 'ifreq:', ifreq

        index = {}                                 # variable for storing freq. diff. indexes
        
        freq_diff = abs(sfreq-freq) 		    # Calc. the frequency difference
        #print "freq_diff=",len(freq_diff)
#  	print ifreq
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
#            print 'yeh', ifreq
#    print 'yeh', ifreq
    print np.array(new_beams).shape
    return np.array(new_beams)
    

for n in [1, 2, 3]:
	if n == 1:
		print 'ant 1 started'
		path = '/home/twum/theo/meerkat_beams/holography'
		s1 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_re.fits'
		s2 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_re.fits'
		s3 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_re.fits'
		s4 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_re.fits'

		s5 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_im.fits'
		s6 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_im.fits'
		s7 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_im.fits'
		s8 = '1487813282_m00_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_im.fits'
	##
	elif n == 3:
		print 'ant 17 started'
		path = '/home/twum/theo/meerkat_beams/holography'
		s1 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_re.fits'
		s2 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_re.fits'
		s3 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_re.fits'
		s4 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_re.fits'

		s5 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_im.fits'
		s6 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_im.fits'
		s7 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_im.fits'
		s8 = '1487813282_m017_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_im.fits'
	##
        else:
        	print 'ant 12 started'
		path = '/home/twum/theo/meerkat_beams/holography'
		s1 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_re.fits'
		s2 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_re.fits'
		s3 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_re.fits'
		s4 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_re.fits'

		s5 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XX_im.fits'
		s6 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_XY_im.fits'
		s7 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YX_im.fits'
		s8 = '1487813282_m012_900MHz_1MHz_500channels_Jones_mb1_scaled_up_to_em_grid_YY_im.fits'

	ss = [s1, s2, s3, s4, s5, s6, s7, s8]

	jns = ['XX_re', 'XY_re', 'YX_re', 'YY_re', 'XX_im', 'XY_im', 'YX_im', 'YY_im' ]

	for jj in xrange(len(ss)):

		f1 = os.path.join(path, '%s' %ss[jj])
		hdr = fits.open('%s' %f1)[0].header
		hdr['CDELT3'] = 10e6
	#	m = [fits.open('%s' %f1)[0].data[i] for i in np.linspace(0,499,50).astype(int)]
		if  (n==1) & (jns[jj] == jns[0]):
			m = [fits.open('%s' %f1)[0].data[i] for i in np.linspace(0,499,500).astype(int)  \
			     if  ((fits.open('%s' %f1)[0].data[i]).max() > 0.98) & ((fits.open('%s' %f1)[0].data[i]).max() < 1.004)]
			     
		        ind = [i for i in np.linspace(0,499,500).astype(int)  \
			     if  ((fits.open('%s' %f1)[0].data[i]).max() > 0.98) & ((fits.open('%s' %f1)[0].data[i]).max() < 1.004)]
	     	else: m = fits.open('%s' %f1)[0].data[ind]; print 'YEBO YEBO'
	
		fits.writeto('hb_meerkat%s_%s.fits' %(str(n).zfill(3), jns[jj]), abs(np.array(m)[:13]),header = hdr, clobber=True)
	
		dt = fits.open('hb_meerkat%s_%s.fits' %(str(n).zfill(3), jns[jj]))[0].data
		hd = fits.open('hb_meerkat%s_%s.fits' %(str(n).zfill(3), jns[jj]))[0].header
		data = np.nan_to_num(dt)
		dataStartFreq = hd['CRVAL3']
		dataEndFreq = 3e8/0.21
		newbeamsFreqStep = hdr['CDELT3'] = 10e6
		interdt = interplotateMissingFreqsbeams(data,dataStartFreq,dataEndFreq,newbeamsFreqStep)
#		print jj
		fits.writeto('hbinterp_meerkat%s_%s.fits' %(str(n).zfill(3), jns[jj]),np.array(interdt),header = hd, clobber=True)
		print 'yoo', jns[jj]


