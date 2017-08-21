#!/usr/bin/python
import cPickle,pylab as plt
import gaussfitsutil
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-clip-0.2.cp', 'rb') as input:
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-holobm2-clip-0.6.cp', 'rb') as input:
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat_scaled-clip-0.6.cp', 'rb') as input: 
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat_scaled_down-clip-0.2.cp', 'rb') as input: 
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat-em-clip-0.2.cp', 'rb') as input: 
#with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat_mb1_scaled_down_to_em_grid-clip-0.2.cp', 'rb') as input: 
with open('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat_mb2_scaled_down_to_em_grid-clip-0.2.cp', 'rb') as input: 
        xx,yy = cPickle.load(input)
#xx,yy = cPickle.load(file('/home/iheanetu/meerkat_beams/holography/beamfits-mearkat_mb1_scaled_down_to_em_grid-clip-0.2.cp'))

count= 1  #   
plt.figure(figsize=(20,20))
#nchan = len(xx[xx.keys()[0].x0)]
#freq = [(900 +i)*0.001 for i in range(nchan)

for ant in xx.keys()[:]: # ['00', '12', '17'][1:2] :
    VAR_THRESHOLD = 0.03
    #ant = '00'
    xx[ant].set_mask((xx[ant].var>VAR_THRESHOLD)|(yy[ant].var<=0))
    # xx[ant].set_mask((xx[ant].var>VAR_THRESHOLD)|(xx[ant].var<=0))
    xx[ant].var[xx[ant].var<=0] = 0
    # xx[ant].var[xx[ant].var<=0] = 0

    plt.subplot(3,2,count)
    #plt.plot(xx[ant].xw)
    #plt.plot(xx[ant].yw)
    plt.plot(xx[ant].xw.data)
    plt.plot(xx[ant].yw.data)
    plt.ylim(1.2,2.3)
    #plt.ylim(1.2,1.6)
    #plt.semilogy()
#     plt.ylim(1.2,2.3)
    count += 1
    plt.subplot(3,2,count)
    plt.plot(xx[ant].x0)
    plt.plot(xx[ant].y0)
    #plt.plot(xx[ant].x0.data)
    #plt.plot(xx[ant].y0.data)
    #plt.semilogy()
    # plt.ylim(-35,5)# 
#     plt.ylim(8.6,8.7)
    count += 1
plt.show()

