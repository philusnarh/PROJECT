import h5py
import numpy as np
import scikits.fitting as fit
import ephem
import kittendal as kd
import re
from pyrap.tables import table
import glob
import copy
import os
#
ntimes = 3140 #785
h5file = 'osk_1444074896.h5'

h5 = kd.kittendal('%s' %h5file)
h5.file.close()

corr_products = h5.corrprods

# ++++++++++++ Select number from  string list +++++++++++ 

s1 = [ map(int, re.findall('\d+', x))  for x in corr_products[:,0] ]
print s1
corr_index = np.unique(np.squeeze(np.array(s1)))
print corr_index
#raise

#msloc = '/home/theonarh/Documents/mario/osk_sim/mkk/mt60'
msfil = 'm60'
msloc = '/home/theonarh/Documents/mario/osk_sim/mkk/m/%s' %msfil
msname =  glob.glob('%s/*.MS' %msloc)
ms_sort = sorted(msname)
#

#
ff = h5file
f3 = h5py.File(ff, 'r+') 
vis_data = f3['Data']['correlator_data'].value

#print vis_data.shape
#print vis_data[:785,:305,0,0].shape

for iter, val in enumerate(corr_index):
    vish = []
    visv = []
    m1 = []
    m2 = []
    for ms in xrange(len(ms_sort)):    

        #ms =  '/home/theonarh/Dropbox/newms/kat7_chan0000.MS'
        tab = table("%s" %ms_sort[ms],readonly=False)
        #tab.colnames()
        mdata = tab.getcol('DATA')
        A0 = tab.getcol('ANTENNA1')
        A1 = tab.getcol('ANTENNA2')
        tab.close()
        vis = mdata[(A0==iter)&(A1==iter)]
        print vis.shape
        #break
        if vis[:,0,3].shape[-1] != ntimes:
        	continue
	else:
		vish.append(vis[:,0,0].real)     #horizontal part
		print np.array(vish).shape
		visv.append(vis[:,0,3].real)     #vertical part
		s = vis[:,0,3].real
		print s
		print s.shape
    
    vish = np.array(vish)
    visv = np.array(visv)    
    ss = vish.shape
    print 'ss0:', ss  
    #raise
    vish = vish.reshape(ss[1],ss[0])
    visv = visv.reshape(ss[1],ss[0])
    #print vish.shape
    ss = vish.shape
    
    print '\n >> ss1:', ss
    print '\n >> iter ', iter
    print '\n *************************************************'
    
#     #for iter, val in enumerate(corr_index1):
#     #for val in corr_index1:
    auto1 = h5.lookup['ant%dh' %val,'ant%dh' %val]    
    vis_data[:ss[0],:ss[1],auto1,0] = vish
    auto2 = h5.lookup['ant%dv' %val,'ant%dv' %val]
    vis_data[:ss[0],:ss[1],auto2,0] = visv
    print vis_data.shape
#raise 
    #
#raise
del f3['Data/correlator_data']
f3.create_dataset('Data/correlator_data', data=vis_data)
#w[...] = np.zeros_like(m1)
#print w
f3.close()
os.system('cp %s %s_%s' %(h5file,msfil,h5file))

print 'Done wai !!!!!!!!!!!'

