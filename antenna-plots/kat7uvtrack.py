import numpy as np
from pyrap.tables import table
from matplotlib import pylab as plt


#

#ms = table("/home/narh/ppKAT/KAT7_12h0s_1400MHz1.MS")
ms = table("/home/narh/kk1.MS")
uvw = ms.getcol("UVW")
A1 = ms.getcol("ANTENNA1")
A2 = ms.getcol("ANTENNA2")
#
for j in xrange(np.max(A1)+1):
	for k in xrange(j+1,np.max(A2)+1):												
		uvloc = (A1==j)&(A2==k) 
		u = uvw[uvloc,0]
		v = uvw[uvloc,1]
		plt.plot(u, v, color="r")
		plt.plot(-u, -v, color="b")
        
plt.axis('image')

plt.xlabel("u /m")
plt.ylabel("v /m")
plt.title("KAT7 UV-COVERAGE")
plt.show()
