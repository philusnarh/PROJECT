
import numpy as np  
f = np.loadtxt('s15degedd.txt')  
flux_ind = np.where(f[:,2]<= f[:,2][1]) 
new_ra = f[:,0][flux_ind[0]] 
new_ra.shape
new_dec = f[:,1][flux_ind[0]]
flux = f[:,2][flux_ind[0]] 
new_ar = np.array([new_ra, new_dec, flux])
np.savetxt('sums15srcs.txt', new_ar.T, delimiter=' ',fmt='%s')

np.savetxt('sums15srcs.txt', new_ar.T, delimiter=' ',fmt='%.4e')

np.savetxt('sums15srcs.txt', new_ar.T, delimiter=' ',fmt='%s')

np.savetxt('sums15srcs.txt', new_ar.T, delimiter=' ',fmt='%.2e')
