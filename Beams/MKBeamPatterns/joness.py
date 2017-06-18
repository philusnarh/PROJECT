from scipy import io
import numpy as np
from matplotlib import pylab as plt
#D =  io.loadmat("MK_GDSatcom_%d.mat" % 1100)
#th = D["th"].squeeze()
#ph = D["ph"].squeeze()
#JHH = D["Jqh"].squeeze()
#JVV = D["Jpv"].squeeze()
#JHV = D["Jqv"].squeeze()
#JVH = D["Jph"].squeeze()
   
#f = plt.figure()
#for B,m in [(abs(JHH[:,0])**2,'b'), (abs(JVV[:,0])**2,'r')]:
#    plt.plot(th, 10*np.log10(B), m)
#    i3dB = np.searchsorted(-B, -0.5, 'left')
#    hbw = th[i3dB] + (th[i3dB+1]-th[i3dB])*(0.5-B[i3dB])/(B[i3dB+1]-B[i3dB])
#    print(hbw, 1/2 * 1.25*(c/freq/1e6)/13.5 * 180.0/np.pi)
#plt.grid(True)
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

def run_test():
	D =  io.loadmat("MK_GDSatcom_%d.mat" % 1100)
	th = D["th"].squeeze()
	ph = D["ph"].squeeze()
	JHH = D["Jqh"].squeeze()
	JVV = D["Jpv"].squeeze()
	JHV = D["Jqv"].squeeze()
	JVH = D["Jph"].squeeze()
	# Create x and y indices
	'''x = np.linspace(0, 200, 201)
	y = np.linspace(0, 200, 201)
	x, y = np.meshgrid(JHH.real.ravel()[:256], JHH.real.ravel()[:256])

	#create data
	data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, 0, 10)
	# plot twoD_Gaussian data generated above
	plt.figure()
	plt.imshow(data.reshape(201, 201))'''
	#plt.colorbar()
	phi,rho = meshgrid(ph*pi/180,sin(th*pi/180))
	f = figure()
	for d,l in [(JHH,(2,2,1)), (JVV,(2,2,4)), (JHV,(2,2,2)), (JVH,(2,2,3))]:
	    ax = f.add_subplot(*l)
	    cs = contour(rho*cos(phi), rho*sin(phi), 10*log10(d*conjugate(d)), range(0, -37, -2))
	    clabel(cs, range(0, -37, -2))
	plt.show()


if __name__ == '__main__':
    run_test()
