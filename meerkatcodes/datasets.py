import sys; sys.path.append('../')
from utilities import *
import glob, katdal, katholog
from holography import load_data
#%load_ext autoreload
%autoreload

path = '/home/asad/data/meerkat/beam/holography/'
h5s = np.sort(glob.glob(path+'*.h5'))
n = len(h5s)
ids = [h5s[i].split('/')[-1] for i in range(n)]

# Load all the datasets using katholog
#d = [katholog.Dataset(h5s[i], 'meerkat') for i in range(n)]

def inspect(h5s):
    global ids, d, n
    # Print list of all datasets with some basic info
    print "There are %i holography datasets available in /net/ike%s :\n"%(n, path)
    print "#, Filename, target, elevation, azimuth, antennae"
    for i in range(n):
        ants = d[i].radialscan_allantenna
        s, t = d[i].scanantennas, d[i].trackantennas
        scanners = [ants[a] for a in s]
        tracker = [ants[a] for a in t]
        el = d[i].env_el[0]
        az = np.mean(d[i].scanaz*(180/np.pi))
        target = d[i].target.name
        print "%i, %s, %s, %.2f, %.2f, %s"%(i, ids[i], target, el, az, ' '.join(scanners))

# Plot the scanning strategy if True
def azel(h5s):
    global ids, d, n
    for i in range(n):
        plt.plot(d[i].scanaz*(180/np.pi), d[i].scanel*(180/np.pi))
        plt.title(ids[i])
        plt.show()

def extract(h5s, i, a, freqs=range(1000, 1010, 1), Npix=256):
    # Read the list and load beamcube from any of the datasets using holography.load_data
    # This will return the beamcube, and also save it into two .fits (real, imag) files
    beamcube = load_data(h5s[i], # Full path of the ith file\
                         ant = a, # the ith scanning antenna [see the list]\
                         freqs = freqs, # frequency list in MHz\
                         Npix = Npix) # Dimension of the beams
