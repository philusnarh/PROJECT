import ephem
import numpy as np
az  = -35.000262416  #deg
el  =  61.0 #51.0910024663 #deg
lat = -30.7214 #deg
lon = 21.4106 #deg
alt = 1038  # m
#ut  = 2455822.20000367 #julian date
date =  "2017/02/23 03:28:04.182" # '2000/01/01 12:00:00'  
# Which Julian Date does Ephem start its own count at?
#J0 = ephem.julian_date(0)

meerkat_observer = ephem.Observer()
meerkat_observer.lon = lon  
meerkat_observer.lat = lat  
meerkat_observer.elevation = alt
meerkat_observer.date = date
meerkat_observer.epoch = ephem.J2000
#meerkat_observer.pressure = 0
#ut - J0

#print observer.date
print meerkat_observer.radec_of(float(ephem.degrees(az)),
	float(ephem.degrees(el)))


