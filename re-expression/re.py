import re
import numpy as np

filename = '/home/theonarh/Desktop/K7s.txt'
data = open("%s" %filename)
dat = data.read()
lst = dat.splitlines()
c1 = np.array(lst)

#  ++++++++++++++ extract a float/ int number from a string in Python  ++++++++++++
#
numeric_const_pattern = r"""
                   [-+]? # optional sign
                   (?:
                 (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                       |
                       (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
                  )
                    # followed by optional exponent part if desired
                           (?: [Ee] [+-]? \d+ ) ?
              
                    """
                   
rx = re.compile(numeric_const_pattern, re.VERBOSE)
#
print '\n >> Generating Numerically-defined antenna element pattern data '

ap = []
for iter in xrange(4,len(c1)):   
    #print rx.findall(c1[iter])
    ap.append(np.array(map(float, rx.findall(c1[iter]) )))
    #print ap
    #break
data_elt =  np.array(ap).T

#    ++++++ Rescale theta +++++++

theta_deg =  12/13.5 * data_elt[0,:]
phi_deg = data_elt[1,:]   
abs_theta = np.sqrt(data_elt[2,:]**2 + data_elt[3,:]**2)
phase_theta_deg = np.arctan2(data_elt[3,:], data_elt[2,:])
abs_phi = np.sqrt(data_elt[4,:]**2 + data_elt[5,:]**2)
phase_phi_deg = np.arctan2(data_elt[5,:], data_elt[4,:])

print '\n >> Generating Embedded Element Pattern Files'

# +++++ header ++++++++
h = ('{a:^8}{b:^8}{c:^8}{d:^8}{e:^8}{f:^8}{g:^8}{h:^8}'.format(a='Theta [deg.]', b='Phi[deg.]', c=' Abs(E   )[ ]',
d='Abs(Theta)[ ] ', e='Phase(Theta)[deg.]', f='Abs(Phi  )[      ]', g='Phase(Phi  )[deg.]', h='Ax.Ratio[      ] '))

np.savetxt('dipole_CST.txt' ,np.array([theta_deg, phi_deg, np.zeros(len(phi_deg)), abs_theta, phase_theta_deg,
          abs_phi, phase_phi_deg, np.zeros(len(phi_deg))]).T, delimiter = ',', header= h )


print 'Done !!!!!'
