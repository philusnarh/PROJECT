import re
import numpy as np

filename = 'asu.tsv'
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
for iter in xrange(45,len(c1)):   
    #print rx.findall(c1[iter])
    ap.append(np.array(map(float, rx.findall(c1[iter]) )))
    #print ap
    #break
data_elt =  np.array(ap).T

print data_elt #.shape
