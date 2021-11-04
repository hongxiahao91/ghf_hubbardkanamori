import numpy as np
import sys

if len(sys.argv) < 2:
    sys.exit("Missing arguments!!! Example: python split_three.py filename")
filename = sys.argv[1]

data = np.loadtxt( filename,  unpack=True, usecols=(0,) )
L = len(data)
lattL =  L/3
if( lattL*3 != L ):
   print "data length problem:", L, lattL

np.savetxt("3_"+filename, np.column_stack( (data[0:lattL], data[lattL:lattL*2], data[lattL*2:lattL*3]) ), fmt=('%26.16e','%26.16e','%26.16e') )
