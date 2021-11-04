import sys
import os
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *

if len(sys.argv) < 2:
    sys.exit("Missing arguments!!! Example: python initialOrderParameter.py rhfRandom")
typ = sys.argv[1]

latt  = Latt_class.read("latt_param")

f =  open("model_param", 'r')
L = int( f.readline() )
N = int( f.readline() )
f.close()

if L != 3*latt.L :
    sys.exit("Lattice size is not consistent in latt_param and model_param!")

if typ=="zero":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

elif typ=="rhfRandom":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)
    htmp = np.random.rand(L,L) + 1j*np.random.rand(L,L)
    htmp = ( htmp + np.conj(htmp.T) ) / 2.0
    densityMatrix[0:L,   0:L  ] = htmp
    densityMatrix[L:2*L, L:2*L] = htmp

elif typ=="uhfRandom":
    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    htmp = np.random.rand(L,L) + 1j*np.random.rand(L,L)
    htmp = ( htmp + np.conj(htmp.T) ) / 2.0
    densityMatrix[0:L,   0:L  ] = htmp

    htmp = np.random.rand(L,L) + 1j*np.random.rand(L,L)
    htmp = ( htmp + np.conj(htmp.T) ) / 2.0
    densityMatrix[L:2*L, L:2*L] = htmp

elif typ=="ghfRandom":

    htmp = np.random.rand(2*L,2*L) + 1j*np.random.rand(2*L,2*L)
    densityMatrix = ( htmp + np.conj(htmp.T) ) / 2.0

elif typ=="afmz":
    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / ( (L/3.) *4 )

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        if coor.sum()%2==0:
           densityMatrix[i,          i]               = n_mean
           densityMatrix[i+latt.L,   i+latt.L]        = n_mean
           densityMatrix[i+latt.L*2, i+latt.L*2]      = n_mean

           densityMatrix[i+L,          i+L]           = n_mean
           densityMatrix[i+latt.L+L,   i+latt.L+L]    = 0
           densityMatrix[i+latt.L*2+L, i+latt.L*2+L]  = 0
        else:
           densityMatrix[i,          i]               = n_mean
           densityMatrix[i+latt.L,   i+latt.L]        = 0
           densityMatrix[i+latt.L*2, i+latt.L*2]      = 0

           densityMatrix[i+L,          i+L]           = n_mean
           densityMatrix[i+latt.L+L,   i+latt.L+L]    = n_mean
           densityMatrix[i+latt.L*2+L, i+latt.L*2+L]  = n_mean

elif typ=="afmx":
    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / ( (L/3.) *4 )

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        if coor.sum()%2==0:
           densityMatrix[i,           i+L         ] = 0;          densityMatrix[i+L,          i         ] = 0;
           densityMatrix[i+latt.L,    i+latt.L+L  ] = n_mean/2;   densityMatrix[i+latt.L+L,   i+latt.L  ] = n_mean/2;
           densityMatrix[i+latt.L*2,  i+latt.L*2+L] = n_mean/2;   densityMatrix[i+latt.L*2+L, i+latt.L*2] = n_mean/2;
        else:
           densityMatrix[i,           i+L         ] = 0;          densityMatrix[i+L,          i         ] = 0;
           densityMatrix[i+latt.L,    i+latt.L+L  ] =-n_mean/2;   densityMatrix[i+latt.L+L,   i+latt.L  ] =-n_mean/2;
           densityMatrix[i+latt.L*2,  i+latt.L*2+L] =-n_mean/2;   densityMatrix[i+latt.L*2+L, i+latt.L*2] =-n_mean/2;

densityMatrix = densityMatrix.ravel()

np.savetxt("densityMatrix.dat", np.column_stack(( np.real(densityMatrix), np.imag(densityMatrix) )),fmt=('%26.16e', '%26.16e'), header="{:26d} \n {:26d} {:26d} \n".format(2, L*2, L*2), comments='')
