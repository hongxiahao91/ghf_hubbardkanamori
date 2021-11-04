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

if L != 12*latt.L :
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

    n_mean = float(N) / L

    for i in range(latt.L):
        densityMatrix[i*12+0, i*12+0] = n_mean
        densityMatrix[i*12+1, i*12+1] = n_mean
        densityMatrix[i*12+2, i*12+2] = n_mean

        densityMatrix[L+i*12+3,  L+i*12+3]  = n_mean
        densityMatrix[L+i*12+4,  L+i*12+4]  = n_mean
        densityMatrix[L+i*12+5,  L+i*12+5]  = n_mean

elif typ=="afmx":
    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        densityMatrix[i*12+0,  L+i*12+0] = n_mean/2.0;   densityMatrix[L+i*12+0, i*12+0]   = n_mean/2.0
        densityMatrix[i*12+1,  L+i*12+1] = n_mean/2.0;   densityMatrix[L+i*12+1, i*12+1]   = n_mean/2.0
        densityMatrix[i*12+2,  L+i*12+2] = n_mean/2.0;   densityMatrix[L+i*12+2, i*12+2]   = n_mean/2.0

        densityMatrix[i*12+3,  L+i*12+3]  = -n_mean/2.0; densityMatrix[L+i*12+3,  i*12+3]  = -n_mean/2.0
        densityMatrix[i*12+4,  L+i*12+4]  = -n_mean/2.0; densityMatrix[L+i*12+4,  i*12+4]  = -n_mean/2.0
        densityMatrix[i*12+5,  L+i*12+5]  = -n_mean/2.0; densityMatrix[L+i*12+5,  i*12+5]  = -n_mean/2.0

elif typ=="afmxy":
    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        densityMatrix[i*12+0,  L+i*12+0] = n_mean/2.0;   densityMatrix[L+i*12+0, i*12+0]   = n_mean/2.0
        densityMatrix[i*12+1,  L+i*12+1] = n_mean/2.0;   densityMatrix[L+i*12+1, i*12+1]   = n_mean/2.0
        densityMatrix[i*12+2,  L+i*12+2] = n_mean/2.0;   densityMatrix[L+i*12+2, i*12+2]   = n_mean/2.0

        densityMatrix[i*12+3,  L+i*12+3]  = -n_mean/2.0; densityMatrix[L+i*12+3,  i*12+3]  = -n_mean/2.0
        densityMatrix[i*12+4,  L+i*12+4]  = -n_mean/2.0; densityMatrix[L+i*12+4,  i*12+4]  = -n_mean/2.0
        densityMatrix[i*12+5,  L+i*12+5]  = -n_mean/2.0; densityMatrix[L+i*12+5,  i*12+5]  = -n_mean/2.0

        densityMatrix[i*12+6,  L+i*12+6] = 1j*n_mean/2.0;   densityMatrix[L+i*12+6, i*12+6]   = -1j*n_mean/2.0
        densityMatrix[i*12+7,  L+i*12+7] = 1j*n_mean/2.0;   densityMatrix[L+i*12+7, i*12+7]   = -1j*n_mean/2.0
        densityMatrix[i*12+8,  L+i*12+8] = 1j*n_mean/2.0;   densityMatrix[L+i*12+8, i*12+8]   = -1j*n_mean/2.0

        densityMatrix[i*12+9,   L+i*12+9]   = -1j*n_mean/2.0; densityMatrix[L+i*12+9,   i*12+9]   = 1j*n_mean/2.0
        densityMatrix[i*12+10,  L+i*12+10]  = -1j*n_mean/2.0; densityMatrix[L+i*12+10,  i*12+10]  = 1j*n_mean/2.0
        densityMatrix[i*12+11,  L+i*12+11]  = -1j*n_mean/2.0; densityMatrix[L+i*12+11,  i*12+11]  = 1j*n_mean/2.0

densityMatrix = densityMatrix.ravel()

np.savetxt("densityMatrix.dat", np.column_stack(( np.real(densityMatrix), np.imag(densityMatrix) )),fmt=('%26.16e', '%26.16e'), header="{:26d} \n {:26d} {:26d} \n".format(2, L*2, L*2), comments='')
