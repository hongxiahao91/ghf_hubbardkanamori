import sys
import os
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *
import random

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

elif typ=="afmx":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        if coor.sum()%2==0:
            densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0
            densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0
            densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0
        else:
            densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.0
            densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.0
            densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.0

elif typ=="afmxy":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        coorxy = coor[0:1]
        coorz  = coor[2]
        if coorz%2==0:
            if coorxy.sum()%2==0:
                densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0
                densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0
            else:
                densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.0
                densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.0
        else:
            if coorxy.sum()%2==0:
                densityMatrix[i*3+0,  L+i*3+0] =  1j*n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -1j*n_mean/2.0
                densityMatrix[i*3+1,  L+i*3+1] =  1j*n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -1j*n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] =  1j*n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -1j*n_mean/2.0
            else:
                densityMatrix[i*3+0,  L+i*3+0] = -1j*n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   =  1j*n_mean/2.0
                densityMatrix[i*3+1,  L+i*3+1] = -1j*n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   =  1j*n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] = -1j*n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   =  1j*n_mean/2.0

elif typ=="colinear":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        if coor[0]%2==0:
            densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0
            densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0
            densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0
        else:
            densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.0
            densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.0
            densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.0

elif typ=="afmxlayer":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        coorx = coor[0]; coorz  = coor[2]
        ranyu = 1+random.random()*0.05
        if coorz%2==0:
            if coorx%2==0:
                densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0*ranyu
                densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0*ranyu
                densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0*ranyu
            else:
                densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.0*ranyu
                densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.0*ranyu
                densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.0*ranyu
        else:
            if coorx%2==0:
                densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.0*ranyu
                densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.0*ranyu
                densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.0*ranyu;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.0*ranyu
            else:
                densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0*ranyu
                densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0*ranyu
                densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0*ranyu;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0*ranyu

elif typ=="afmxylayer":

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean = float(N) / L

    for i in range(latt.L):
        coor=np.array( latt.coor(i) )
        coorx = coor[0]; coorz  = coor[2]
        if coorz%2==0:
            if coorx%2==0:
                densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0 - 1j*n_mean/2.1
                densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0 - 1j*n_mean/2.1
                densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0 - 1j*n_mean/2.1
            else:
                densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.1 + 1j*n_mean/2.0
                densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.1 + 1j*n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.1 + 1j*n_mean/2.0
        else:
            if coorx%2==0:
                densityMatrix[i*3+1,  L+i*3+1] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean/2.1 + 1j*n_mean/2.0
                densityMatrix[i*3+2,  L+i*3+2] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean/2.1 + 1j*n_mean/2.0
                densityMatrix[i*3+0,  L+i*3+0] = -n_mean/2.1 - 1j*n_mean/2.0;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean/2.1 + 1j*n_mean/2.0
            else:
                densityMatrix[i*3+0,  L+i*3+0] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean/2.0 - 1j*n_mean/2.1
                densityMatrix[i*3+1,  L+i*3+1] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean/2.0 - 1j*n_mean/2.1
                densityMatrix[i*3+2,  L+i*3+2] =  n_mean/2.0 + 1j*n_mean/2.1;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean/2.0 - 1j*n_mean/2.1

elif typ=="afmfcc":

    if latt.n[0]%2 !=0 or latt.n[1]%2 !=0 or latt.n[2]%2 !=0:
        print "Error, lattice rank should be even!"
        exit()

    densityMatrix = np.zeros((2*L, 2*L), dtype=complex)

    n_mean_up = 2*float(N)/L  
    n_mean_dn = 2*float(N)/(3*L)
    print n_mean_up, n_mean_dn
    dex0 = np.array([0,0,0]); dex1 = np.array([1,1,1])

    for i in range(latt.L):
        coor_i = latt.coor(i)
        if( np.array_equal( ( coor_i - dex0 ) %2, dex0 ) or np.array_equal( ( coor_i - dex1 ) %2, dex0 ) ):
            print i
            densityMatrix[i*3+0,  L+i*3+0] =  n_mean_up;   densityMatrix[L+i*3+0, i*3+0]   =  n_mean_up
            densityMatrix[i*3+1,  L+i*3+1] =  n_mean_up;   densityMatrix[L+i*3+1, i*3+1]   =  n_mean_up
            densityMatrix[i*3+2,  L+i*3+2] =  n_mean_up;   densityMatrix[L+i*3+2, i*3+2]   =  n_mean_up
        else:
            densityMatrix[i*3+0,  L+i*3+0] = -n_mean_dn;   densityMatrix[L+i*3+0, i*3+0]   = -n_mean_dn
            densityMatrix[i*3+1,  L+i*3+1] = -n_mean_dn;   densityMatrix[L+i*3+1, i*3+1]   = -n_mean_dn
            densityMatrix[i*3+2,  L+i*3+2] = -n_mean_dn;   densityMatrix[L+i*3+2, i*3+2]   = -n_mean_dn

densityMatrix = densityMatrix.ravel()

np.savetxt("densityMatrix.dat", np.column_stack(( np.real(densityMatrix), np.imag(densityMatrix) )),fmt=('%26.16e', '%26.16e'), header="{:26d} \n {:26d} {:26d} \n".format(2, L*2, L*2), comments='')
