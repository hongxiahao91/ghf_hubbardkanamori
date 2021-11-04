import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Model Parameter
latt_n     = [4,3]
ktwist     = [0.012,0.034]
t1         = 1.0
U          = 2.0
mu         = 0.0
Ntot       = 10
UpDnFlag   = 0    # 0 up=dn, 1 up=conj(dn) ==> different twist

numberkana = 2
site_i     = [0, 1]
site_j     = [1, 2]
U1         = [1.5, 1.5]
U2         = [1.0, 1.0]
J          = [0.5, 0.5]

#Set lattice information
latt = Latt_class( latt_n )
latt.write("latt_param")
up_i, up_j, up_K = HubbardNearestNeighborHopping(latt, ktwist, t1)
if UpDnFlag == 0:
    dn_i = up_i; dn_j = up_j; dn_K = up_K
elif UpDnFlag == 1:
    dn_i, dn_j, dn_K = HubbardNearestNeighborHopping(latt, -np.array(ktwist), t1)
else:
    print( "WRONG!!! Do not know UpDnFlag!!!" )
    sys.exit(1)

Tmatrix = np.zeros((2 * latt.L, 2 * latt.L), dtype='complex', order='F')
for i in range( len(up_K) ):
    Tmatrix[up_i[i], up_j[i]] += up_K[i]
for i in range( len(dn_K) ):
    Tmatrix[dn_i[i] + latt.L, dn_j[i] + latt.L] += dn_K[i]


f = open("model_param", 'w')
f.write( '{:16d} \n'.format(latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format("charge") )
for i in range( 2*latt.L ):
    for j in range( 2*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( U ) )
for i in range( numberkana ):
    f.write( '{:16d} \n'.format(site_i[i]) )
for i in range( numberkana ):
    f.write( '{:16d} \n'.format(site_j[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(U1[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(U2[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(J[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(0.0) )
f.close()

#Method Parameter
initialType               =  "readOrderParameter"  # "setFromModel", "readWaveFunction", "readOrderParameter"
convergeType              =  "energy"        # "energy", "orderParameter"
convergeTolerance         =  1e-10
maxIterateStep            =  10000
annealMagnitude           =  0.3
annealStep                =  10
relaxMagnitude            =  0.4          # 1.0 fully relax to new order paramter, 0.0 not update
seed                      =  985456376    # 1. read file, 0. random, else is seeds
initOrderParamType        =  "rhfRandom"       # zero, rhfRandom, uhfRandom, ghfRandom, afmz, afmx

#write method_param
f = open('ghf_param', 'w')
f.write(   '{:>26} \n'.format(initialType       ) )
f.write(   '{:>26} \n'.format(convergeType      ) )
f.write('{:26.18e} \n'.format(convergeTolerance ) )
f.write(   '{:26d} \n'.format(maxIterateStep    ) )
f.write('{:26.18e} \n'.format(annealMagnitude   ) )
f.write(   '{:26d} \n'.format(annealStep        ) )
f.write('{:26.18e} \n'.format(relaxMagnitude    ) )
f.write(   '{:26d} \n'.format(seed              ) )
f.write(   '{:>26} \n'.format(initOrderParamType) )
f.close()

#Run script to generate orderParameter file.
pythonFileDir =  os.path.dirname(os.path.realpath(__file__))
if initOrderParamType is not None:
    subprocess.call( "python {0}/initialOrderParameter.py {1}".format(pythonFileDir, initOrderParamType),shell=True)