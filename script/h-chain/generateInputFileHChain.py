import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Model Parameter
latt_n     = [6]
ktwist     = [0]
t1         = 1.0
e1         = 0.0
e2         = 8.0
e3         = 16.0
U          = 8.0
UOne       = 0.70*U
UTwo       = 0.55*U
Ntot       = latt_n[0]
UpDnFlag   = 0    # 0 up=dn, 1 up=conj(dn) ==> different twist

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

#Set Kanamori operator
numberkana = latt.L*3
site_i = []; site_j = []
for i in range(latt.L):
    site_i.append( i       ); site_i.append( i           ); site_i.append( i + latt.L   )
    site_j.append( i+latt.L); site_j.append( i + latt.L*2); site_j.append( i + latt.L*2 )

U1     = [UOne]*numberkana
U2     = [UTwo]*numberkana
J      = [0]*numberkana

Tmatrix = np.zeros((6 * latt.L, 6 * latt.L), dtype='complex', order='F')
#Hopping
for i in range( len(up_K) ):
    Tmatrix[up_i[i] + 0*latt.L, up_j[i] + 0*latt.L] += up_K[i]
    Tmatrix[up_i[i] + 1*latt.L, up_j[i] + 1*latt.L] += up_K[i]
    Tmatrix[up_i[i] + 2*latt.L, up_j[i] + 2*latt.L] += up_K[i]
for i in range( len(dn_K) ):
    Tmatrix[dn_i[i] + 3*latt.L, dn_j[i] + 3*latt.L] += dn_K[i]
    Tmatrix[dn_i[i] + 4*latt.L, dn_j[i] + 4*latt.L] += dn_K[i]
    Tmatrix[dn_i[i] + 5*latt.L, dn_j[i] + 5*latt.L] += dn_K[i]
#Chemical potential
for i in range(latt.L):
    Tmatrix[i+0*latt.L, i+0*latt.L] += e1
    Tmatrix[i+3*latt.L, i+3*latt.L] += e1

    Tmatrix[i+1*latt.L, i+1*latt.L] += e2
    Tmatrix[i+4*latt.L, i+4*latt.L] += e2

    Tmatrix[i+2*latt.L, i+2*latt.L] += e3
    Tmatrix[i+5*latt.L, i+5*latt.L] += e3

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(3*latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format("charge") )
for i in range( 6*latt.L ):
    for j in range( 6*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( 3*latt.L ):
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
annealStep                =  0
relaxMagnitude            =  0.4          # 1.0 fully relax to new order paramter, 0.0 not update
seed                      =  985456376    # 1. read file, 0. random, else is seeds
initOrderParamType        =  "afmx"       # zero, rhfRandom, uhfRandom, ghfRandom, afmz, afmx

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
    subprocess.call( "python {0}/initialOrderParameterThreeBand.py {1}".format(pythonFileDir, initOrderParamType),shell=True)
