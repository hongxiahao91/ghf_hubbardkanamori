import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *
from rham295K import *

#Model Parameter
latt_n  = [1,1,1]
latt    = Latt_class( latt_n )
U_one   = 2.3
J_one   = 0.35
Ntot    = latt.L*2*4   # 2 Ru per unit cell, 4 electrons per Ru.

numberkana = latt.L*2*3 # 4 Ru per unit cell, 3 kanamori per Ru
site_i = []; site_j = []
for i in range(latt.L*2):
    site_i.append( 0+i*3 ); site_i.append( 0+i*3 ); site_i.append( 1+i*3 )
    site_j.append( 1+i*3 ); site_j.append( 2+i*3 ); site_j.append( 2+i*3 )

U1  = [U_one-2*J_one]*latt.L*6
U2  = [U_one-3*J_one]*latt.L*6
J   = [J_one]*latt.L*6

#Set lattice information
latt = Latt_class( latt_n )
latt.write("latt_param")

#Total state is 2*3*latt.L, additional 2 for soc
Tmatrix = np.zeros(( 12 * latt.L, 12 * latt.L), dtype='complex', order='F')
jumpx = [0] if latt.n[0]==1 else [-1,0,1]
jumpy = [0] if latt.n[1]==1 else [-1,0,1]
jumpz = [0] if latt.n[2]==1 else [-1,0,1]

if latt.n[0]==2: jumpx = [0, 1]
if latt.n[1]==2: jumpy = [0, 1]
if latt.n[2]==2: jumpz = [0, 1]

for i in range(latt.L):
    coor_i = latt.coor(i)
    coor_j = [0,0,0]
    for dx in jumpx:
        coor_j[0] = latt.bound( coor_i[0]+dx, latt.n[0]  )
        for dy in jumpy:
            coor_j[1] = latt.bound( coor_i[1]+dy, latt.n[1]  )
            for dz in jumpz:
                coor_j[2] = latt.bound( coor_i[2]+dz, latt.n[2]  )
                j = latt.index( coor_j )
                Tmatrix[ 6*i:6*(i+1), 6*j:6*(j+1) ] += np.array( Hopping[(dx, dy, dz)] )[0:6,0:6]

Tmatrix[ 6*latt.L:12*latt.L, 6*latt.L:12*latt.L ] = Tmatrix[ 0:6*latt.L, 0:6*latt.L ]


f = open("model_param", 'w')
f.write( '{:16d} \n'.format(6*latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format("charge") )
for i in range( 12*latt.L ):
    for j in range( 12*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( 6*latt.L ):
    f.write( '{:26.18e} \n'.format( U_one ) )
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
maxIterateStep            =  1000
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
    subprocess.call( "python {0}/initialOrderParameterCa2RuO4Onelayer.py {1}".format(pythonFileDir, initOrderParamType),shell=True)
