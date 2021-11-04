import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *

def readHopping( fullerFile ):
    f = open( fullerFile+".txt", 'r')

    fullAll = f.read()
    fullAll = fullAll.split()

    if( len(fullAll) != 1500 ):
        print( "Error!!! Fullerenes file read length is wrong!" )
        exit()

    Hopping={}
    for i in range(125):
        dist    = tuple( map(int, fullAll[(0+i*12):(3+i*12)] ) )
        matrix  = np.array( map(float, fullAll[(3+i*12):(12+i*12)] ) ).reshape(3,3)
        Hopping[dist] = matrix
    f.close()

    return Hopping

#Model Parameter
fullerFile = "Cs3C60_804"
latt_n  = [4,2,2]
latt    = Latt_class( latt_n )
U_one   = 0.5
J_one   = -0.04
Ntot    = latt.L*3   # 4 electrons per fullerene.

gen_path= os.path.dirname(os.path.realpath(__file__))
Hopping = readHopping( gen_path+'/'+fullerFile )

numberkana = latt.L*3 # 3 kanamori per fullerene
site_i = []; site_j = []
for i in range(latt.L):
    site_i.append( 0+i*3 ); site_i.append( 0+i*3 ); site_i.append( 1+i*3 )
    site_j.append( 1+i*3 ); site_j.append( 2+i*3 ); site_j.append( 2+i*3 )

U1  = [U_one-2*J_one]*latt.L*3
U2  = [U_one-3*J_one]*latt.L*3
J   = [J_one]*latt.L*3

#Set lattice information
latt = Latt_class( latt_n )
latt.write("latt_param")

#Total state is 3*latt.L, additional 2 for soc
Tmatrix = np.zeros( (6*latt.L, 6*latt.L), dtype='complex', order='F')
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
                Tmatrix[ 3*i:3*(i+1), 3*j:3*(j+1) ] += Hopping[(dx, dy, dz)]

Tmatrix[ 3*latt.L:6*latt.L, 3*latt.L:6*latt.L ] = Tmatrix[ 0:3*latt.L, 0:3*latt.L ]

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(3*latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format("charge") )
for i in range( 6*latt.L ):
    for j in range( 6*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( 3*latt.L ):
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
relaxMagnitude            =  0.28         # 1.0 fully relax to new order paramter, 0.0 not update
seed                      =  985456376    # 1. read file, 0. random, else is seeds
initOrderParamType        =  "afmfcc"       # zero, rhfRandom, uhfRandom, ghfRandom, afmx, afmxy, colinear, afmxlayer, afmxylayer, afmfcc

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
    subprocess.call( "python {0}/initialOrderParameterFullerenes.py {1}".format(pythonFileDir, initOrderParamType),shell=True)
