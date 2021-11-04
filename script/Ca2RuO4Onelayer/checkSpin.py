import numpy as np  

Sx, Sy = np.loadtxt("Splus_1_Average.dat", usecols=(0,1), unpack=True)
nup = np.loadtxt("Nup_1_Average.dat",usecols=(0,), unpack=True)
ndn = np.loadtxt("Ndn_1_Average.dat",usecols=(0,), unpack=True)
Sz = (nup - ndn)/2.0

print( "Nup-Ndn={}".format( np.sum(nup)-np.sum(ndn) ) )
print( "SxTot={}, SyTot={}, SzTot={}".format( np.sum(Sx), np.sum(Sy), np.sum(Sz) ) )
print( "{:<3}  {:<9} {:<9} {:<9}".format("Ru", "Sx", "Sy", "Sz") )
for i in range(0, len(Sz), 3):
    print( "{:<3} {:< 5f} {:< 5f} {:< 5f}".format( i//3, np.sum( Sx[i:i+3] ), np.sum( Sy[i:i+3] ), np.sum( Sz[i:i+3] ) )  )
