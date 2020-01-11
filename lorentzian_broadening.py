'''
Please input the "filepath".
This program can give the Lorentzian broadening spectrum.
By mhf.


Date: 202001108

'''
#=========================================================================
# Frequency: the first dimension of usecols
# Intensity: the second dimension of usecols.
# Please modify them!
# output_array: column 0 is frequency, column 1 is intensity
#==============================================================================
import numpy as np

filepath=("morse.dat")
fre,intensity=np.loadtxt(filepath,comments='#',converters=None,usecols=(0,4),unpack=True)

x_min=0.5
x_max=3.0
Npoints=1000
g=0.09

X=np.linspace(x_min, x_max, Npoints)
Y=np.zeros(Npoints)

for x0,y0 in zip(fre,intensity):
    Y= Y + y0*g/(2*np.pi)/((X-x0)**2+0.25*g**2)
output_array=np.zeros([Npoints,2])
output_array[:,0]=X
output_array[:,1]=Y


np.savetxt('broadened.dat',output_array,fmt='%.8f')