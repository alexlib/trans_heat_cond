# Calculate roots for sphere, cylinder, and slab analytical solution
# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# modules

import numpy as np
import scipy.special as sp
import matplotlib.pyplot as py
from funcRoots import roots

py.close('all')

# Sphere Roots Function: f(z) = 1 - z*cot(z) - Bi
# note that cot = 1/tan

def funcZetaSph(z):
    Bi = 5
    f = 1 - z*(1 / np.tan(z)) - Bi
    return f


# Cylinder Roots Function: f(z) = z*J1(z) - Bi*J0(z)

def funcZetaCyl(z):
    Bi = 5
    f = z*sp.j1(z) - Bi*sp.j0(z)
    return f
    

# Slab Roots Function: f(z) = z*tan(z) - Bi

def funcZetaSlab(z):
    Bi = 5
    f = z*np.tan(z) - Bi
    return f
    

# Define shape factor and range of values to evaluate roots function    
b = 0
z = np.linspace(0, 15, num=100) # range of z-values from 0 to 15 as 100 samples

if b==2:
    fz = funcZetaSph(z)     # evaluate sphere function at z-values
    f = funcZetaSph         # declare sphere function as roots function
elif b==1:
    fz = funcZetaCyl(z)     # evaluate cylinder function at z-values
    f = funcZetaCyl         # declare cylinder function as roots function
elif b==0:
    fz = funcZetaSlab(z)
    f = funcZetaSlab

# Calculate positive roots of the function for a range of z-values
roots = roots(z, fz, f, b)
print('roots\n', roots)

# plot function and roots
py.figure(1)
py.plot(z, fz, lw=2)
py.scatter(roots, np.zeros(len(roots)), c='r', s=60, edgecolor='none')
py.axhline(y=0, c='k', ls='--')
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.ylim([-20, 20])
py.xlim(0)
py.grid()
py.show()
