"""
Plot the zeta, Bi equation and its positive roots.
For use in the analytical solution of a 1D sphere, cylinder, or slab transient 
heat conduction.

References: 
1) Recktenwald 2006
2) Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304
"""

# Modules and Required Functions
#------------------------------------------------------------------------------
 
import numpy as np
import matplotlib.pyplot as py
from funcRoots import roots
from funcZeta import funcZetaSph, funcZetaCyl, funcZetaSlab

py.close('all')

# Parameters
#------------------------------------------------------------------------------
b = 2   # shape factor where 2 sphere, 1 cylinder, 0 slab
Bi = 5  # Biot number

# range of z-values from 0 to 15 as 100 samples
# the range at which the zeta, Bi equation (function) will be evaluated 
z = np.linspace(0, 15, num=100)

# Positive Roots of the zeta, Bi Equation for Sphere, Cylinder, or Slab
#------------------------------------------------------------------------------

roots = roots(z, b, Bi)     # positive roots of the zeta, Bi equation
print('roots\n', roots)     # print roots to console

# Plot Function and its Positive Roots
#------------------------------------------------------------------------------
if b==2:
    fz = funcZetaSph(z, Bi)     # evaluate sphere function at z-values    
    title = "Sphere"            # title for plot
elif b==1:
    fz = funcZetaCyl(z, Bi)     # evaluate cylinder function at z-values
    title = "Cylinder"          # title for plot
elif b==0:
    fz = funcZetaSlab(z, Bi)    # evaluate slab function at the z-values
    title = "Slab"              # title for plot
    
py.figure(1)
py.plot(z, fz, lw=2)
py.scatter(roots, np.zeros(len(roots)), c='r', s=60, edgecolor='none')
py.axhline(y=0, c='k', ls='--')
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.title(title)
py.ylim([-20, 20])
py.xlim(0)
py.grid()
py.show()
