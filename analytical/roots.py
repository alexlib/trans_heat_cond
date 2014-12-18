"""
Plot the zeta, Bi equation and its positive roots. For use in the analytical
solution of a 1D sphere, cylinder, or slab transient heat conduction.

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

Bi = 5  # Biot number

# range of z-values from 0 to 15 as 100 samples
# the range at which the zeta, Bi equation (function) will be evaluated 
z = np.linspace(0, 15, num=100)

# Sphere Positive Roots and Zeta, Bi Equation
#------------------------------------------------------------------------------

b = 2   # shape factor where 2 sphere, 1 cylinder, 0 slab

rts_sph = roots(z, b, Bi)   # positive roots of the zeta, Bi equation
print('roots of sphere zeta, Bi equation \n', rts_sph)

fz_sph = funcZetaSph(z, Bi) # evaluate zeta, Bi sphere equation at z-values 

# Cylinder Positive Roots and Zeta, Bi Equation
#------------------------------------------------------------------------------

b = 1   # shape factor where 2 sphere, 1 cylinder, 0 slab

rts_cyl = roots(z, b, Bi)   # positive roots of the zeta, Bi equation
print('roots of cylinder zeta, Bi equation \n', rts_cyl)

fz_cyl = funcZetaCyl(z, Bi) # evaluate zeta, Bi cylinder equation at z-values 

# Slab Positive Roots and Zeta, Bi Equation
#------------------------------------------------------------------------------

b = 0   # shape factor where 2 sphere, 1 cylinder, 0 slab

rts_slab = roots(z, b, Bi)   # positive roots of the zeta, Bi equation
print('roots of slab zeta, Bi equation \n', rts_slab)

fz_slab = funcZetaSlab(z, Bi) # evaluate zeta, Bi slab equation at z-values

# Plot Function and its Positive Roots
#------------------------------------------------------------------------------
    
py.figure(1)
py.plot(z, fz_sph, lw=2)
py.scatter(rts_sph, np.zeros(len(rts_sph)), c='r', s=60, edgecolor='none')
py.axhline(y=0, c='k', ls='--')
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.title('Sphere')
py.ylim([-20, 20])
py.xlim(0)
py.grid()
py.show()

py.figure(2)
py.plot(z, fz_cyl, lw=2)
py.scatter(rts_cyl, np.zeros(len(rts_cyl)), c='r', s=60, edgecolor='none')
py.axhline(y=0, c='k', ls='--')
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.title('Cylinder')
py.ylim([-20, 20])
py.xlim(0)
py.grid()
py.show()

py.figure(3)
py.plot(z, fz_slab, lw=2)
py.scatter(rts_slab, np.zeros(len(rts_slab)), c='r', s=60, edgecolor='none')
py.axhline(y=0, c='k', ls='--')
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.title('Slab')
py.ylim([-20, 20])
py.xlim(0)
py.grid()
py.show()
