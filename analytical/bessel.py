"""
Example of using SciPy for bessel functions of the first kind. The calculated
values are printed to the console.

Reference:
Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304, pg.1017
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import scipy.special as sp

# Bessel Function Example
#------------------------------------------------------------------------------

# header for the data printed in the console
print('step\t J0\t J1')

# calculate and print the Bessel function values in a given range
# the values should agree with Bergman 2011, Table B.4, pg.1017
for i in np.arange(0, 2.5, 0.1):
    j0 = sp.j0(i)
    j1 = sp.j1(i)
    print('{:.4f}\t {:.4f}\t {:.4f}'.format(i, j0, j1))
