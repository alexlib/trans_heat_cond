"""
Example of using SciPy for bessel functions of the first kind.

References:
Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304, pg.1017
"""

# modules

import numpy as np
import scipy.special as sp


# check the Bessel functions with Incropera 2011, Table B.4, pg.1017
 
for i in np.arange(0, 2.5, 0.1):
    if i == 0:
        print('step\t J0\t J1')
    j0 = sp.j0(i)
    j1 = sp.j1(i)
    print('{:.4f}\t {:.4f}\t {:.4f}'.format(i, j0, j1))
