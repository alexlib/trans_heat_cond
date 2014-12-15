# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# modules

import numpy as np
import scipy.optimize as op
import scipy.special as sp
import matplotlib.pyplot as py

py.close('all')

# function f(z) = 1 - z*cot(z) - Bi
# note that cot = 1/tan

def funcZeta(z):
    Bi = 5
    f = 1 - z*(1 / np.tan(z)) - Bi
    return f


# solve function for a range of z-values

z = np.linspace(0, 15, num=100) # roots from 0 to 15 for 100 samples
fz = funcZeta(z)                # evaluate function a z values

# find roots of the function

sign = np.sign(fz)          # sign change of values in fz
diff = np.diff(sign)        # calculate difference between fz values
where = np.where(diff>0)    # index of positive value change, returns a tuple

roots = np.zeros(len(where[0])) # setup array to store roots

for i, j in enumerate(where[0]):
    roots[i] = op.brentq(funcZeta, z[j], z[j+1]) # roots for a given interval

print('where\n', where)
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
