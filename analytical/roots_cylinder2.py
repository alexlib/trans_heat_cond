# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# modules

import numpy as np
import scipy.optimize as op
import scipy.special as sp
import matplotlib.pyplot as py

py.close('all')

# check the Bessel functions with Incropera2011, Table B.4, pg.1017
 
for i in np.arange(0, 2.5, 0.1):
    if i == 0:
        print('step\t J0\t J1')
    j0 = sp.j0(i)
    j1 = sp.j1(i)
    print('{:.4f}\t {:.4f}\t {:.4f}'.format(i, j0, j1))


# function f(z) = z*J1(z) - Bi*J0(z)

def funcZeta(z):
    Bi = 5
    f = z*sp.j1(z) - Bi*sp.j0(z)
    return f


# solve function for a range of z-values

z = np.linspace(0, 15, num=100) # roots from 0 to 15 for 100 samples
fz = funcZeta(z)                # evaluate function a z values

# find roots of the function

sign = np.sign(fz)      # sign change of values in fz
diff = np.diff(sign)    # calculate difference between fz values
where = np.where(diff)  # determine index where values change, returns a tuple

roots = np.zeros(len(where[0])) # setup array to store roots

for i, j in enumerate(where[0]):
    roots[i] = op.brentq(funcZeta, z[j], z[j+1]) # roots for a given interval

print('where\n', where)
print('roots\n', roots)

# plot function

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
