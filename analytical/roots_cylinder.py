# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries and packages
import numpy as np
import scipy.optimize as opt
import scipy.special as sp
import matplotlib.pyplot as py

# check the Bessel functions with Incropera2011, Table B.4, pg.1017 
for i in np.arange(0, 2.5, 0.1):
    j0 = sp.j0(i)
    j1 = sp.j1(i)
    print('{:.4f} - {:.4f} - {:.4f}'.format(i, j0, j1))

# range of roots
step = 0.1
zi = 0
zf = 15 + step

z = np.arange(zi, zf, step)

# function f(z) = z*J1(z) - Bi*J0(z)
Bi = 5

def funcZeta(z):
    f = z*sp.j1(z) - Bi*sp.j0(z)
    return f
    
fz = funcZeta(z)

# find sign change of function f(z) and positive roots
root = []

k = 0
s = 0

for i in z:
    k = k + 1
    s = s + 0.1
    if k < len(z) and fz[k-1] < 0 and fz[k] > 0:
        root.append(opt.brentq(funcZeta, s-0.1, s))
        print('s\t', s)        
    elif k < len(z) and fz[k-1] > 0 and fz[k] < 0:
        print('s\t', s)
        root.append(opt.brentq(funcZeta, s-0.1, s))

print('roots \n', root)

# plot
py.figure(1)
py.plot(z, fz)
py.ylabel('f ($\zeta$)')
py.xlabel('$\zeta$')
py.ylim([-20, 20])
py.axhline(y=0,color='k',linestyle='--')
py.rcParams['xtick.major.pad'] = 6
py.rcParams['ytick.major.pad'] = 6
py.grid()
py.show()
