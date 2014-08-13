# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries and packages
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as py

# Biot number,  Bi = h*r / k
Bi = 5

# range of roots
step = 0.1
zi = 0
zf = 15 + step  # add step to end on integer

# function f(z) = 1 - z*cot(z) - Bi
# note that cot = 1/tan
def funcZeta(z):
    f = 1 - z*(1 / np.tan(z)) - Bi
    return f

# find sign change and roots of function
f = []
root = []

for i in np.arange(zi, zf, step):
    f.append(funcZeta(i))

k = 0
s = 0

for i in range(0, len(f)):
    print(i, k, s, f[k])
    k = k + 1
    s = s + 0.1
    if f[k-1] < 0 and f[k] > 0:
        root.append(opt.brentq(funcZeta, s-0.1, s))
        print('* sign change')

# print all positive roots
for i in root:
    print('root = %.4f' % i)

# plot
z = np.arange(zi, zf, step)
fz = funcZeta(z)

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
