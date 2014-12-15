# Roots of the function f(z) = 1 - z*cot(z) - Bi
# where z is the positive root

# modules

import numpy as np
import scipy.optimize as op

def roots(z, fz, func, b):
    
    sign = np.sign(fz)          # sign change of values in fz
    diff = np.diff(sign)        # calculate difference between fz values
    
    if b == 2 or b==0:
        where = np.where(diff>0)    # index of positive value change, returns a tuple
        print('sphere or slab')
    elif b==1:
        where = np.where(diff)
        print('cylinder')
    
    roots = np.zeros(len(where[0])) # setup array to store roots
    
    for i, j in enumerate(where[0]):
        roots[i] = op.brentq(func, z[j], z[j+1]) # roots for a given interval
        
    return roots

