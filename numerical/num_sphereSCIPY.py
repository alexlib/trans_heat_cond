# Implicit numerical solution of 1D transient heat conduction
# b = 1 cylinder, b = 2 sphere
# convection at surface & no heat of reaction
# see Ozisik1993, Ch.12, pg.459
# tridiagonal solution using SciPy solve_banded function

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# modules
import numpy as np
import scipy.linalg as sp
import matplotlib.pyplot as py
import time

ti = time.clock()   # start time

# Parameters from Papadikis2010a Table 1
# -------------------------------------------------------------------------
rho = 700       # density of wood, kg/m^3
d = 0.035e-2    # wood particle diameter, m
cpw = 1500      # biomass specific heat capacity, J/kg*K
kw = 0.105      # biomass thermal conductivity, W/m*K
h = 375         # heat transfer coefficient, W/m^2*K
Ti = 300        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Implicit Numerical Model where b = 1 cylinder and b = 2 sphere
# -------------------------------------------------------------------------
nt = 1000       # number of time steps
tmax = 0.8      # max time, s
dt = tmax/nt    # time step, s
t = np.arange(0,tmax+dt,dt)

nr = 99     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

b = 2   # run model as a cylinder (b = 1) or as a sphere (b = 2)

alpha = kw/(rho*cpw)    # thermal diffusivity, alfa = kw / rho*cp, m^2/s
Fo = alpha*dt/(dr**2)   # Fourier number, Fo = alfa*dt / dr^2, (-)
Bi = h*dr/kw            # Biot numbmer, Bi = h*dr / kw, (-)

# create array [TT] to store temperature values
# note that row = time step, column = node
TT = np.zeros((len(t), m))
TT[0, :] = Ti   # first row is initial temperature of sphere or cylinder

# Build tridiagonal band array [ab] and initial column vector {bb} 
# -------------------------------------------------------------------------

# create banded matrix to hold tridiagonal matrix
ab = np.zeros((3, m))

# create range for internal nodes
i = np.arange(1, m-1)

# upper diagonal
ab[0, 1] = -2*(1+b)*Fo              # center node T1
ab[0, 2:] = -Fo*(1 + b/(2*(i+1)))   # internal nodes Tm+1

# center diagonal
ab[1, 0] = 1 + 2*(1+b)*Fo                       # center node T0
ab[1, 1:m-1] = 1 + 2*Fo                         # internal nodes Tm
ab[1, m-1] = 1 + 2*Fo*(1 + Bi + (b/(2*m))*Bi)   # surface node Tr

# lower diagonal
ab[2, 0:m-2] = -Fo*(1 - b/(2*(i+1)))    # internal nodes Tm-1
ab[2, m-2] = -2*Fo                      # surface node Tr-1

# create column vector {bb}
bb = np.zeros(m)
bb[0] = Ti
bb[1:m-1] = Ti
bb[m-1] = Ti + 2*Fo*Bi*(1 + b/(2*m))*Tinf

# check: ab and b
#print('ab \n', ab)  # banded matrix of tridiagonal coefficient matrix [A]
#print('bb \n', bb)  # vector {b}

# Solve using scipy.sparse.linalg.lsqr
# -------------------------------------------------------------------------

for i in range(1, nt+1):
    T = sp.solve_banded((1, 1), ab, bb)
    bb = T.copy()
    bb[m-1] = T[m-1] + 2*Fo*Bi*(1 + b/(2*m))*Tinf
    TT[i, :] = T
    
# check: display final T, best if nr is small number like nr=5
#print('T \n', T)

# plot results
py.figure(1)
py.plot(t, TT[:, m-1], '-k', label='surface')
py.plot(t, TT[:, 0], '--k', label='center')
py.axhline(Tinf, color='r', linestyle='-.')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.ylim([Ti-20, Tinf+20])
py.xlim([0, tmax])
py.grid()
py.show()

# elapsed time
print('scipy time', time.clock()-ti, 'seconds')
