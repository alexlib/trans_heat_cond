# Implicit numerical solution of 1D transient heat conduction
# b = 1 cylinder, b = 2 sphere
# convection at surface & no heat of reaction
# see Ozisik1993, Ch.9
# nonlinear version, variable k term, d/dr(k*r^2*dT/dr)
# just a check to make sure it matches num_sphere.py results

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# modules
import numpy as np
import matplotlib.pyplot as py

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
TT = np.zeros((1,m))

TT[0, :] = Ti   # first row is initial temperature of sphere or cylinder

# build coefficient matrix [A] and initial column vector {C}
A = np.zeros((m,m))     # pre-allocate [A] array
C = np.zeros((m,1))     # pre-allocate {C} vector

v = dt / (rho*cpw)

A[0,0] = 1 + (2*v*kw*(1+b)) / (dr**2)
A[0,1] = -(2*v*kw*(1+b)) / (dr**2)
C[0,0] = Ti
    
for k in range(1,m-1):
        
    ri = (k*dr)**b
    rminus12 = ((k-0.5)*dr)**b
    rplus12 = ((k+0.5)*dr)**b
    kminus12 = (kw+kw)/2
    kplus12 = (kw+kw)/2
    
    w = dt/(rho*cpw*ri*(dr**2))
    z = dt/(rho*cpw)
    
    A[k,k-1] = -w*rminus12*kminus12
    A[k,k] = 1 + w*rminus12*kminus12 + w*rplus12*kplus12
    A[k,k+1] = -w*rplus12*kplus12
    C[k,0] = Ti

ww = dt/(rho*cpw)
krminus12 = (kw+kw)/2

A[m-1,m-2] = -(2*ww/(dr**2))*krminus12
A[m-1,m-1] = 1 + (2*ww/(dr**2))*krminus12 + ww*((2/dr)+(b/r))*h
C[m-1,0] = Ti + ww*((2/dr)+(b/r))*h*Tinf

# check: print [A] and [C] to console, best viewed at small nr like nr=3
print('A \n', A)
print('C \n', C)

# solve system of equations [A]{T} = {C} for column vector {T}
for i in range(1, nt+1):
    T = np.linalg.solve(A,C)
    C = T.copy()
    C[m-1, 0] = T[m-1, 0] + ww*((2/dr)+(b/r))*h*Tinf
    TT = np.vstack((TT, T.T))

# check: display final T, best if nr = 3
print('T \n', T)

# plot results
py.figure(1)
py.plot(t, TT[:, m-1], '-k')
py.plot(t, TT[:, 0], '--k')
py.axhline(Tinf, color='r', linestyle='-.')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.ylim([Ti-20, Tinf+20])
py.xlim([0, tmax])
py.grid()
py.show()