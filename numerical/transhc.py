# function for transient heat conduction in a sphere

def sphere(k, rho, c, dt, dr, h, T, Tinf, i, j, lt, m, A, C):
    
    b = 2   # run model as a cylinder (b = 1) or as a sphere (b = 2)

    alpha = k/(rho*c)       # thermal diffusivity, alfa = kw / rho*cp, m^2/s
    Fo = alpha*dt/(dr**2)   # Fourier number, Fo = alfa*dt / dr^2, (-)
    Bi = h*dr/k             # Biot numbmer, Bi = h*dr / kw, (-)
    
    # center nodes T0 and T1
    A[0, 0] = 1 + 2*(1+b)*Fo[0]    # node T0
    A[0, 1] = -2*(1+b)*Fo[0]       # node T1
    C[0, 0] = T[i-1, 0]
    
    # internal nodes Tm-1, Tm, Tm+1
    A[j, j-1] = -Fo[j]*(1 - b/(2*(j+1)))   # Tm-1
    A[j, j] = 1 + 2*Fo[j]                  # Tm
    A[j, j+1] = -Fo[j]*(1 + b/(2*(j+1)))   # Tm+1
    C[j, 0] = T[i-1, j]
    
    # surface nodes Tr-1 and Tr
    A[m-1, m-2] = -2*Fo[m-1]                             # node Tr-1
    A[m-1, m-1] = 1 + 2*Fo[m-1]*(1 + Bi + (b/(2*m))*Bi)  # node Tr
    C[m-1, 0] = T[i-1, m-1] + 2*Fo[m-1]*Bi*(1 + b/(2*m))*Tinf
    
    # return array [A] and column vector {C}
    return A, C