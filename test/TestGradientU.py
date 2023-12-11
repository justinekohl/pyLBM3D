import numpy as np
import Util
import QS.ComputeGrad as Grad


def lattice(L,n):
  '''
  computes lattice for cube L x L x L with n x n x n segments
  '''
  dx = L/n
  x = np.zeros((n+1, n+1, n+1, 3), dtype=np.double)
  for i in range(0, n+1):
    for j in range(0, n+1):
      for k in range(0, n+1):
        x[i,j,k] = np.array([i*dx, j*dx, k*dx])
  return [dx, x]

def uAndGradientUAnalytic(x,L, uFunction, gradUFunction):
    n = len(x)-1 # 
    kk = len(x[0]) -1
    l = len(x[0][0])-1
    u = np.zeros((n+1, kk+1, l+1,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
          for k in range(0, l+1):
             u[i,j,k] = uFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2]) 
    gradientOfUAnalytic =  np.zeros((n+1, n+1, n+1,3,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
        for k in range(0, l+1):
          gradientOfUAnalytic[i,j,k] = gradUFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2])
    return [u, gradientOfUAnalytic]


def globalError(gradientOfU, gradientOfUAnalytic):
  n = len(gradientOfU)-1
  kk = len(gradientOfU[0]) -1
  l = len(gradientOfU[0][0])-1
  e_local = np.zeros((n+1, kk+1, l+1), dtype=np.double)
  e_global = 0.0
  for i in range(0, n+1):
    for j in range(0, kk+1):
      for k in range(0, l+1):
        e_local[i,j,k] = np.sqrt( (gradientOfU[i,j,k,0,0] - gradientOfUAnalytic[i,j,k,0,0]) ** 2 +
                                    (gradientOfU[i,j,k,0,1] - gradientOfUAnalytic[i,j,k,0,1]) ** 2 +
                                    (gradientOfU[i,j,k,0,2] - gradientOfUAnalytic[i,j,k,0,2]) ** 2 +
                                    (gradientOfU[i,j,k,1,0] - gradientOfUAnalytic[i,j,k,1,0]) ** 2 +
                                    (gradientOfU[i,j,k,1,1] - gradientOfUAnalytic[i,j,k,1,1]) ** 2 +
                                    (gradientOfU[i,j,k,1,2] - gradientOfUAnalytic[i,j,k,1,2]) ** 2 +
                                    (gradientOfU[i,j,k,2,0] - gradientOfUAnalytic[i,j,k,2,0]) ** 2 +
                                    (gradientOfU[i,j,k,2,1] - gradientOfUAnalytic[i,j,k,2,1]) ** 2 +
                                    (gradientOfU[i,j,k,2,2] - gradientOfUAnalytic[i,j,k,2,2]) ** 2  )
        e_global = e_global +  e_local[i,j,k]
  N = (n+1) * (kk+1) * (l+1)
  return e_global / ( N )


def uSin(x,y,z):
      ux = np.sin(2.0*np.pi * x /  L ) + np.sin(2.0*np.pi * y /  L )  + np.sin(2.0*np.pi * z /  L )
      uy = np.sin(2.0*np.pi * x /  L ) + np.sin(2.0*np.pi * y /  L )  + np.sin(2.0*np.pi * z /  L )
      uz = np.sin(2.0*np.pi * x /  L ) + np.sin(2.0*np.pi * y /  L )  + np.sin(2.0*np.pi * z /  L )
      return np.array([ux,uy,uz])

def gradUSin(x,y,z):
      duxDx = 2.0*np.pi / L * np.cos(2.0*np.pi * x / ( L ))
      duxDy = 2.0*np.pi / L * np.cos(2.0*np.pi * y/ ( L ))
      duxDz = 2.0*np.pi / L * np.cos(2.0*np.pi * z/ ( L ))
      duyDx = 2.0*np.pi / L * np.cos(2.0*np.pi * x / ( L ))
      duyDy = 2.0*np.pi / L * np.cos(2.0*np.pi * y / ( L ))
      duyDz = 2.0*np.pi / L * np.cos(2.0*np.pi * z / ( L ))
      duzDx = 2.0*np.pi / L * np.cos(2.0*np.pi * x / ( L ))
      duzDy = 2.0*np.pi / L * np.cos(2.0*np.pi * y / ( L ))
      duzDz = 2.0*np.pi / L * np.cos(2.0*np.pi * z / ( L ))
      return np.array([ [duxDx, duxDy, duxDz],
                       [duyDx, duyDy, duyDz],
                       [duzDx, duzDy, duzDz]
         
      ])

def uCubic(x,y,z):
      a = 1.4
      b = 2.3
      c = -5.6
      ux = a * x ** 3
      uy = a * x ** 2 + b * y ** 3
      uz = c * z ** 3
      return np.array([ux,uy,uz])

def gradUCubic(x,y,z):
      a = 1.4
      b = 2.3
      c = -5.6
      duxDx = 3.0 * a * x ** 2
      duxDy = 0.0
      duxDz = 0.0
      duyDx = 2.0 * a * x
      duyDy = 3.0 * b * y ** 2
      duyDz = 0.0
      duzDx = 0.0
      duzDy = 0.0
      duzDz = 3.0 * c * z ** 2
      return np.array([ [duxDx, duxDy, duxDz],
                       [duyDx, duyDy, duyDz],
                       [duzDx, duzDy, duzDz]
         
      ])

### Processing starts here
#    
n_test = [2,4,8,16,32,64]
e_global = np.zeros((len(n_test)),dtype=np.double)

# test cube L x L x L 
L=1.0

nn=0 # index for global error
for n in n_test:
    [dx, x] = lattice(L, n)

    [u, gradientOfUAnalytic] = uAndGradientUAnalytic(x,L, uCubic, gradUCubic)
    #gradientOfU = Util.computeGradientU(u,dxArg=dx)
    gradientOfU = Grad.ComupteGradientU(u,dx)
    e_global[nn] = globalError(gradientOfU, gradientOfUAnalytic)
    nn = nn+1

import matplotlib.pyplot as plt
plt.plot(n_test, e_global, label='error')
#plt.plot(x, y, label='sin(x)')
#plt.plot(x, dy_dx, label="derivative of sin(x)")
plt.legend()
plt.xlabel('nodes')
plt.ylabel('error')
plt.yscale('log')
plt.xscale('log',base=2)
plt.grid()
plt.title('Computing gradient of U on 3D Grid - Error')
plt.savefig('E:\git\PA_png\gradientU_test2.png')
plt.show()
    

  