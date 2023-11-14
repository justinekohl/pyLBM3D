import numpy as np
import Util


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
    u = np.empty((n+1, kk+1, l+1,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
          for k in range(0, l+1):
             u[i,j,k] = uFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2]) 
    gradientOfUAnalytic =  np.empty((n+1, n+1, n+1,3,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
        for k in range(0, l+1):
          gradientOfUAnalytic[i,j,k] = gradUFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2])
    return [u, gradientOfUAnalytic]


def globalError(gradientOfU, gradientOfUAnalytic):
  n = len(gradientOfU)-1
  kk = len(gradientOfU[0]) -1
  l = len(gradientOfU[0][0])-1
  e_global = 0.0

  print("global error ...", end='\r')
  for i in range(0, n+1):
    for j in range(0, kk+1):
      for k in range(0, l+1):
        e_local = np.sqrt( (gradientOfU[i,j,k,0,0] - gradientOfUAnalytic[i,j,k,0,0]) ** 2 +
                           (gradientOfU[i,j,k,0,1] - gradientOfUAnalytic[i,j,k,0,1]) ** 2 +
                           (gradientOfU[i,j,k,0,2] - gradientOfUAnalytic[i,j,k,0,2]) ** 2 +
                           (gradientOfU[i,j,k,1,0] - gradientOfUAnalytic[i,j,k,1,0]) ** 2 +
                           (gradientOfU[i,j,k,1,1] - gradientOfUAnalytic[i,j,k,1,1]) ** 2 +
                           (gradientOfU[i,j,k,1,2] - gradientOfUAnalytic[i,j,k,1,2]) ** 2 +
                           (gradientOfU[i,j,k,2,0] - gradientOfUAnalytic[i,j,k,2,0]) ** 2 +
                           (gradientOfU[i,j,k,2,1] - gradientOfUAnalytic[i,j,k,2,1]) ** 2 +
                           (gradientOfU[i,j,k,2,2] - gradientOfUAnalytic[i,j,k,2,2]) ** 2  )
        e_global += e_local
  print("global error done")

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
# n_test = [8,16,32,64,128]     # 4th order FD needs more than 5 nodes
n_test = np.logspace(3, 7, num=5, base=2, dtype=int)     # 4th order FD needs more than 5 nodes

# test cube L x L x L 
L=1.0

gradientOfU_np = None       # numpy
gradientOfU_fd = None       # findiff

e_global_np = []
e_global_fd = []
for n in n_test:
    print(n, "--")
    [dx, x] = lattice(L, n)

    print("analytical ...", end='\r')
    [u, gradientOfUAnalytic] = uAndGradientUAnalytic(x,L, uSin, gradUSin)
    # [u, gradientOfUAnalytic] = uAndGradientUAnalytic(x,L, uCubic, gradUCubic)
    print("analytical done")
    
    # compute 2nd order NumPy gradient
    print("numpy ...", end='\r')
    gradientOfU_np = Util.computeGradientU(u, dxArg=dx)
    print("numpy done")
    e_global_np.append(globalError(gradientOfU_np, gradientOfUAnalytic))
    
    # compute 4th order FinDiff gradient
    print("findiff ...", end='\r')
    gradientOfU_fd = Util.computeGradientU(u, dxArg=dx, use_np=False)
    print("findiff done")
    e_global_fd.append(globalError(gradientOfU_fd, gradientOfUAnalytic))

import matplotlib.pyplot as plt

# plot global error
fig, ax = plt.subplots()
if gradientOfU_np is not None:
    ax.plot(n_test, e_global_np, label='error NumPy')
if gradientOfU_fd is not None:
    ax.plot(n_test, e_global_fd, label='error FinDiff')
ax.set_xscale('log', base=2) 
ax.set_yscale('log') 
#plt.plot(x, y, label='sin(x)')
#plt.plot(x, dy_dx, label="derivative of sin(x)")
ax.legend()
ax.set_xlabel('# nodes')
ax.set_ylabel('error')
ax.set_title('Computing gradient of U on 3D Grid - Error')
ax.grid()

nh = n_test[2], n_test[-2]
eh = [None, None]
eh[0] = e_global_np[2]*0.9
eh[1] = eh[0] * 10**-2

ax.plot(nh, eh, ':', c='tab:blue')

fig.show()
