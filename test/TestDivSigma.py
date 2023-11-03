import numpy as np
import QS.QuasiStatic as QS
import QS.ComputeDiv as div

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

def sigmaAndDivSigmaAnalytic(x,L, sigmaFunction, divSigmaFunction):
    n = len(x)-1 # 
    kk = len(x[0]) -1
    l = len(x[0][0])-1
    sigma = np.zeros((n+1, kk+1, l+1,3,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
          for k in range(0, l+1):
             sigma[i,j,k] = sigmaFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2]) 
    divOfSigmaAnalytic =  np.zeros((n+1, n+1, n+1,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, kk+1):
        for k in range(0, l+1):
          divOfSigmaAnalytic[i,j,k] = divSigmaFunction(x[i,j,k,0], x[i,j,k,1], x[i,j,k,2])
    return [sigma, divOfSigmaAnalytic]


def globalError(divOfSigma, divOfSigmaAnalytic):
  n = len(divOfSigma)-1
  kk = len(divOfSigma[0]) -1
  l = len(divOfSigma[0][0])-1
  e_local = np.zeros((n+1, kk+1, l+1), dtype=np.double)
  e_global = 0.0
  for i in range(0, n+1):
    for j in range(0, kk+1):
      for k in range(0, l+1):
        e_local[i,j,k] = np.sqrt( (divOfSigma[i,j,k,0] - divOfSigmaAnalytic[i,j,k,0]) ** 2 +
                                    (divOfSigma[i,j,k,1] - divOfSigmaAnalytic[i,j,k,1]) ** 2 +
                                    (divOfSigma[i,j,k,2] - divOfSigmaAnalytic[i,j,k,2]) ** 2  )
        e_global = e_global +  e_local[i,j,k]
  N = (n+1) * (kk+1) * (l+1)
  return e_global / ( N )


def sigmaCubic(x,y,z):
      a = 1.4
      b = 2.3
      c = -5.6
      sxx = a * x ** 3 
      sxy = a * x ** 3 + b * y ** 3 + c * z ** 3
      sxz = 0.0
      syy = b * y ** 2
      syz = c * z
      szz = x * z ** 3

      return np.array([ [sxx, sxy, sxz],
                       [sxy, syy, syz],
                       [sxz, syz, szz]])

def divSigmaCubic(x,y,z):
      a = 1.4
      b = 2.3
      c = -5.6
      dSxxDx = 3.0 * a * x ** 2
      dSxyDy = 3.0 * b * y ** 2
      dSxzDz = 0.0

      dSyxDx = 3.0 * a * x ** 2
      dSyyDy = 2.0 * b * y
      dSyzDz = c

      dSzxDx = 0.0
      dSzyDy = 0.0
      dSzzDz = 3.0 * x * z ** 2
      return np.array([ dSxxDx + dSxyDy + dSxzDz,
                        dSyxDx + dSyyDy + dSyzDz,
                        dSzxDx + dSzyDy + dSzzDz
         
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

    [sigma, divOfSigmaAnalytic] = sigmaAndDivSigmaAnalytic(x,L, sigmaCubic, divSigmaCubic)
    #divOfSigma= QS.divOfSigma(sigma,dx)                    ##### with gradient-function
    divOfSigma= div.ComupteDivergence(sigma,dx)   ##### with finite difference schemes
    e_global[nn] = globalError(divOfSigma, divOfSigmaAnalytic)
    nn = nn+1

import matplotlib.pyplot as plt
plt.plot(n_test, e_global, label='error')
#plt.plot(x, y, label='sin(x)')
#plt.plot(x, dy_dx, label="derivative of sin(x)")
#plt.yscale('log')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Computing divergence of Sigma on 3D Grid - Error')
plt.show()
    

  