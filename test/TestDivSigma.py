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

  return np.array([[sxx, sxy, sxz],
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
  return np.array([
    dSxxDx + dSxyDy + dSxzDz,
    dSyxDx + dSyyDy + dSyzDz,
    dSzxDx + dSzyDy + dSzzDz,
  ])

def sigma10(x,y,z):
  """sigma with max 10th order"""
  a = 1.4
  b = 2.3
  c = -5.6

  sxx = a * x ** 10
  sxy = a * x ** 10 + b * y ** 10 + c * z ** 10
  sxz = 0.0
  
  syy = b * y ** 8
  syz = c * z ** 5
  szz = x * z ** 3

  return np.array([
    [sxx, sxy, sxz],
    [sxy, syy, syz],
    [sxz, syz, szz],
  ])


def divSigma10(x,y,z):
  a = 1.4
  b = 2.3
  c = -5.6

  dSxxDx = 10 * a * x ** 9
  dSxyDy = 10 * b * y ** 9
  dSxzDz = 0.0

  dSyxDx = 10 * a * x ** 9
  dSyyDy = 8 * b * y ** 7
  dSyzDz = 5 * c * z ** 4

  dSzxDx = 0.0
  dSzyDy = 0.0
  dSzzDz = 3.0 * x * z ** 2
  return np.array([
    dSxxDx + dSxyDy + dSxzDz,
    dSyxDx + dSyyDy + dSyzDz,
    dSzxDx + dSzyDy + dSzzDz,
  ])


### Processing starts here
# n_test = [8,16,32]     # 4th order FD needs more than 5 nodes
n_test = np.logspace(3, 7, num=5, base=2, dtype=int)     # 4th order FD needs more than 5 nodes

# test cube L x L x L 
L = 1.0
# L = 2*np.pi

e_global_fd = []
e_global_np = []
e_global_div = []

divOfSigma_fd = None
divOfSigma_np = None
divOfSigma_div = None


for n in n_test:
  print(n, "--")
  [dx, x] = lattice(L, n)

  # [sigma, divOfSigmaAnalytic] = sigmaAndDivSigmaAnalytic(x,L, sigmaCubic, divSigmaCubic)
  [sigma, divOfSigmaAnalytic] = sigmaAndDivSigmaAnalytic(x, L, sigma10, divSigma10)
  
  divOfSigma_fd = QS.divOfSigma(sigma, dx)                    ##### with gradient-function
  e_global_fd.append(globalError(divOfSigma_fd, divOfSigmaAnalytic))
  
  divOfSigma_np = QS.divOfSigmaNumpy(sigma, dx)               ##### with findiff Divergence
  e_global_np.append(globalError(divOfSigma_np, divOfSigmaAnalytic))
  
  # divOfSigma_div = div.ComupteDivergence(sigma,dx)                ##### with finite difference schemes
  # e_global_div.append(globalError(divOfSigma_div, divOfSigmaAnalytic))



import matplotlib.pyplot as plt
# plot global error
fig, ax = plt.subplots()
if divOfSigma_np is not None:
  ax.plot(n_test, e_global_np, label='error NumPy')
if divOfSigma_fd is not None:
  ax.plot(n_test, e_global_fd, label='error FinDiff')
# if divOfSigma_div is not None:
#   ax.plot(n_test, e_global_div, label='error FD scheme')

ax.set_xscale('log', base=2) 
ax.set_yscale('log') 
ax.legend()
ax.set_xlabel('# nodes')
ax.set_ylabel('error')
ax.set_title(r'Computing divergence of $\sigma$ on 3D Grid - Error')
ax.grid()

fig.show()
