import numpy as np
n_test = [2,4,8,16,32,64]
e_global = np.zeros((len(n_test)),dtype=np.double)

# test cube L x L x L 
L=1.0

nn=0 # index for global error
for n in n_test:
    x = np.zeros((n+1, n+1, n+1, 3), dtype=np.double)
    y = np.zeros((n+1, n+1, n+1), dtype=np.double)
    dx = L/n
    for i in range(0, n+1):
        for j in range(0, n+1):
            for k in range(0, n+1):
                x[i,j,k] = np.array([i*dx, j*dx, k*dx])
                y[i,j,k] = np.sin(2.0*np.pi * x[i,j,k,0] / ( L )) + np.sin(2.0*np.pi * x[i,j,k,1] / ( L ))  + np.sin(2.0*np.pi * x[i,j,k,2] / ( L ))
                #y_ref[i,j,k] = np.cos()
    
    # compute gradient
    gradientOfY = np.gradient(y,dx,dx,dx,edge_order=2)
    # build reference solution
    gradientOfYAnalytic =  np.zeros((n+1, n+1, n+1,3), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, n+1):
        for k in range(0, n+1):
          gradientOfYAnalytic[i,j,k] = np.array([ 2.0*np.pi / L * np.cos(2.0*np.pi * x[i,j,k,0] / ( L )),  2.0*np.pi / L * np.cos(2.0*np.pi * x[i,j,k,1] / ( L )), 2.0*np.pi / L * np.cos(2.0*np.pi * x[i,j,k,2] / ( L ))])

    # compute error
    e_local = np.zeros((n+1, n+1, n+1), dtype=np.double)
    for i in range(0, n+1):
      for j in range(0, n+1):
        for k in range(0, n+1):
          e_local[i,j,k] = np.sqrt( (gradientOfY[0][i,j,k] - gradientOfYAnalytic[i,j,k,0]) ** 2 +
                                    (gradientOfY[1][i,j,k] - gradientOfYAnalytic[i,j,k,1]) ** 2 +
                                    (gradientOfY[2][i,j,k] - gradientOfYAnalytic[i,j,k,2]) ** 2 )
          e_global[nn] = e_global[nn] +  e_local[i,j,k]
          
    N = (n+1) * (n+1) * (n+1)
    e_global[nn] = e_global[nn] / N
    nn = nn+1

import matplotlib.pyplot as plt
plt.plot(n_test, e_global, label='error')
#plt.plot(x, y, label='sin(x)')
#plt.plot(x, dy_dx, label="derivative of sin(x)")
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Differentiation of Sinusoidal Function on 3D grid')
plt.show()
    

  