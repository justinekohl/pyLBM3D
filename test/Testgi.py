import numpy as np


'''
    v: mass force field, m x m x m array with 1x3 vector
    g: kinematic viscosity (external force), m x m x m array with 1x3 vector
    gi: forcing term projected over the velocity space

'''
a = np.array([3,3,3])
b = np.array([2,2,2])

m = 2

a_array = np.zeros((m,m,m,3))
b_array = np.zeros((m,m,m,3))

for i in range(0, len(a_array)):
        for j in range(0, len(a_array[0])):
            for k in range(0, len(a_array[0][0])):
                a_array[i,j,k] = 3
                
for i in range(0, len(a_array)):
        for j in range(0, len(a_array[0])):
            for k in range(0, len(a_array[0][0])):
                b_array[i,j,k] = 2

ein_1 = np.einsum('ijkl,ijkl', a_array, b_array)
ein_2 = np.einsum('i,i', a, b)

print(ein_1)    #((3*2)*3)*8
print(ein_2)

# g = QS.g(rho, divSigma)
# v = QS.v(rho,j)
# gi = QS.gi(g,cc,w,rho,cs,v) # TODO this is cs = 1/sqrt(3)