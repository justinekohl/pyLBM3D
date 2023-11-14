import numpy as np
from findiff import Gradient, Divergence


def checkIfIndicesInArrayBounds(iArg, jArg, kArg, arrayArg):
    return iArg < len(arrayArg)  and iArg >= 0 and jArg < len(arrayArg[0]) and jArg >= 0 and kArg < len(arrayArg[0][0]) and kArg >= 0


def computeDivergenceUFromDisplacementField(uArg, dxArg):
    divUOut = np.zeros((len(uArg), len(uArg[0]), len(uArg[0][0])), dtype=np.double)
    uX = np.zeros((len(uArg), len(uArg[0]), len(uArg[0][0])), dtype=np.double)
    uY = np.zeros((len(uArg), len(uArg[0]), len(uArg[0][0])), dtype=np.double)
    uZ = np.zeros((len(uArg), len(uArg[0]), len(uArg[0][0])), dtype=np.double)


    for i in range(0, len(uArg)):
        for j in range(0, len(uArg[0])):
            for k in range(0, len(uArg[0][0])):
                uX[i][j][k] = uArg[i][j][k][0]
                uY[i][j][k] = uArg[i][j][k][1]
                uZ[i][j][k] = uArg[i][j][k][2]

    gradientUx = np.gradient(uX, dxArg, dxArg, dxArg, edge_order=2)
    gradientUy = np.gradient(uY, dxArg, dxArg, dxArg, edge_order=2)
    gradientUz = np.gradient(uZ, dxArg, dxArg, dxArg, edge_order=2)

    for i in range(0, len(uArg)):
        for j in range(0, len(uArg[0])):
            for k in range(0, len(uArg[0][0])):
                divUOut[i][j][k] = gradientUx[0][i][j][k] + gradientUy[1][i][j][k] + gradientUz[2][i][j][k]

    return divUOut

def computeGradientU(uArg, dxArg, *, use_np=True):
    uX = uArg[...,0]
    uY = uArg[...,1]
    uZ = uArg[...,2]
    nx, ny, nz = uX.shape
    gradU = np.empty((nx, ny, nz, 3, 3), dtype=np.double)

    if not use_np:
        grad = Gradient(h=[dxArg, dxArg, dxArg], acc=4)

        gradientUx = grad(uX)
        gradientUy = grad(uY)
        gradientUz = grad(uZ)
    else:
        gradientUx = np.gradient(uX, dxArg, dxArg, dxArg, edge_order=2)
        gradientUy = np.gradient(uY, dxArg, dxArg, dxArg, edge_order=2)
        gradientUz = np.gradient(uZ, dxArg, dxArg, dxArg, edge_order=2)

    # nx, ny, nz, *_ = uArg.shape
    gradu_tmp = np.array((gradientUx, gradientUy, gradientUz))
    for i in range(3):
        for j in range(3):
            gradU[...,i,j] = gradu_tmp[i,j,...]

    return gradU


def dJyDy(jArg, dxArg):
    j1d1 = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)
    j2d2 = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)
    j3d3 = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)
    uX = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)
    uY = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)
    uZ = np.zeros((len(jArg), len(jArg[0]), len(jArg[0][0])), dtype=np.double)


    for i in range(0, len(jArg)):
        for j in range(0, len(jArg[0])):
            for k in range(0, len(jArg[0][0])):
                uX[i][j][k] = jArg[i][j][k][0]
                uY[i][j][k] = jArg[i][j][k][1]
                uZ[i][j][k] = jArg[i][j][k][2]

    gradientUx = np.gradient(uX, dxArg, dxArg, dxArg, edge_order=2)
    gradientUy = np.gradient(uY, dxArg, dxArg, dxArg, edge_order=2)
    gradientUz = np.gradient(uZ, dxArg, dxArg, dxArg, edge_order=2)

    for i in range(0, len(jArg)):
        for j in range(0, len(jArg[0])):
            for k in range(0, len(jArg[0][0])):
                j1d1[i][j][k] = gradientUx[0][i][j][k] + gradientUy[1][i][j][k] + gradientUz[2][i][j][k]
                j2d2[i][j][k] = gradientUx[0][i][j][k] + gradientUy[1][i][j][k] + gradientUz[2][i][j][k]
                j3d3[i][j][k] = gradientUx[0][i][j][k] + gradientUy[1][i][j][k] + gradientUz[2][i][j][k]

    return [j1d1, j2d2, j3d3]

def trace(eps):
    trEps = np.zeros((len(eps), len(eps[0]), len(eps[0][0])), dtype=np.double)
    for i in range(0, len(eps)):
        for j in range(0, len(eps[0])):
            for k in range(0, len(eps[0][0])):
                trEps[i,j,k] = eps[i,j,k,0,0] +  eps[i,j,k,1,1] + eps[i,j,k,2,2]
    return trEps