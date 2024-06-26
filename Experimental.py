import numpy as np
import copy
import Util 
import Settings

#dx = Settings.getLatticeInformation()[0]


def intitialize(rho0Arg, csArg, ccArg, wArg, mArg, nArg, oArg, lamArg, mueArg):
    '''
    returns the initial distribution functions (as equilibrium distributions), moments and displacement field
    :param rho0Arg: the initial density, a scalar
    :param csArg: wave speed of shear waves, a scalar
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param mArg: dimension of lattice in x direction
    :param nArg: dimension of lattice in y direction
    :param oArg: dimension of lattice in z direction
    :param lamArg: Lame parameter
    :param mueArg: Lame parameter
    :return: [fOut, j0, sigma0, u0]
    '''
    sigma0 = np.zeros((mArg, nArg, oArg, 3, 3), dtype=np.double)
    j0 = np.zeros((mArg, nArg, oArg, 3), dtype=np.double)
    u0 = np.zeros((mArg, nArg, oArg, 3), dtype=np.double)
    rho = np.zeros((mArg, nArg, oArg), dtype=np.double)
    rho.fill(rho0Arg)
    fOut = equilibriumDistribution(rho, j0, sigma0, ccArg, wArg, csArg, lamArg, mueArg, rho0Arg)
    return [fOut, j0, sigma0, u0]


def equilibriumDistribution(rhoArg, jArg, sigmaArg, ccArg, wArg, csArg, lamArg, mueArg, rho0Arg):
    '''
    :param rhoArg: the density in lattice dimensions (m,n,o)
    :param jArg: the momentum in lattice dimensions (m,n,o,3)
    :param sigmaArg: the stress in lattice dimensions (m,n,o,3,3)
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param csArg: wave speed of shear waves, a scalar
    :param lamArg: Lame parameter
    :param mueArg: Lame parameter
    :param rho0Arg: the initial density, a scalar
    :return: feqOut the equilibrium distribution in lattice dimensions (m,n,o)
    '''
    feqOut = np.zeros((len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]), 27), dtype=np.double)
    I = np.identity(3, dtype=np.double)
    II = np.zeros((3, 3, 3), dtype=np.double)
    II[0, 0, 0] = 1.0
    II[1, 1, 1] = 1.0
    II[2, 2, 2] = 1.0
    for i in range(0, len(feqOut)):
        for j in range(0, len(feqOut[0])):
            for k in range(0, len(feqOut[0][0])):
                for l in range(0, len(feqOut[0][0][0])):
                    tmp1 = np.tensordot(ccArg[l], jArg[i][j][k], axes=1)
                    tmp2 = np.tensordot(
                        (-sigmaArg[i][j][k] - rhoArg[i][j][k] * csArg ** 2 * I),
                        (np.outer(ccArg[l], ccArg[l].transpose()) - csArg ** 2 * I),
                        axes=2)

                    leftTmp3 = np.einsum('a,bc->abc', jArg[i, j, k], I) - np.einsum('c,abc->abc', jArg[i, j, k], II)

                    rightTmp3 = np.einsum('a,b,c->abc', ccArg[l], ccArg[l], ccArg[l]) - csArg ** 2 * (
                        np.einsum('a,bc->abc', ccArg[l], I) + np.einsum('b,ac->abc', ccArg[l], I) +
                        np.einsum('c,ab->abc', ccArg[l], I))

                    tmp3 = (np.einsum('abc,abc', leftTmp3, rightTmp3))
                    feqOut[i][j][k][l] = wArg[l] * (
                                rhoArg[i][j][k] + 1.0 / (csArg ** 2) * tmp1 + 1.0 / (2.0 * csArg ** 4) * tmp2 + 1.0 / (
                                    6.0 * csArg ** 6) * (lamArg - mueArg) / rho0Arg * tmp3)

    return feqOut


def sourceTermPsi(bArg, rho0Arg, dJyDy, ccArg, wArg, csArg, mueArg, laArg):
    '''
    :param bArg: the body force per unit mass in lattice dimensions (m,n,o,3)
    :param rho0Arg: the original density, a scalar
    :param dJyDy: an array of the derivatives of the momentum [dJxDx, dJyDy, dJzDz] in lattice dimensions (3,m,n,o)
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param csArg: wave speed of shear waves, a scalar
    :param mueArg: Lame parameter
    :param laArg: Lame parameter
    :return: the source term in the LBE in lattice dimensions (m,n,o,27)
    '''
    psiOut = np.zeros((len(bArg), len(bArg[0]), len(bArg[0][0]), 27), dtype=np.double)
    I = np.identity(3, dtype=np.double)
    for i in range(0, len(psiOut)):
        for j in range(0, len(psiOut[0])):
            for k in range(0, len(psiOut[0][0])):
                for l in range(0, len(psiOut[0][0][0])):
                    left = np.array([[dJyDy[0][i, j, k], 0, 0], [0, dJyDy[1][i, j, k], 0], [0, 0, dJyDy[2][i, j, k]]], dtype=float)
                    right = np.outer(ccArg[l], ccArg[l].transpose()) - csArg ** 2 * I
                    tmp = np.einsum('ab,ab', left, right)
                    psiOut[i][j][k][l] = wArg[l] * (rho0Arg * 1.0 / (csArg ** 2) * np.tensordot(ccArg[l], bArg[i][j][k], axes=1) + 1.0 / (csArg ** 4) * (mueArg - laArg) / rho0Arg * tmp)
    return psiOut

def firstSource(bArg, rho0Arg):
    '''
    :param bArg: the body force per unit mass in lattice dimensions (m,n,o,3)
    :param rho0Arg: the original density, a scalar
    :return: the source term in the first moment (m,n,o,3)
    '''
    SOut = np.zeros((len(bArg), len(bArg[0]), len(bArg[0][0]), 3), dtype=np.double)
    for i in range(0, len(SOut)):
        for j in range(0, len(SOut[0])):
            for k in range(0, len(SOut[0][0])):
                SOut[i,j,k] = bArg[i,j,k] * rho0Arg
    return SOut

def secondSource(dJyDy, laArg, mueArg, rho0Arg):
    '''
    :param dJyDy: an array of the derivatives of the momentum [dJxDx, dJyDy, dJzDz] in lattice dimensions (3,m,n,o)
    :param laArg: Lame parameter
    :param mueArg: Lame parameter
    :param rho0Arg: the original density, a scalar
    :return: the source term in the second moment (m,n,o,3,3)
    '''
    SOut = np.zeros((len(dJyDy[0]), len(dJyDy[0][0]), len(dJyDy[0][0][0]), 3, 3), dtype=np.double)
    #I = np.identity(3, dtype=np.double)
    for i in range(0, len(SOut)):
        for j in range(0, len(SOut[0])):
            for k in range(0, len(SOut[0][0])):
                SOut[i, j, k] = (mueArg - laArg) / rho0Arg * np.array([[dJyDy[0][i,j,k], 0, 0], [0, dJyDy[1][i,j,k], 0], [0, 0, dJyDy[2][i,j,k]]], dtype=float)
    return SOut


def rho(zerothMomentArg):
    '''
    :param zerothMomentArg: the zeroth moment in lattice dimensions (m,n,o)
    :return: the current density in lattice dimensions (m,n,o)
    '''
    return copy.deepcopy(zerothMomentArg)


def j(firstMomentArg,dtArg,firstSourceArg):
    '''
    :param firstMomentArg: the first moment in lattice dimensions (m,n,o,3)
    :param dtArg: the time step size, a scalar
    :param firstSourceArg: the first order source term in lattice dimensions (m,n,o,3)
    :return: the momentum in lattice dimensions (m,n,o,3)
    '''
    jOut = np.zeros(firstMomentArg.shape, dtype=np.double)
    for i in range(0, len(jOut)):
        for j in range(0, len(jOut[0])):
            for k in range(0, len(jOut[0][0])):
                jOut[i,j,k] = firstMomentArg[i,j,k] + dtArg/2.0 * firstSourceArg[i,j,k]
    return jOut


def sigma(secondMomentArg,dtArg,secondSourceArg):
    '''
    :param secondMomentArg: the second moment in lattice dimensions (m,n,o,3,3)
    :param dtArg: the time step size, a scalar
    :param secondSourceArg: the second order source term in lattice dimensions (m,n,o,3,3)
    :return: the Cauchy stress in lattice dimensions (m,n,o,3,3)
    '''
    sigmaOut = np.zeros(secondMomentArg.shape, dtype=np.double)
    for i in range(0, len(sigmaOut)):
        for j in range(0, len(sigmaOut[0])):
            for k in range(0, len(sigmaOut[0][0])):
                sigmaOut[i, j, k] = -secondMomentArg[i, j, k] + dtArg/2.0 * secondSourceArg[i, j, k]
    return sigmaOut


def computeU(uOldArg, rho0Arg, jArg, jOldArg, dtArg):
    '''
    :param uOldArg: the displacement from the last time step in lattice dimensions (m,n,o,3)
    :param rho0Arg: the original density, a scalar
    :param jArg: the current momentum in lattice dimensions (m,n,o,3)
    :param jOldArg: the momentum of the last time step in lattice dimensions (m,n,o,3)
    :param dtArg: the time step size
    :return: the displacement for the current time step in lattice dimensions (m,n,o,3)
    '''
    uNew = np.zeros((len(uOldArg), len(uOldArg[0]), len(uOldArg[0][0]), 3), dtype=np.double)
    for i in range(0, len(uOldArg)):  # f0
        for j in range(0, len(uOldArg[0])):
            for k in range(0, len(uOldArg[0][0])):
                uNew[i][j][k] = uOldArg[i][j][k] + (jArg[i][j][k] + jOldArg[i][j][k]) / rho0Arg / 2.0 * dtArg
    return uNew

def computeJForU(uDesiredArg,uOldArg,dtArg,rho0Arg,jOldArg):
    '''
    Solves the integration rule for j if a desired u is given
    :param uDesiredArg:
    :param uOldArg:
    :param dtArg:
    :param rho0Arg:
    :param jOldArg:
    :return:
    '''
    j = (uDesiredArg-uOldArg) / dtArg * rho0Arg * 2.0 -jOldArg
    return j


import BoundaryConditions as BC

def dirichletBoundaryRule(jBdArg, csArg, ccArg, wArg, fCollRelevantArg):
    fBouncedBack = - fCollRelevantArg + 2.0 / (csArg ** 2) * wArg * np.tensordot(ccArg,jBdArg,axes=1)
    return fBouncedBack

def neumannBoundaryRule(sigmaBdArg, rhoBdArg, csArg, ccArg, wArg, fCollRelevantArg):
    tmp1 = -sigmaBdArg - rhoBdArg * csArg ** 2 * np.identity(3, dtype=np.double)
    tmp2 = np.outer(ccArg, ccArg.transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
    fBouncedBack = - fCollRelevantArg + 2.0 * wArg * (rhoBdArg + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1, tmp2, axes=2)))
    return fBouncedBack

def bounceBackFsAtPlane(fArg, fCollArg, coordinateArg, coordinateValueArg):
    fCollRelevant = BC.selectAtCoordinate(fCollArg, coordinateArg, coordinateValueArg)
    fRelevant = BC.selectAtCoordinate(fArg, coordinateArg, coordinateValueArg)
    indicesMissing = BC.getMissingDistributionFunctionIndices(fArg, coordinateArg, coordinateValueArg)
    return [fCollRelevant, fRelevant, indicesMissing]


def identity(sigmaBdArg, latticePointLocationArg):
    '''
        A default transform function for sigma 
    '''
    return sigmaBdArg

def applyNeumannBoundaryConditions(fArg, fCollArg, rhoArg, csArg, ccArg, wArg, sigmaBdArg, sigmaArg, coordinateArg='x', coordinateValueArg=0, boundaryRule = neumannBoundaryRule, sigmaTransformFunction = identity, dx=0):
    '''
    :param fArg: the distribution function before the boundary conditions have been applied at the given plane in lattice dimensions (m,n,o)
    :param fCollArg: the distribution function after collision has been applied in lattice dimensions (m,n,o)
    :param rhoArg: the density in lattice dimensions (m,n,o)
    :param csArg: csArg: wave speed of shear waves, a scalar
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param sigmaBdArg: the prescribed stress at this plane with nan's for undefined values (3,3)
    :param sigmaArg: the stress field in lattice dimensions (m,n,o,3,3)
    :param coordinateArg: 'x', 'y', 'z' the coordinate direction identifying the plane
    :param coordinateValueArg: the index identifying the plane (either 0 or max in the respective direction)
    :return: the distribution function after the boundary conditions have been applied at the given plane in lattice dimensions (m,n,o)
    '''
    def computeFieldsForBounceBack(rhoArg, ccArg, coordinateArg, coordinateValueArg):
        rhoBd = BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg, coordinateValueArg)
        sigmaBd = computeSigmaBd(sigmaBdArg, sigmaArg, ccArg, coordinateArg, coordinateValueArg, sigmaBdTransformFunction=sigmaTransformFunction,dx=dx)
        return [rhoBd, sigmaBd]



    #rhoBd = BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg, coordinateValueArg)
    #fCollRelevant = BC.selectAtCoordinate(fCollArg, coordinateArg, coordinateValueArg)
    fOut = copy.deepcopy(fArg)
    #fRelevant = BC.selectAtCoordinate(fOut, coordinateArg, coordinateValueArg)
    #indicesMissing = BC.getMissingDistributionFunctionIndices(fArg, coordinateArg, coordinateValueArg)

    [fCollRelevant, fRelevant, indicesMissing] = bounceBackFsAtPlane(fOut, fCollArg, coordinateArg, coordinateValueArg)

    #sigmaBd = computeSigmaBd(sigmaBdArg, sigmaArg, ccArg, coordinateArg, coordinateValueArg)

    [rhoBd, sigmaBd] = computeFieldsForBounceBack(rhoArg, ccArg, coordinateArg, coordinateValueArg)

    for i in range(0, len(fRelevant)):
        for j in range(0, len(fRelevant[0])):
            for l in indicesMissing:
                oL = BC.getOppositeLatticeDirection(l)
                #[i1, i2, i3] = getCurrentIndicesInXYZOrder(coordinateArg, coordinateValueArg, i, j)
                #oL = BC.getOppositeLatticeDirection(
                #    BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]), i1,
                #                                        i2, i3))

                #tmp1 = -sigmaBd[i, j, l] - rhoBd[i, j, l] * csArg ** 2 * np.identity(3, dtype=np.double)
                #tmp2 = np.outer(ccArg[oL], ccArg[oL].transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
                #fRelevant[i, j, l] = - fCollRelevant[i, j, oL] + 2.0 * wArg[oL] * (rhoBd[i, j, l] + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1, tmp2, axes=2)))
                fRelevant[i, j, l] = boundaryRule(sigmaBd[i, j, l], rhoBd[i, j, l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, j, oL])

    return fOut


def applyNeumannBoundaryConditionsAtEdge(fArg, fCollArg,  rhoArg,  csArg, ccArg,  wArg, sigmaBdArg1, sigmaBdArg2, sigmaArg, coordinateArg1='x', coordinateValueArg1=0, coordinateArg2='y', coordinateValueArg2=0, boundaryRule=neumannBoundaryRule, sigmaTransformFunction = identity, dx=0):
    '''

    :param fArg: the distribution function before the boundary conditions have been applied at the given plane in lattice dimensions (m,n,o)
    :param fCollArg: the distribution function after collision has been applied in lattice dimensions (m,n,o)
    :param rhoArg: the density in lattice dimensions (m,n,o)
    :param csArg: csArg: wave speed of shear waves, a scalar
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param sigmaBdArg1: the prescribed stress at this plane with nan's for undefined values (3,3) at the first plane
    :param sigmaBdArg2: the prescribed stress at this plane with nan's for undefined values (3,3) at the second plane
    :param sigmaArg: the stress field in lattice dimensions (m,n,o,3,3)
    :param coordinateArg1: 'x', 'y', 'z' the coordinate direction identifying the plane 1
    :param coordinateValueArg1: the index identifying the plane 1 (either 0 or max in the respective direction)
    :param coordinateArg2: 'x', 'y', 'z' the coordinate direction identifying the plane 2
    :param coordinateValueArg2: the index identifying the plane 2 (either 0 or max in the respective direction)
    :return: the distribution function after the boundary conditions have been applied at the given edge in lattice dimensions (m,n,o)
    '''
    rhoBd = BC.reduceSurfaceToEdge(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg1, coordinateValueArg1),coordinateArg1,coordinateArg2,coordinateValueArg2)
    fCollRelevant = BC.selectAtEdge(fCollArg,coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    fOut = copy.deepcopy(fArg)
    fRelevant = BC.selectAtEdge(fOut, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    indicesMissing = BC.getMissingDistributionFunctionIndicesAtEdge(fArg, coordinateArg1,coordinateValueArg1, coordinateArg2, coordinateValueArg2)

    sigmaBd = 1.0/2.0 * (sigmaBdArg1 + sigmaBdArg2)
    sigmaBd = BC.reduceSurfaceToEdge(computeSigmaBd(sigmaBd, sigmaArg, ccArg, coordinateArg1, coordinateValueArg1, sigmaBdTransformFunction=sigmaTransformFunction,dx=dx), coordinateArg1, coordinateArg2,coordinateValueArg2)

    for i in range(0, len(fRelevant)):
        for l in indicesMissing:
            oL = BC.getOppositeLatticeDirection(l)
            #[i1, i2, i3] = getCurrentIndicesAtEdgeInXYZOrder(coordinateArg1, coordinateValueArg1, coordinateArg2,
            #                                                 coordinateValueArg2, i)
            #oL = BC.getOppositeLatticeDirection(
            #    BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]), i1, i2,
            #                                        i3))

            #tmp1 = -sigmaBd[i, l] - rhoBd[i, l] * csArg ** 2 * np.identity(3,dtype=np.double)
            #tmp2 = np.outer(ccArg[oL], ccArg[oL].transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
            #fRelevant[i,l] = - fCollRelevant[i,  oL] + 2.0 * wArg[oL] * (rhoBd[i,l] + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1,tmp2, axes=2)))

            fRelevant[i, l] = boundaryRule(sigmaBd[i,  l], rhoBd[i,  l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, oL])
    return fOut


def applyNeumannBoundaryConditionsAtCorner(fArg, fCollArg, rhoArg,  csArg, ccArg,  wArg, sigmaBdArg1, sigmaBdArg2, sigmaBdArg3, sigmaArg,  coordinateValueArg1=0, coordinateValueArg2=0, coordinateValueArg3=0, boundaryRule=neumannBoundaryRule, sigmaTransformFunction = identity, dx=0):
    '''

    :param fArg: the distribution function before the boundary conditions have been applied at the given plane in lattice dimensions (m,n,o)
    :param fCollArg: the distribution function after collision has been applied in lattice dimensions (m,n,o)
    :param rhoArg: the density in lattice dimensions (m,n,o)
    :param csArg: csArg: wave speed of shear waves, a scalar
    :param ccArg: array of lattice speeds, (27,3)
    :param wArg: weights of lattice directions (27)
    :param sigmaBdArg1: the prescribed stress at this plane with nan's for undefined values (3,3) at the first plane, x-direction
    :param sigmaBdArg2: the prescribed stress at this plane with nan's for undefined values (3,3) at the second plane, y-direction
    :param sigmaBdArg3: the prescribed stress at this plane with nan's for undefined values (3,3) at the third plane, z-direction
    :param sigmaArg: the stress field in lattice dimensions (m,n,o,3,3)
    :param coordinateValueArg1: the index identifying the plane 1 (either 0 or max in the x direction)
    :param coordinateValueArg2: the index identifying the plane 2 (either 0 or max in the y direction)
    :param coordinateValueArg3: the index identifying the plane 3 (either 0 or max in the z direction)
    :return: the distribution function after the boundary conditions have been applied at the given corner in lattice dimensions (m,n,o)
    '''
    rhoBd = BC.reduceSurfaceToCorner(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, 'x', coordinateValueArg1),coordinateValueArg2, coordinateValueArg3)
    fCollRelevant = fCollArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    fOut = copy.deepcopy(fArg)
    fRelevant = fOut[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    indicesMissing = BC.getMissingDistributionFunctionIndicesAtCorner(fArg, coordinateValueArg1, coordinateValueArg2, coordinateValueArg3)

    sigmaBd = 1.0/3.0 * (sigmaBdArg1 + sigmaBdArg2 + sigmaBdArg3) # TODO average here okay?
    sigmaBd = BC.reduceSurfaceToCorner(
        computeSigmaBd(sigmaBd, sigmaArg, ccArg, 'x', coordinateValueArg1, sigmaBdTransformFunction=sigmaTransformFunction,dx=dx), coordinateValueArg2, coordinateValueArg3)

    for l in indicesMissing:
        oL = BC.getOppositeLatticeDirection(l)
        #oL = BC.getOppositeLatticeDirection(
        #    BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]),
        #                                        coordinateValueArg1, coordinateValueArg2, coordinateValueArg3))

        #tmp1 = -sigmaBd[l] - rhoBd[l] * csArg ** 2 * np.identity(3,dtype=np.double)
        #tmp2 = np.outer(ccArg[oL], ccArg[oL].transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
        #fRelevant[l] = - fCollRelevant[oL] + 2.0 * wArg[oL] * (rhoBd[l] + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1,tmp2, axes=2)))
        fRelevant[l] = boundaryRule(sigmaBd[l], rhoBd[l], csArg, ccArg[oL], wArg[oL], fCollRelevant[oL])

    return fOut


def computeSigmaBd(sigmaBC, sigmaArg, ccArg, coordinateArg='x', coordinateValueArg=0, sigmaBdTransformFunction=identity, dx=0):
    '''

    :param sigmaBC: the prescribed stress at this plane with nan's for undefined values (3,3)
    :param sigmaBdTransformFunction: a function that allows to compute location dependent stress values
    :param sigmaArg:  the stress field in lattice dimensions (m,n,o,3,3)
    :param ccArg: array of lattice speeds, (27,3)
    :param coordinateArg: 'x', 'y', 'z' the coordinate direction identifying the plane
    :param coordinateValueArg: the index identifying the plane (either 0 or max in the respective direction)
    :return: the stress field for all lattice links at all lattice points of a given plane, accounting for boundary conditions; in plane dimensions (k,l,27,3,3)
    '''
    def applyTractionBoundaryConditions(sigmaBd, sigmaBC):
        '''

        :param sigmaBd: the current stress field at the plane in plane dimensions (k,l,27,3,3)
        :param sigmaBC: the prescribed stress at this plane with nan's for undefined values (3,3)
        :return: the stress field for all lattice links at all lattice points of a given plane, accounting for boundary conditions; in plane dimensions (k,l,27,3,3)
        '''
        sigmaBd = copy.deepcopy(sigmaBd)
        for i in range(0, len(sigmaBd)):
            for j in range(0, len(sigmaBd[i])):
                ## this location can be used to compute stress location dependent (similiarly for edges and corners)
                latticePointLocation = Util.getLocationSurface(coordinateArg=coordinateArg, coordinateValueArg=coordinateValueArg,i=i, j=j, dx=dx)
                sigmaBCTransformed = sigmaBdTransformFunction(sigmaBC, latticePointLocation)
                for l in range(0, len(sigmaBd[i][j])):
                    for ii in range(0, len(sigmaBd[i][j][l])):
                        for jj in range(0, len(sigmaBd[i][j][l][ii])):
                            if (not np.isnan(sigmaBC[ii,jj])):    
                                sigmaBd[i,j,l,ii,jj] = sigmaBCTransformed[ii,jj]
        return sigmaBd

    def computeSigmaBdWithoutExtrapolation(sigmaArg, ccArg, coordinateArg, coordinateValueArg):
        '''

        :param sigmaArg: the stress field in lattice dimensions (m,n,o,3,3)
        :param ccArg: array of lattice speeds, (27,3)
        :param coordinateArg: 'x', 'y', 'z' the coordinate direction identifying the plane
        :param coordinateValueArg: the index identifying the plane (either 0 or max in the respective direction)
        :return: the current stress field at the plane in plane dimensions (k,l,27,3,3)
        '''
        sigmaAtCoordinate = BC.selectAtCoordinate(sigmaArg,coordinateArg,coordinateValueArg)
        sigmaBd = np.zeros((len(sigmaAtCoordinate),len(sigmaAtCoordinate[0]),len(ccArg),3,3), dtype=np.double)

        for i in range(0, len(sigmaBd)):
            for j in range(0, len(sigmaBd[i])):
                for l in range(0, len(sigmaBd[i][j])):
                    sigmaBd[i,j,l] = sigmaAtCoordinate[i,j]
        return sigmaBd

    sigmaBd = applyTractionBoundaryConditions(computeSigmaBdWithoutExtrapolation(sigmaArg,ccArg,coordinateArg,coordinateValueArg), sigmaBC)
    return sigmaBd


def applyDirichletBoundaryConditions(fArg, fCollArg, rho0Arg, csArg, ccArg, wArg, uBdArg, uOldArg, jOldArg, dtArg, coordinateArg='x', coordinateValueArg=0, boundaryRule = dirichletBoundaryRule):
    def computeFieldsForBounceBack(rho0Arg, uOldArg, uBdArg, jOldArg, dtArg, ccArg, coordinateArg, coordinateValueArg):
        uOldBd = BC.selectAtCoordinate(uOldArg,coordinateArg,coordinateValueArg)
        jOldBd = BC.selectAtCoordinate(jOldArg,coordinateArg,coordinateValueArg)
        jBd = np.zeros((len(jOldBd), len(jOldBd[0]),len(ccArg),3))
        for i in range(0,len(jOldBd)):
            for j in range(0,len(jOldBd[i])):
                for l in range(0,len(ccArg)):
                    jBd[i,j,l] = computeJForU(uBdArg, uOldBd[i,j], dtArg, rho0Arg, jOldBd[i,j])
                    #jBd[i, j, l] = np.array([0.0, 0.0, 0.0])
        return jBd

    fOut = copy.deepcopy(fArg)

    [fCollRelevant, fRelevant, indicesMissing] = bounceBackFsAtPlane(fOut, fCollArg, coordinateArg, coordinateValueArg)

    jBd = computeFieldsForBounceBack(rho0Arg, uOldArg, uBdArg, jOldArg, dtArg, ccArg, coordinateArg, coordinateValueArg)

    for i in range(0, len(fRelevant)):
        for j in range(0, len(fRelevant[0])):
            for l in indicesMissing:
                oL = BC.getOppositeLatticeDirection(l)
                fRelevant[i, j, l] = boundaryRule(jBd[i, j, l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, j, oL])

    return fOut


def applyDirichletNeumannBoundaryConditionsAtEdge(fArg, fCollArg,  rhoArg,  csArg, ccArg,  wArg, uBdArg, rho0Arg, uOldArg, jOldArg, dtArg, sigmaBdArg, sigmaArg,  coordinateArg1='x', coordinateValueArg1=0, coordinateArg2='y', coordinateValueArg2=0):
    def computeFieldsForBounceBack():
        rhoBd = BC.reduceSurfaceToEdge(
            BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg1, coordinateValueArg1), coordinateArg1,
            coordinateArg2, coordinateValueArg2)
        sigmaBd = sigmaBdArg
        sigmaBd = BC.reduceSurfaceToEdge(computeSigmaBd(sigmaBd, sigmaArg, ccArg, coordinateArg1, coordinateValueArg1),
                                         coordinateArg1, coordinateArg2, coordinateValueArg2)

        uOldBd = BC.reduceSurfaceToEdge(BC.selectAtCoordinate(uOldArg,coordinateArg1,coordinateValueArg1),coordinateArg1, coordinateArg2, coordinateValueArg2)
        jOldBd = BC.reduceSurfaceToEdge(BC.selectAtCoordinate(jOldArg,coordinateArg1,coordinateValueArg1),coordinateArg1, coordinateArg2, coordinateValueArg2)
        jBd = np.zeros((len(jOldBd),len(ccArg),3))
        for i in range(0,len(jOldBd)):
            for l in range(0,len(ccArg)):
                jBd[i,l] = computeJForU(uBdArg, uOldBd[i], dtArg, rho0Arg, jOldBd[i])
                    #jBd[i, j, l] = np.array([0.0, 0.0, 0.0])
        return [jBd, rhoBd, sigmaBd]

    def computeFsForBounceBack():
        fCollRelevant = BC.selectAtEdge(fCollArg, coordinateArg1, coordinateValueArg1, coordinateArg2,
                                        coordinateValueArg2)
        fOut = copy.deepcopy(fArg)
        fRelevant = BC.selectAtEdge(fOut, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
        indicesMissing = BC.getMissingDistributionFunctionIndicesAtEdge(fArg, coordinateArg1, coordinateValueArg1,
                                                                        coordinateArg2, coordinateValueArg2)
        return [fCollRelevant, fOut, fRelevant, indicesMissing]


    #rhoBd = BC.reduceSurfaceToEdge(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg1, coordinateValueArg1),coordinateArg1,coordinateArg2,coordinateValueArg2)
    #fCollRelevant = BC.selectAtEdge(fCollArg,coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    #fOut = copy.deepcopy(fArg)
    #fRelevant = BC.selectAtEdge(fOut, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    #indicesMissing = BC.getMissingDistributionFunctionIndicesAtEdge(fArg, coordinateArg1,coordinateValueArg1, coordinateArg2, coordinateValueArg2)

    #sigmaBd = sigmaBdArg
    #sigmaBd = BC.reduceSurfaceToEdge(computeSigmaBd(sigmaBd, sigmaArg, ccArg, coordinateArg1, coordinateValueArg1), coordinateArg1, coordinateArg2,coordinateValueArg2)

    [jBd, rhoBd, sigmaBd] = computeFieldsForBounceBack()
    [fCollRelevant, fOut, fRelevant, indicesMissing] = computeFsForBounceBack()

    for i in range(0, len(fRelevant)):
        for l in indicesMissing:
            oL = BC.getOppositeLatticeDirection(l)

            #tmp1 = -sigmaBd[i, l] - rhoBd[i, l] * csArg ** 2 * np.identity(3,dtype=np.double)
            #tmp2 = np.outer(ccArg[oL], ccArg[oL].transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
            #fRelevant[i,l] = - fCollRelevant[i,  oL] + 2.0 * wArg[oL] * (rhoBd[i,l] + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1,tmp2, axes=2)))
            #fRelevant[i, l] = (neumannBoundaryRule(sigmaBd[i,  l], rhoBd[i,  l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, oL]) + dirichletBoundaryRule(jBd[i, l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, oL])) / 2.0
            fRelevant[i, l] = dirichletBoundaryRule(jBd[i, l], csArg,
                                                    ccArg[oL], wArg[oL],
                                                    fCollRelevant[
                                                        i, oL])
    return fOut


def applyDirichletNeumannBoundaryConditionsAtCorner(fArg, fCollArg, rhoArg,  csArg, ccArg,  wArg, uBdArg, rho0Arg, uOldArg, jOldArg, dtArg, sigmaBdArg2, sigmaBdArg3, sigmaArg,  coordinateValueArg1=0, coordinateValueArg2=0, coordinateValueArg3=0):
    def computeFieldsForBounceBack():
        rhoBd = BC.reduceSurfaceToCorner(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, 'x', coordinateValueArg1),
                                         coordinateValueArg2, coordinateValueArg3)
        sigmaBd2 = BC.reduceSurfaceToCorner(
        computeSigmaBd(sigmaBdArg2, sigmaArg, ccArg, 'x', coordinateValueArg1), coordinateValueArg2, coordinateValueArg3)
        sigmaBd3 = BC.reduceSurfaceToCorner(
            computeSigmaBd(sigmaBdArg3, sigmaArg, ccArg, 'x', coordinateValueArg1), coordinateValueArg2,
            coordinateValueArg3)

        uOldBd = uOldArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
        jOldBd = jOldArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
        jBd = np.zeros((len(ccArg), 3))
        for l in range(0, len(ccArg)):
                jBd[l] = computeJForU(uBdArg, uOldBd, dtArg, rho0Arg, jOldBd)
                # jBd[i, j, l] = np.array([0.0, 0.0, 0.0])
        return [rhoBd, sigmaBd2, sigmaBd3, jBd]

    def computeFsForBounceBack():
        fCollRelevant = fCollArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
        fOut = copy.deepcopy(fArg)
        fRelevant = fOut[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
        indicesMissing = BC.getMissingDistributionFunctionIndicesAtCorner(fArg, coordinateValueArg1,
                                                                          coordinateValueArg2, coordinateValueArg3)
        return [fCollRelevant, fOut, fRelevant, indicesMissing]


    #rhoBd = BC.reduceSurfaceToCorner(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, 'x', coordinateValueArg1),coordinateValueArg2, coordinateValueArg3)
    #fCollRelevant = fCollArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    #fOut = copy.deepcopy(fArg)
    #fRelevant = fOut[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    #indicesMissing = BC.getMissingDistributionFunctionIndicesAtCorner(fArg, coordinateValueArg1, coordinateValueArg2, coordinateValueArg3)

    [rhoBd, sigmaBd2, sigmaBd3, jBd] = computeFieldsForBounceBack()

    [fCollRelevant, fOut, fRelevant, indicesMissing] = computeFsForBounceBack()

    #sigmaBd = 1.0/3.0 * (sigmaBdArg1 + sigmaBdArg2 + sigmaBdArg3) # TODO average here okay?
    #sigmaBd = BC.reduceSurfaceToCorner(
    #    computeSigmaBd(sigmaBd, sigmaArg, ccArg, 'x', coordinateValueArg1), coordinateValueArg2, coordinateValueArg3)

    for l in indicesMissing:
        oL = BC.getOppositeLatticeDirection(l)

        #tmp1 = -sigmaBd[l] - rhoBd[l] * csArg ** 2 * np.identity(3,dtype=np.double)
        #tmp2 = np.outer(ccArg[oL], ccArg[oL].transpose()) - csArg ** 2 * np.identity(3, dtype=np.double)
        #fRelevant[l] = - fCollRelevant[oL] + 2.0 * wArg[oL] * (rhoBd[l] + 1.0 / (2.0 * csArg ** 4) * (np.tensordot(tmp1,tmp2, axes=2)))
        #fRelevant[l] = (neumannBoundaryRule(sigmaBd2[l], rhoBd[l], csArg, ccArg[oL], wArg[oL], fCollRelevant[oL]) + neumannBoundaryRule(sigmaBd3[l], rhoBd[l], csArg, ccArg[oL], wArg[oL], fCollRelevant[oL]) + dirichletBoundaryRule(jBd[ l], csArg, ccArg[oL], wArg[oL], fCollRelevant[ oL])) / 3.0
        fRelevant[l] = dirichletBoundaryRule(
            jBd[l], csArg, ccArg[oL], wArg[oL], fCollRelevant[oL])
    return fOut


def getCurrentIndicesInXYZOrder( coordinateArg, coordinateValueArg,  index1, index2):
    if(coordinateArg == 'x'):
        return [coordinateValueArg, index1, index2]
    elif (coordinateArg == 'y'):
        return [index1, coordinateValueArg, index2]
    elif (coordinateArg == 'z'):
        return [index1, index2, coordinateValueArg]
    else:
        raise Exception("Invalid input")


def applyDirichletBoundaryConditions(fArg, fCollArg, rhoArg, csArg, ccArg, cArg, wArg, jBdArg,
                                     coordinateArg='x', coordinateValueArg=0):
    # print(rhoArg)
    rhoBd = BC.extrapolateScalarToBd(rhoArg, ccArg, cArg, coordinateArg, coordinateValueArg)
    # dotu = np.zeros(uOldArg.shape, dtype=np.double)
    # for i in range(0, len(uOldArg)):
    #     for j in range(0, len(uOldArg[0])):
    #         for k in range(0, len(uOldArg[0][0])):
    #             dotu[i, j, k] = (uBdArg - uOldArg[i, j, k]) / dtArg
    # dotuRelevant = BC.selectAtCoordinate(dotu, coordinateArg, coordinateValueArg)
    # print(dotuRelevant)
    # print(rhoBd)
    jBd = np.zeros((len(rhoBd), len(rhoBd[0]), 27, 3), dtype=np.double)  # len(rhoArg[0][0])
    # print(jBd)
    for i in range(0, len(jBd)):
        for j in range(0, len(jBd[0])):
            # for k in range(0, len(uBdArg[0][0])):
            # uBd = uBdArg[i,j,k]
            for l in range(0, 27):
                # jBd = dotu * rhoBd
                #jBd[i, j, l] = dotuRelevant[i, j] * rhoBd[i, j, l]
                jBd[i, j, l] = jBdArg

    fCollRelevant = BC.selectAtCoordinate(fCollArg, coordinateArg, coordinateValueArg)
    fOut = copy.deepcopy(fArg)
    fRelevant = BC.selectAtCoordinate(fOut, coordinateArg, coordinateValueArg)
    indicesMissing = BC.getMissingDistributionFunctionIndices(fArg, coordinateArg, coordinateValueArg)
    for i in range(0, len(fRelevant)):
        for j in range(0, len(fRelevant[0])):
            for l in indicesMissing:
                #[i1, i2, i3] = getCurrentIndicesInXYZOrder(coordinateArg, coordinateValueArg, i, j)
                #oL = BC.getOppositeLatticeDirection(
                #    BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]), i1,
                #                                        i2, i3))
                oL = BC.getOppositeLatticeDirection(l)
                fRelevant[i, j, l] = fCollRelevant[i, j, oL] - 2.0 / csArg ** 2 * wArg[oL] * np.tensordot(ccArg[oL],
                                                                                                          jBd[i, j, l],
                                                                                                          axes=1)

    return fOut


def getCurrentIndicesAtEdgeInXYZOrder( coordinateArg1, coordinateValueArg1, coordinateArg2,
                                           coordinateValueArg2, index):
    if(coordinateArg1 == 'x' and coordinateArg2 == 'y'):
        return [coordinateValueArg1, coordinateValueArg2, index]
    elif (coordinateArg1 == 'x' and coordinateArg2 == 'z'):
        return [coordinateValueArg1, index, coordinateValueArg2 ]
    elif (coordinateArg1 == 'y' and coordinateArg2 == 'z'):
        return [ index, coordinateValueArg1, coordinateValueArg2]
    else:
        raise Exception("Invalid input")

def applyDirichletBoundaryConditionsAtEdge(fArg, fCollArg, rhoArg, csArg, ccArg, cArg, wArg, jBdArg,
                                           coordinateArg1='x', coordinateValueArg1=0, coordinateArg2='y',
                                           coordinateValueArg2=0):
    rhoBd = BC.reduceSurfaceToEdge(BC.extrapolateScalarToBd(rhoArg, ccArg, cArg, coordinateArg1, coordinateValueArg1),
                                   coordinateArg1, coordinateArg2, coordinateValueArg2)
    # print(rhoBd.shape)
    # 2D?
    # dotu = np.zeros(uOldArg.shape, dtype=np.double)
    # for i in range(0, len(uOldArg)):
    #     for j in range(0, len(uOldArg[0])):
    #         for k in range(0, len(uOldArg[0][0])):
    #             dotu[i, j, k] = (uBdArg - uOldArg[i, j, k]) / dtArg
    #
    # dotuRelevant = BC.selectAtEdge(dotu, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    # # print(dotuRelevant)
    # # 2D?
    jBd = np.zeros((len(rhoBd), 27, 3), dtype=np.double)  # len(rhoArg[0][0])
    for i in range(0, len(jBd)):
        for l in range(0, 27):
            #jBd[i, l] = dotuRelevant[i] * rhoBd[i, l]
            jBd[i, l] = jBdArg

    fCollRelevant = BC.selectAtEdge(fCollArg, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    fOut = copy.deepcopy(fArg)
    fRelevant = BC.selectAtEdge(fOut, coordinateArg1, coordinateValueArg1, coordinateArg2, coordinateValueArg2)
    indicesMissing = BC.getMissingDistributionFunctionIndicesAtEdge(fArg, coordinateArg1,coordinateValueArg1, coordinateArg2, coordinateValueArg2) ###Alte Version
    #indicesMissing = getMissingDistributionFunctionIndicesAtEdgeDirichlet(fArg, ccArg, cArg, coordinateArg1,
    #                                                                      coordinateValueArg1, coordinateArg2,
    #coordinateValueArg2)
    for i in range(0, len(fRelevant)):
        for l in indicesMissing:
            oL = BC.getOppositeLatticeDirection(l)
            #[i1, i2, i3] = getCurrentIndicesAtEdgeInXYZOrder(coordinateArg1, coordinateValueArg1, coordinateArg2,
            #                                coordinateValueArg2, i)
            #oL = BC.getOppositeLatticeDirection(
            #     BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]),i1,i2,i3))
            fRelevant[i, l] = fCollRelevant[i, oL] - 2.0 / csArg ** 2 * wArg[oL] * np.tensordot(ccArg[oL], jBd[i, l],
                                                                                                axes=1)

    return fOut


def applyDirichletBoundaryConditionsAtCorner(fArg, fCollArg, rhoArg, csArg, ccArg,  wArg, jBdArg,
                                              coordinateValueArg1=0,
                                             coordinateValueArg2=0, coordinateValueArg3=0):
    #rhoBd = BC.reduceSurfaceToCorner(BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, 'x', coordinateValueArg1),
    #                                 coordinateValueArg2, coordinateValueArg3)
    # print(rhoBd)
    # dotu = np.zeros(uOldArg.shape, dtype=np.double)
    # for i in range(0, len(uOldArg)):
    #     for j in range(0, len(uOldArg[0])):
    #         for k in range(0, len(uOldArg[0][0])):
    #             dotu[i, j, k] = (uBdArg - uOldArg[i, j, k]) / dtArg

    #dotuRelevant = dotu[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    # print(dotuRelevant)
    jBd = np.zeros((27, 3), dtype=np.double)  # len(rhoArg[0][0])
    for l in range(0, 27):
        #jBd[l] = dotuRelevant * rhoBd[l]
        jBd[l] = jBdArg

    fCollRelevant = fCollArg[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    fOut = copy.deepcopy(fArg)
    fRelevant = fOut[coordinateValueArg1, coordinateValueArg2, coordinateValueArg3]
    indicesMissing = BC.getMissingDistributionFunctionIndicesAtCorner(fArg, coordinateValueArg1, coordinateValueArg2, coordinateValueArg3) ###Alte Version
    #indicesMissing = getMissingDistributionFunctionIndicesAtCornerDirichlet(fArg, ccArg, cArg, coordinateValueArg1,
    #                                                                        coordinateValueArg2, coordinateValueArg3)
    for l in indicesMissing:
        oL = BC.getOppositeLatticeDirection(l)
        #oL = BC.getOppositeLatticeDirection(BC.getLatticeDirectionForBounceBack(l, ccArg, len(rhoArg), len(rhoArg[0]), len(rhoArg[0][0]), coordinateValueArg1, coordinateValueArg2, coordinateValueArg3))
        fRelevant[l] = fCollRelevant[oL] - 2.0 / csArg ** 2 * wArg[oL] * np.tensordot(ccArg[oL], jBd[l], axes=1)

    return fOut