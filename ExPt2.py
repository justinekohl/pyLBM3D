import numpy as np
import copy
import BoundaryConditions as BC
import Experimental as Ex

#Sigma with SigmaBC
#Corner and Edges not included
def SigmaOnBoundaryConditions(sigma, sigmaBC, SetOfTuples):
    #sigmawithBd = copy.deepcopy(sigma)
    sigmawithBd = np.zeros((len(sigma),len(sigma[0]),len(sigma[0][0]),3,3), dtype=np.double)
    
    sigmaBCX = np.array([[0.0, 0.0, 0.0],
                        [0.0, np.nan, np.nan],
                        [0.0, np.nan, np.nan]])
    
    sigmaBCY = np.array([[np.nan, 0.0, np.nan],
                        [0.0, 0.0, 0.0],
                        [np.nan, 0, np.nan]])

    sigmaBCZ = np.array([[np.nan, np.nan, 0.0],
                        [np.nan, np.nan, 0.0],
                        [0.0, 0.0, 0.0]])

    #edges
    sigmaBCEdgeXY = np.array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, np.nan]])

    sigmaBCEdgeYZ = np.array([[np.nan, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0]])

    sigmaBCEdgeXZ = np.array([[0.0, 0.0, 0.0],
                              [0.0, np.nan, 0.0],
                              [0.0, 0.0, 0.0]])

    # corner
    sigmaBCCorner = np.array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])
    
    # for i in range(0, len(sigmawithBd)):
    #     for j in range(0, len(sigmawithBd[i])):
    #         for k in range(0, len(sigmawithBd[i][j])):
    #             if (i == 0 or i == len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCX)):
    #                     for jj in range(0, len(sigmaBCX[ii])):
    #                         if (not np.isnan(sigmaBCX[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCX[ii,jj]
    #             elif (i != 0 or i != len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCY)):
    #                     for jj in range(0, len(sigmaBCY[ii])):
    #                         if (not np.isnan(sigmaBCY[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCY[ii,jj]
    #             elif (i != 0 or i != len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCZ)):
    #                     for jj in range(0, len(sigmaBCZ[ii])):
    #                         if (not np.isnan(sigmaBCZ[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCZ[ii,jj]
    #             elif (i == 0 or i == len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCEdgeXY)):
    #                     for jj in range(0, len(sigmaBCEdgeXY[ii])):
    #                         if (not np.isnan(sigmaBCEdgeXY[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCEdgeXY[ii,jj]
    #             elif (i == 0 or i == len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCEdgeXZ)):
    #                     for jj in range(0, len(sigmaBCEdgeXZ[ii])):
    #                         if (not np.isnan(sigmaBCEdgeXZ[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCEdgeXZ[ii,jj]
    #             elif (i != 0 or i != len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCEdgeYZ)):
    #                     for jj in range(0, len(sigmaBCEdgeYZ[ii])):
    #                         if (not np.isnan(sigmaBCEdgeYZ[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCEdgeYZ[ii,jj]
    #             elif (i == 0 or i == len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
    #                 for ii in range(0, len(sigmaBCCorner)):
    #                     for jj in range(0, len(sigmaBCCorner[ii])):
    #                         if (not np.isnan(sigmaBCCorner[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBCCorner[ii,jj]
    #             if (i,j,k) in SetOfTuples:
    #                 for ii in range(0, len(sigmaBC)):
    #                     for jj in range(0, len(sigmaBC[ii])):
    #                         if (not np.isnan(sigmaBC[ii,jj])):
    #                             sigmawithBd[i,j,k,ii,jj] = sigmaBC[ii,jj]
    for i in range(0, len(sigmawithBd)):
        for j in range(0, len(sigmawithBd[i])):
            for k in range(0, len(sigmawithBd[i][j])):
                if (i == 0 or i == len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCX
                elif (i != 0 or i != len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCY
                elif (i != 0 or i != len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCZ
                elif (i == 0 or i == len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k != 0 or k != len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCEdgeXY
                elif (i == 0 or i == len(sigmawithBd)) and (j != 0 or j != len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCEdgeXZ
                elif (i != 0 or i != len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCEdgeYZ
                elif (i == 0 or i == len(sigmawithBd)) and (j == 0 or j == len(sigmawithBd[i])) and (k == 0 or k == len(sigmawithBd[i][j])):
                    sigmawithBd[i,j,k] = sigmaBCCorner
                if (i,j,k) in SetOfTuples:
                    sigmawithBd[i,j,k] = sigmaBC
    return sigmawithBd

def applyNeumannBoundaryConditions(fArg, fCollArg, rhoArg, csArg, ccArg, wArg, sigmaBdArg, sigmaArg, coordinateArg='x', coordinateValueArg=0, boundaryRule = Ex.neumannBoundaryRule):
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
        sigmaBd = computeSigmaBd(sigmaBdArg, sigmaArg, ccArg, coordinateArg, coordinateValueArg)
        return [rhoBd, sigmaBd]

    #rhoBd = BC.computeRhoBdWithoutExtrapolation(rhoArg, ccArg, coordinateArg, coordinateValueArg)
    #fCollRelevant = BC.selectAtCoordinate(fCollArg, coordinateArg, coordinateValueArg)
    fOut = copy.deepcopy(fArg)
    #fRelevant = BC.selectAtCoordinate(fOut, coordinateArg, coordinateValueArg)
    #indicesMissing = BC.getMissingDistributionFunctionIndices(fArg, coordinateArg, coordinateValueArg)

    [fCollRelevant, fRelevant, indicesMissing] = Ex.bounceBackFsAtPlane(fOut, fCollArg, coordinateArg, coordinateValueArg)

    #sigmaBd = computeSigmaBd(sigmaBdArg, sigmaArg, ccArg, coordinateArg, coordinateValueArg)

    [rhoBd, sigmaBd] = computeFieldsForBounceBack(rhoArg, ccArg, coordinateArg, coordinateValueArg)

    for i in range(0, len(fRelevant)):
        for j in range(0, len(fRelevant[0])):
            for l in indicesMissing:
                oL = BC.getOppositeLatticeDirection(l)
                fRelevant[i, j, l] = boundaryRule(sigmaBd[i, j, l], rhoBd[i, j, l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, j, oL])

    return fOut

def computeSigmaBd(sigmaBC, sigmaArg, ccArg, coordinateArg='x', coordinateValueArg=0):
    '''

    :param sigmaBC: the prescribed stress at this plane with nan's for undefined values (3,3)
    :param sigmaArg:  the stress field in lattice dimensions (m,n,o,3,3)
    :param ccArg: array of lattice speeds, (27,3)
    :param coordinateArg: 'x', 'y', 'z' the coordinate direction identifying the plane
    :param coordinateValueArg: the index identifying the plane (either 0 or max in the respective direction)
    :return: the stress field for all lattice links at all lattice points of a given plane, accounting for boundary conditions; in plane dimensions (k,l,27,3,3)
    '''
    def applyTractionBoundaryConditions(sigmaBd, sigmaBC, coordinateArg, coordinateValueArg):
        '''

        :param sigmaBd: the current stress field at the plane in plane dimensions (k,l,27,3,3)
        :param sigmaBC: the prescribed stress at this plane with nan's for undefined values (3,3)
        :return: the stress field for all lattice links at all lattice points of a given plane, accounting for boundary conditions; in plane dimensions (k,l,27,3,3)
        '''
        sigmaBCCoordinate = BC.selectAtCoordinate(sigmaBC,coordinateArg,coordinateValueArg)
        
        sigmaBd = copy.deepcopy(sigmaBd)
        for i in range(0, len(sigmaBd)):
            for j in range(0, len(sigmaBd[i])):
                for l in range(0, len(sigmaBd[i][j])):
                    for ii in range(0, len(sigmaBd[i][j][l])):
                        for jj in range(0, len(sigmaBd[i][j][l][ii])):
                            if (not np.isnan(sigmaBCCoordinate[i,j,ii,jj])):
                                sigmaBd[i,j,l,ii,jj] = sigmaBCCoordinate[i,j,ii,jj]
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

    sigmaBd = applyTractionBoundaryConditions(computeSigmaBdWithoutExtrapolation(sigmaArg,ccArg,coordinateArg,coordinateValueArg), sigmaBC, coordinateArg, coordinateValueArg)
    return sigmaBd

def applyNeumannBoundaryConditionsAtEdge(fArg, fCollArg,  rhoArg,  csArg, ccArg,  wArg, sigmaBdArg1, sigmaBdArg2, sigmaArg, coordinateArg1='x', coordinateValueArg1=0, coordinateArg2='y', coordinateValueArg2=0, boundaryRule=Ex.neumannBoundaryRule):
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
    sigmaBd = BC.reduceSurfaceToEdge(computeSigmaBd(sigmaBd, sigmaArg, ccArg, coordinateArg1, coordinateValueArg1), coordinateArg1, coordinateArg2,coordinateValueArg2)

    for i in range(0, len(fRelevant)):
        for l in indicesMissing:
            oL = BC.getOppositeLatticeDirection(l)
            fRelevant[i, l] = boundaryRule(sigmaBd[i,  l], rhoBd[i,  l], csArg, ccArg[oL], wArg[oL], fCollRelevant[i, oL])
    return fOut

def applyNeumannBoundaryConditionsAtCorner(fArg, fCollArg, rhoArg,  csArg, ccArg,  wArg, sigmaBdArg1, sigmaBdArg2, sigmaBdArg3, sigmaArg,  coordinateValueArg1=0, coordinateValueArg2=0, coordinateValueArg3=0, boundaryRule=Ex.neumannBoundaryRule):
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
    sigmaBd = BC.reduceSurfaceToCorner(computeSigmaBd(sigmaBd, sigmaArg, ccArg, 'x', coordinateValueArg1), coordinateValueArg2, coordinateValueArg3)

    for l in indicesMissing:
        oL = BC.getOppositeLatticeDirection(l)
        fRelevant[l] = boundaryRule(sigmaBd[l], rhoBd[l], csArg, ccArg[oL], wArg[oL], fCollRelevant[oL])

    return fOut