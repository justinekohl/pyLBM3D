import numpy as np
import math
import Settings as SettingsModule
import Experimental as Ex
import PostProcessing
import Core
import Util
import os
import ExPt2

nameOfSimulation = "Block3D"
pathToVTK = "./vtk4/"

lam = 1.0
mue = 1.0

rho0 = 1.0
P0 = np.zeros((3, 3))
j0 = np.zeros(3)

ax = 1.0
maxX = 30
maxY = 8
maxZ = 8
dx = ax/maxX  # spacing
xx = np.zeros((maxX, maxY, maxZ, 3), dtype=np.double)
for i in range(0, len(xx)):
    for j in range(0,len(xx[0])):
        for k in range(0, len(xx[0][0])):
            xx[i, j, k] = np.array([np.double(i) * dx, np.double(j) * dx, np.double(k) * dx], dtype=np.double)

cs = math.sqrt(mue/rho0)
dt = 1.0 / math.sqrt(3.0) * dx / cs
c = dx / dt
tau = 0.55 * dt

#print(np.__version__)

#f = np.zeros((m, n, o, 27), dtype = np.double)
#fNew = np.zeros((m, n, o, 27), dtype = np.double)
#fEq = np.zeros((m, n, o, 27), dtype = np.double)


[cc, w] = SettingsModule.getLatticeVelocitiesWeights(c)


[f, j, sigma, u] = Ex.intitialize(rho0, cs, cc, w, maxX, maxY, maxZ, lam, mue)
b = np.zeros((maxX, maxY, maxZ, 3), dtype=np.double)


tMax = 2.0
t = 0.0
m = 0 #gradient for sigma
start = 0.01 #Start Karft drauf geben
end = 2.0 #Ende Karft drauf geben
k = int(0)

ListOfTuples = {(14,7,0), (14,7,1), (14,7,2), (14,7,3), (14,7,4), (14,7,5), (14,7,6), (14,7,7)}

outputFile = None
pointIndices = [[maxX-1, maxY-1, maxZ-1], [maxX-1, 15, 15], [10,10,10]]

while t <= tMax:
    fNew = np.zeros((maxX, maxY, maxZ, 27), dtype=np.double)
    fNew.fill(np.nan)

    rho = Core.computeRho(f)
    jOld = j
    j = Ex.j(Core.firstMoment(f, cc), dt, Ex.firstSource(b, rho0))
    sigma = Ex.sigma(Core.secondMoment(f, cc), dt, Ex.secondSource(Util.dJyDy(j, dx), lam, mue, rho0))
    u = Ex.computeU(u, rho0, j, jOld, dt)

    PostProcessing.writeVTKMaster(k, nameOfSimulation, pathToVTK, t, xx, u, sigma)

    # uOut = []
    # for indices in pointIndices:
    #     uOut.append(u[indices[0], indices[1], indices[2]])

    # outputFile = PostProcessing.writeToFile(outputFile,"./Neumann.dis",uOut,t)

    fEq = Ex.equilibriumDistribution(rho, j, sigma, cc, w, cs, lam, mue, rho0)
    psi = Ex.sourceTermPsi(b, rho0, Util.dJyDy(j, dx), cc, w, cs, mue, lam)

    fColl = Core.collide(f, fEq, psi, dt, tau)
    #print(fNew)
    fNew = Core.stream(fColl, cc, c)
    #print(fNew)

    #Sigma for each lattice point
    if start < t < end:
        sigmaBCxx = np.array([[np.nan, 0, np.nan],
                    [0, -0.0001, 0],
                    [np.nan, 0, np.nan]])
        sigmaBC = ExPt2.SigmaOnBoundaryConditions(sigma, sigmaBCxx, ListOfTuples)
        #m += 1
    else: 
        sigmaBCxx = np.array([[np.nan, 0.0, np.nan],
            [0.0, 0.0, 0.0],
            [np.nan, 0.0, np.nan]])
        sigmaBC =  ExPt2.SigmaOnBoundaryConditions(sigma, sigmaBCxx, ListOfTuples)

    # apply BC at z=0
    fNew = ExPt2.applyNeumannBoundaryConditions(fNew, fColl, rho, cs, cc, w, sigmaBC, sigma, 'z', 0)
    
    # apply BC at z=max
    fNew = ExPt2.applyNeumannBoundaryConditions(fNew, fColl, rho, cs, cc, w, sigmaBC, sigma,
                                             'z', maxZ - 1)
    
    # apply BC at y=0
    fNew = ExPt2.applyNeumannBoundaryConditions(fNew, fColl, rho,  cs, cc, w, sigmaBC, sigma,
                                                             'y', 0)
    
    # apply BC at y=max
    fNew = ExPt2.applyNeumannBoundaryConditions(fNew, fColl, rho, cs, cc, w, sigmaBC, sigma,
                                                             'y', maxY - 1)
    
    # apply BC at x=0
    jBd = np.array([0, 0, 0])
    fNew = Ex.applyDirichletBoundaryConditions(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', 0)
    
    # apply BC at x=max
    fNew = Ex.applyDirichletBoundaryConditions(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', maxX -1)

    #Edges
    
    ### Dirichlet B.C.
    # apply BC at edge x = min, y = min
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', 0, 'y',
                                                     0)

    # apply BC at edge x = min, y = max
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', 0, 'y',
                                                     maxY - 1)

    # apply BC at edge x = min, z = min
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', 0, 'z',
                                                     0)

    # apply BC at edge x = min, z = max
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', 0, 'z',
                                                     maxZ - 1)
    
    # apply BC at edge x = max, y = min
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', maxX -1, 'y', 0)

    # apply BC at edge x = max, y = max
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', maxX -1, 'y', maxY - 1)

    # apply BC at edge x = max, z = min
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', maxX -1, 'z', 0)

    # apply BC at edge x = max, z = max
    fNew = Ex.applyDirichletBoundaryConditionsAtEdge(fNew, fColl, rho, cs, cc, c, w, jBd, 'x', maxX -1, 'z', maxZ - 1)
    
    
    ### Neumann B.C.
    # apply BC at edge y = min, z = min
    fNew = ExPt2.applyNeumannBoundaryConditionsAtEdge(fNew, fColl,  rho,  cs, cc,  w, sigmaBC, sigmaBC, sigma,  'y', 0, 'z', 0)

    # apply BC at edge y = min, z = max
    fNew = ExPt2.applyNeumannBoundaryConditionsAtEdge(fNew, fColl,  rho,  cs, cc,  w, sigmaBC, sigmaBC, sigma,  'y', 0, 'z', maxZ-1)

    # apply BC at edge y = max, z = min
    fNew = ExPt2.applyNeumannBoundaryConditionsAtEdge(fNew, fColl,  rho,  cs, cc,  w, sigmaBC, sigmaBC, sigma,  'y', maxY-1, 'z', 0)

    # apply BC at edge y = max, z = max
    fNew = ExPt2.applyNeumannBoundaryConditionsAtEdge(fNew, fColl,  rho,  cs, cc,  w, sigmaBC, sigmaBC, sigma, 'y', maxY-1, 'z', maxZ-1)
    
    #Corner

    #xmin, ymin, zmin
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  0, 0, 0)
    
    # xmax, ymin, zmin
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  maxX-1, 0, 0)

    # xmax, ymax, zmin
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  maxX-1, maxY-1, 0)

    # xmin, ymax, zmin
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  0, maxY-1, 0)


    #xmin, ymin, zmax
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  0, 0, maxZ - 1)

    # xmax, ymin, zmax
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  maxX-1, 0, maxZ - 1)

    # xmax, ymax, zmax
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  maxX-1, maxY-1, maxZ - 1)

    # xmin, ymax, zmax
    fNew = Ex.applyDirichletBoundaryConditionsAtCorner(fNew, fColl, rho, cs, cc,  w, jBd,  0, maxY - 1, maxZ - 1)
    


    f = fNew

    k = k + 1
    print(k)
    t = t + dt

if outputFile is not None:
    outputFile.close()
print("End")
