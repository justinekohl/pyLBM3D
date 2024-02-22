# Beispiel von
# LBMS/examples/solid lattice boltzmann dynamics/soliddynamics.cpp

import numpy as np
import math
import Settings as SettingsModule
import Experimental as Ex
import PostProcessing
import Core
import Util

import QS.QuasiStatic as QS
import QS.QuasiStaticBC as QSBC
import QS.ComputeGrad as Grad
import QS.ComputeDiv as Div

nameOfSimulation = "Block3D"
pathToVTK = "E:/git/vtk/"

dt = 0.000001 #Weicht ab vom Original, um Zeit zu sparen
tMax = 0.0002
nodesPerMeter = 50 # erster Versuch mit nur 50
lx = 1.0
ly = 0.2
lz = 0.2
omega = 0.001
rho0 = 7860 # kg/m^3
moduleYoung = 200.0*pow(10.0, 9.0) # E-Modul N/m^2
nu = 0.3 #poisson's ratio

lam = (moduleYoung*nu)/((1.0+nu)*(1.0-2.0*nu))
mue = moduleYoung/(2.0*(1.0+nu))
#constitutive law for isotropic solids: hooke's law with lame coefficients
# moduleYoungLatticeUnits = moduleYoung*dx*dt*dt
# lam = (moduleYoungLatticeUnits*nu)/((1.0+nu)*(1.0-2.0*nu))
# mue = moduleYoungLatticeUnits/(2.0*(1.0+nu))

P0 = np.zeros((3, 3))
j0 = np.zeros(3)

# setting up lattice point positions
dx = 1.0/nodesPerMeter # spacing
maxX = int(lx/dx)
maxY = int(ly/dx)
maxZ = int(lz/dx)
xx = np.zeros((maxX, maxY, maxZ, 3), dtype=np.double)
for i in range(0, len(xx)):
    for j in range(0,len(xx[0])):
        for k in range(0, len(xx[0][0])):
            xx[i, j, k] = np.array([np.double(i) * dx, np.double(j) * dx, np.double(k) * dx], dtype=np.double)

cs = math.sqrt(mue/rho0)
c = dx / dt

tau = 1.0 / omega
#tau = 0.55 * dt

[cc, w] = SettingsModule.getLatticeVelocitiesWeights(c)

[f, j, sigma, u] = QS.intitialize(rho0, w, maxX, maxY, maxZ)
divSigma = Div.ComupteDivergence(sigma,dx)
#divSigma = QS.divOfSigma(sigma,dx)
rho = Core.computeRho(f)
#b = np.zeros((maxX, maxY, maxZ, 3), dtype=np.double) # TODO not needed

t = 0.0
k = int(0)


outputFile = None
pointIndices = [[maxX-1, maxY-1, maxZ-1], [maxX-1, maxY-1, maxZ-1], [0,0,0]]

while t <= tMax:
    fNew = np.zeros((maxX, maxY, maxZ, 27), dtype=np.double)
    fNew.fill(np.nan)

    # Lattice forces
    g = QS.g(rho, divSigma)
    v = QS.v(rho,j)
    gi = QS.gi(g,cc,w,rho,cs,v) # TODO this is cs = 1/sqrt(3)

    # BC ####################################
    def uBdFromCoordinates(coordinates):
        clocal = 0.001  # TODO maybe ramping
        #nu = lam / (2.0 * (lam + mue))
        #nu = 0.15 #poisson's ratio
        return np.array([ clocal * coordinates[0], - nu * clocal * coordinates[1], - nu * clocal * coordinates[2] ])

    visited = np.zeros((maxX, maxY, maxZ), dtype=bool)
    # TODO for fixed BC use bounce back, halfway?
    #xmin
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0,rho, w,u,uBdFromCoordinates,visited,'x',0)

    #xmax
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0, rho, w,u,uBdFromCoordinates,visited,'x',maxX-1)

    #ymin
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0, rho, w,u,uBdFromCoordinates,visited,'y',0)

    #ymax
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0, rho, w,u,uBdFromCoordinates,visited,'y',maxY-1)

    #zmin
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0, rho, w,u,uBdFromCoordinates,visited,'z',0)

    #zmax
    [f,u, visited] = QSBC.applyDirichletBoundaryConditions(f,dx,dt,rho0, rho, w,u,uBdFromCoordinates,visited,'z',maxZ-1)
    
    # End BC #################################



    # collision and streaming
    fEq = QS.equilibriumDistribution(rho,w)

    fColl = QS.collide(f,fEq,gi,dt,omega)
    fNew = Core.stream(fColl, cc, c)
    f = fNew

    # compute displacement
    jOld = j
    j = Ex.j(Core.firstMoment(f, cc), dt, g)
     # compute rho
    rho = Core.computeRho(f)
    u = QS.computeU(u, rho0, j, jOld, dt, v, rho) # TODO rho0 more stable or use rho?
    # TODO prevent overwriting at Boundary, why ux not constant for x = const.?
   
    # compute strain, stress, stress divergence
    #gradU = Util.computeGradientU(u,dx)
    gradU = Grad.ComupteGradientU(u,dx)
    sigma = QS.sigmaFromDisplacement(gradU,lam,mue)
    #divSigma = QS.divOfSigma(sigma,dx)
    divSigma = Div.ComupteDivergence(sigma,dx)

    # Postprocessing ################
    PostProcessing.writeVTKMaster(k, nameOfSimulation, pathToVTK, t, xx, u, sigma)
    uOut = []
    for indices in pointIndices:
        uOut.append(u[indices[0], indices[1], indices[2]])

    outputFile = PostProcessing.writeToFile(outputFile,"./DirichletQS.dis",uOut,t)
    ######################

    k = k + 1
    print(k)
    t = t + dt

if outputFile is not None:
    outputFile.close()
print("End")

