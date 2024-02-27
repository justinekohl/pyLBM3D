import numpy as np


def getLatticeVelocitiesWeights(c):
    cc = np.array([[0.0, 0.0, 0.0],  [c, 0.0, 0.0], [-c, 0, 0], [0, c, 0], [0, -c, 0], [0, 0, c], [0, 0, -c],  # 0  - 6
               [c, c, 0], [-c, -c, 0], [c, 0, c], [-c, 0, -c], [0, c, c], [0, -c, -c],  # 7-12
               [c, -c, 0], [-c, c, 0], [c, 0, -c], [-c, 0, c], [0, c, -c], [0, -c, c],  # 13-18
               [c, c, c], [-c, -c, -c], [c, c, -c], [-c, -c, c], [c, -c, c], [-c, c, -c], [-c, c, c], [c, -c, -c]  # 19 - 26
               ], dtype = float)
    w = np.array([8.0/27.0,  # 0
              2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0,  # 1-6
              1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0,  # 7  - 18
              1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0], dtype = float)
    return [cc, w]


# def getLatticeInformation():
#     ax = 1.0
#     maxX = 20
#     maxY = maxX
#     maxZ = maxX
#     dx = ax/maxX  # spacingcl
#     return [dx, maxX, maxY, maxZ]

