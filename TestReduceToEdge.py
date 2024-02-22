import unittest
import numpy as np
import BoundaryConditions as BC

# Create a 3D array with integer numbers
three_dimensional_array = np.array([[[[1], [2], [3]],
                                     [[4], [5], [6]],
                                     [[7], [8], [9]]],

                                    [[[10], [11], [12]],
                                     [[13], [14], [15]],
                                     [[16], [17], [18]]],

                                    [[[19], [20], [21]],
                                     [[22], [23], [24]],
                                     [[25], [26], [27]]]])

def compare_edges(dataAtEdge, expectedDataAtEdge):
  return dataAtEdge[0] == expectedDataAtEdge[0] and dataAtEdge[1] == expectedDataAtEdge[1] and dataAtEdge[2] == expectedDataAtEdge[2]

# Front left
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","z",0)

fl = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[1,4,7])


# Front right
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","z",2)

fr = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[3,6,9])


# Front top
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","y",0)

ft = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[1,2,3])


# Front bottom
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","y",2)

fb = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[7,8,9])

# Back left
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","z",0)

bl = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[19,22,25])

# Back right
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","z",2)

br = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[21,24,27])


# Back top
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","y",0)

bt = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[19,20,21])

# Back bottom
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"x",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"x","y",2)

bb = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[25,26,27])



# Left Bottom
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"y",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"y","z",0)

lb = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[7,16,25])

# Left Top
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"y",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"y","z",0)

lt = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[1,10,19])



# Right Bottom
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"y",2)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"y","z",2)

rb = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[9,18,27])

# Right Top
dataAtSurface = BC.selectAtCoordinate(three_dimensional_array,"y",0)
dataAtEdge = BC.reduceSurfaceToEdge(dataAtSurface,"y","z",2)

rt = compare_edges(dataAtEdge=dataAtEdge, expectedDataAtEdge=[3,12,21])

assert(fl and fr and ft and fb and bl and br and bt and bb and lb and lt and rb and rt)
print("<<< Test passed >>>")





