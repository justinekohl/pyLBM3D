import numpy as np
import math
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt


#Only for x == MaxX, y, z
def twist_angle(uArg, xxArg, m):

    punktez = []
    punktey = []
    
    middle_y = (len(xxArg[0]) - 1)/2
    
    maxX = len(xxArg) - 1
    
    if middle_y % 2:
        rangey = [math.floor(middle_y), math.ceil(middle_y)]
        for j in rangey:
            for k in range(0, len(xxArg[0][0])):
                real_position = xxArg[maxX][j][k]
                u_InThisPoint = uArg[maxX][j][k]
                punktey.append(real_position[1] + u_InThisPoint[1])
                punktez.append(real_position[2] + u_InThisPoint[2])
                
    else:
        for k in range(0, len(xxArg[0][0])):
            real_position = xxArg[maxX][middle_y][k]
            u_InThisPoint = uArg[maxX][middle_y][k]
            punktey.append(real_position[1] + u_InThisPoint[1])
            punktez.append(real_position[2] + u_InThisPoint[2])
        
    
    z = np.array(punktez).reshape(-1,1)
    y = np.array(punktey)

    model = LinearRegression().fit(z, y)

    alpha  = np.arctan(model.coef_[0]) * (180/np.pi)
    
    plt.scatter(z,y)
    plt.plot(z,(model.coef_)*z+(model.intercept_),color='r')
    plt.title('Verdrehwinkel:'+str(round(alpha,3)))
    #plt.show()
    plt.savefig('E:\git2.0\savepngress\Regressionsgerade'+str(m)+'.png')
    plt.close()
    m+=1

    return m, alpha