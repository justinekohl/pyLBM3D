import numpy as np

def ComupteGradientU(UArg, deltaxArg):
    gradient = np.zeros((len(UArg), len(UArg[0]), len(UArg[0][0]),3,3), dtype=np.double)

    for XArg in range(0, len(UArg)):
        for YArg in range(0, len(UArg[0])):
            for ZArg in range(0, len(UArg[0][0])):

                if (XArg != 0 and YArg != 0 and ZArg != 0 and XArg != len(UArg)-1 and YArg != len(UArg[0])-1 and ZArg != len(UArg[0][0])-1):
                    
                    #highly accurate central difference scheme + O(DeltaX^4)
                    if (XArg != 1 and YArg != 1 and ZArg != 1 and XArg != len(UArg)-2 and YArg != len(UArg[0])-2 and ZArg != len(UArg[0][0])-2):
                        
                        gradU_x = (8.0*UArg[XArg+1][YArg][ZArg] - 8.0*UArg[XArg-1][YArg][ZArg] - UArg[XArg+2][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(12*deltaxArg)
                        gradU_y = (8.0*UArg[XArg][YArg+1][ZArg] - 8.0*UArg[XArg][YArg-1][ZArg] - UArg[XArg][YArg+2][ZArg] + UArg[XArg][YArg-2][ZArg])/(12*deltaxArg)
                        gradU_z = (8.0*UArg[XArg][YArg][ZArg+1] - 8.0*UArg[XArg][YArg][ZArg-1] - UArg[XArg][YArg][ZArg+2] + UArg[XArg][YArg][ZArg-2])/(12*deltaxArg)

                    #central difference scheme + O(DeltaX^2)
                    else:
                        gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                        gradU_y = (UArg[XArg][YArg+1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                        gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                else:
                    # depending on boundary conditions

                    #x=0
                    if XArg == 0:
                        if (YArg != 0 and YArg != len(UArg[0])-1 and ZArg != 0 and ZArg != len(UArg[0][0])-1):
                            
                            gradU_x = (4.0*UArg[XArg+1][YArg][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (UArg[XArg][YArg+1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                        #y=0
                        elif (YArg == 0 and YArg != len(UArg[0])-1 and ZArg != 0 and ZArg != len(UArg[0][0])-1):

                            gradU_x = (4.0*UArg[XArg+1][YArg][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                        #z=0
                        elif (YArg != 0 and YArg != len(UArg[0]) - 1 and ZArg == 0 and ZArg != len(UArg[0][0]) - 1):

                            gradU_x = (4.0*UArg[XArg+1][YArg][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (UArg[XArg][YArg+1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)

                        #y=ymax
                        elif (YArg != 0 and YArg == len(UArg[0]) - 1 and ZArg != 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (4.0*UArg[XArg+1][YArg][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                        
                        #z=zmax
                        elif (YArg != 0 and YArg != len(UArg[0]) - 1 and ZArg != 0 and ZArg == len(UArg[0][0]) - 1):
                            
                            gradU_x = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg-1][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (UArg[XArg][YArg+1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                        ## Ecken
                        #x=max and y=z=0
                        elif (YArg == 0 and YArg != len(UArg[0]) - 1 and ZArg == 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg-1][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                        
                        #x=max=y and z=0
                        elif (YArg != 0 and YArg == len(UArg[0]) - 1 and ZArg == 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg-1][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #x=max=z and y=0
                        elif (YArg == 0 and YArg != len(UArg[0]) - 1 and ZArg != 0 and ZArg == len(UArg[0][0]) - 1):
                            
                            gradU_x = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg-1][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                        #x=y=z=max 
                        elif (YArg != 0 and YArg == len(UArg[0]) - 1 and ZArg != 0 and ZArg == len(UArg[0][0]) - 1):
                            
                            gradU_x = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg-1][YArg][ZArg] + UArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                    
                    #y=0
                    elif YArg == 0:  
                         
                        if (XArg != 0 and XArg != len(UArg) - 1 and ZArg != 0 and ZArg != len(UArg[0][0]) - 1):
                             
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                            #divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                        
                        #z=0
                        elif (XArg != 0 and XArg != len(UArg) - 1 and ZArg == 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #z=zmax    
                        elif (XArg != 0 and XArg != len(UArg) - 1 and ZArg != 0 and ZArg == len(UArg[0][0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (4.0*UArg[XArg][YArg+1][ZArg] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)

                    #y=ymax
                    elif YArg == len(UArg[0]) - 1:
                        
                        if (XArg != 0 and XArg != len(UArg) - 1 and ZArg != 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (UArg[XArg][YArg][ZArg+1] - UArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                         
                        #z=0
                        elif (XArg != 0 and XArg != len(UArg) - 1 and ZArg == 0 and ZArg != len(UArg[0][0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #z=zmax    
                        elif (XArg != 0 and XArg != len(UArg) - 1 and ZArg != 0 and ZArg == len(UArg[0][0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg-1][ZArg] + UArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                    
                    #z=0      
                    elif ZArg == 0:    
                        if (XArg != 0 and XArg != len(UArg) - 1 and YArg != 0 and YArg != len(UArg[0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (UArg[XArg][YArg+1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            gradU_z = (4.0*UArg[XArg][YArg][ZArg+1] - 3.0*UArg[XArg][YArg][ZArg] - UArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                    
                    #z=zmax    
                    elif ZArg == len(UArg[0][0]):
                        if (XArg != 0 and XArg != len(UArg) - 1 and YArg != 0 and YArg != len(UArg[0]) - 1):
                            
                            gradU_x = (UArg[XArg+1][YArg][ZArg] - UArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            gradU_y = (UArg[XArg][YArg-1][ZArg] - UArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            gradU_z = (3.0*UArg[XArg][YArg][ZArg] - 4.0*UArg[XArg][YArg][ZArg-1] + UArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                        
                        
                gradient[XArg][YArg][ZArg][0][0] = gradU_x[0]
                gradient[XArg][YArg][ZArg][0][1] = gradU_y[0]
                gradient[XArg][YArg][ZArg][0][2] = gradU_z[0]
                gradient[XArg][YArg][ZArg][1][0] = gradU_x[1]
                gradient[XArg][YArg][ZArg][1][1] = gradU_y[1]
                gradient[XArg][YArg][ZArg][1][2] = gradU_z[1]
                gradient[XArg][YArg][ZArg][2][0] = gradU_x[2]
                gradient[XArg][YArg][ZArg][2][1] = gradU_y[2]
                gradient[XArg][YArg][ZArg][2][2] = gradU_z[2]
                
    return gradient