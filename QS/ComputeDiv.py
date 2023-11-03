import numpy as np

def ComupteDivergence(stressStartArg, deltaxArg):
    
    divergence = np.zeros((len(stressStartArg), len(stressStartArg[0]), len(stressStartArg[0][0]),3), dtype=np.double)
    
    # vicinity check

    for XArg in range(0, len(stressStartArg)):
        for YArg in range(0, len(stressStartArg[0])):
            for ZArg in range(0, len(stressStartArg[0][0])):
                
                if (XArg != 0 and YArg != 0 and ZArg != 0 and XArg != len(stressStartArg)-1 and YArg != len(stressStartArg[0])-1 and ZArg != len(stressStartArg[0][0])-1):
                    
                    # highly accurate central difference scheme + O(DeltaX^4)
                    if (XArg != 1 and YArg != 1 and ZArg != 1 and XArg != len(stressStartArg)-2 and YArg != len(stressStartArg[0])-2 and ZArg != len(stressStartArg[0][0])-2):
                        divStressStartCell0First = (8.0*stressStartArg[XArg+1][YArg][ZArg] - 8.0*stressStartArg[XArg-1][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(12*deltaxArg)
                        divStressStartCell0Second = (8.0*stressStartArg[XArg][YArg+1][ZArg] - 8.0*stressStartArg[XArg][YArg-1][ZArg] - stressStartArg[XArg][YArg+2][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(12*deltaxArg)
                        divStressStartCell0Third = (8.0*stressStartArg[XArg][YArg][ZArg+1] - 8.0*stressStartArg[XArg][YArg][ZArg-1] - stressStartArg[XArg][YArg][ZArg+2] + stressStartArg[XArg][YArg][ZArg-2])/(12*deltaxArg)
                        
                        #divergence[XArg+1][YArg][ZArg][0] = divStressStartCell0First[0][0] + divStressStartCell0Second[0][1] + divStressStartCell0Third[0][2]
                        #divergence[XArg+1][YArg][ZArg][1] = divStressStartCell0First[1][0] + divStressStartCell0Second[1][1] + divStressStartCell0Third[1][2]
                        #divergence[XArg+1][YArg][ZArg][2] = divStressStartCell0First[2][0] + divStressStartCell0Second[2][1] + divStressStartCell0Third[2][2]

                    else: 
                    #central difference scheme + O(DeltaX^2)
                        divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                        divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                        divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)


                else:
                    # depending on boundary conditions

                    #x=0
                    if XArg == 0:
                        if (YArg != 0 and YArg != len(stressStartArg[0])-1 and ZArg != 0 and ZArg != len(stressStartArg[0][0])-1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                            #divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                            #divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)

                        #y=0
                        elif (YArg == 0 and YArg != len(stressStartArg[0])-1 and ZArg != 0 and ZArg != len(stressStartArg[0][0])-1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                            
                            #divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point forward Nx Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Nx Ny scheme O(DeltaX^2)

                        #z=0
                        elif (YArg != 0 and YArg != len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                        
                        #y=ymax
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                            
                            #divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX][globalY-1]+stressStart[2][globalX][globalY-2])/(2*DeltaX); // three-point forward Nx backward Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (3.0*stressStart[1][globalX][globalY]-4.0*stressStart[1][globalX][globalY-1]+stressStart[1][globalX][globalY-2])/(2*DeltaX); // three-point forward Nx backward Ny scheme O(DeltaX^2)

                        #z=zmax
                        elif (YArg != 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                        
                        ## Ecken
                        #x=y=z=0
                        elif (YArg == 0 and YArg != len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                        
                        #x=z=0 and y=max 
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #x=y=0 and z=max 
                        elif (YArg == 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                        #x=0 and y=z=max 
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (4.0*stressStartArg[XArg+1][YArg][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg+2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                    
                    #x=xmax    
                    elif XArg == len(stressStartArg):
                        
                        if(YArg != 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                            #divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                            #divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                        
                        #y=0
                        elif (YArg == 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                            #divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point backward Nx forward Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point backward Nx forward Ny scheme O(DeltaX^2)
                            
                        #z=0
                        elif (YArg != 0 and YArg != len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #y=ymax
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                            
                            #divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX][globalY-1]+stressStart[2][globalX][globalY-2])/(2*DeltaX); // three-point backward Nx Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (3.0*stressStart[1][globalX][globalY]-4.0*stressStart[1][globalX][globalY-1]+stressStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Nx Ny scheme O(DeltaX^2)
                        
                        #z=zmax
                        elif (YArg != 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                        ## Ecken
                        #x=max and y=z=0
                        elif (YArg == 0 and YArg != len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                        
                        #x=max=y and z=0
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #x=max=z and y=0
                        elif (YArg == 0 and YArg != len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                        #x=y=z=max 
                        elif (YArg != 0 and YArg == len(stressStartArg[0]) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg-1][YArg][ZArg] + stressStartArg[XArg-2][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                    
                    #y=0
                    elif YArg == 0:  
                         
                        if (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                             
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)

                            #divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                            #divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                        
                        #z=0
                        elif (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #z=zmax    
                        elif (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (4.0*stressStartArg[XArg][YArg+1][ZArg] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg+2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)

                    #y=ymax
                    elif YArg == len(stressStartArg[0]) - 1:
                        
                        if (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg != 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (stressStartArg[XArg][YArg][ZArg+1] - stressStartArg[XArg][YArg][ZArg-1])/(2*deltaxArg)
                         
                        #z=0
                        elif (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg == 0 and ZArg != len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                            
                        #z=zmax    
                        elif (XArg != 0 and XArg != len(stressStartArg) - 1 and ZArg != 0 and ZArg == len(stressStartArg[0][0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg-1][ZArg] + stressStartArg[XArg][YArg-2][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                    
                    #z=0      
                    elif ZArg == 0:    
                        if (XArg != 0 and XArg != len(stressStartArg) - 1 and YArg != 0 and YArg != len(stressStartArg[0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg+1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (4.0*stressStartArg[XArg][YArg][ZArg+1] - 3.0*stressStartArg[XArg][YArg][ZArg] - stressStartArg[XArg][YArg][ZArg+2])/(2*deltaxArg)
                    
                    #z=zmax    
                    elif ZArg == len(stressStartArg[0][0]):
                        if (XArg != 0 and XArg != len(stressStartArg) - 1 and YArg != 0 and YArg != len(stressStartArg[0]) - 1):
                            
                            divStressStartCell0First = (stressStartArg[XArg+1][YArg][ZArg] - stressStartArg[XArg-1][YArg][ZArg])/(2*deltaxArg)
                            divStressStartCell0Second = (stressStartArg[XArg][YArg-1][ZArg] - stressStartArg[XArg][YArg-1][ZArg])/(2*deltaxArg)
                            divStressStartCell0Third = (3.0*stressStartArg[XArg][YArg][ZArg] - 4.0*stressStartArg[XArg][YArg][ZArg-1] + stressStartArg[XArg][YArg][ZArg-2])/(2*deltaxArg)
                            
                divergence[XArg][YArg][ZArg][0] = divStressStartCell0First[0][0] + divStressStartCell0Second[0][1] + divStressStartCell0Third[0][2]
                divergence[XArg][YArg][ZArg][1] = divStressStartCell0First[1][0] + divStressStartCell0Second[1][1] + divStressStartCell0Third[1][2]
                divergence[XArg][YArg][ZArg][2] = divStressStartCell0First[2][0] + divStressStartCell0Second[2][1] + divStressStartCell0Third[2][2]

    
    return divergence
