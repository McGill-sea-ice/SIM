"""
@author: fiorellab
this function generates ice block 
entry: dimension of domain nx by ny and bc for boundary conditions
output: array +boundary conditions
open water 
"""
import numpy as np
import warnings

def mask_generator(nx, ny, bc):
    
    #adding boundary conditions (suggeré 10 pt de maille à gauche et a droit)
    nx=nx+bc
    ny=ny+bc
    mask=np.zeros((nx,ny))
    mask[1:,bc:-bc] = 1
    
    #message to the user on  mask ratio validity
    if nx/ny!=2.5:
        warnings.warn('choose a ratio nx/ny=2.5')
    return mask

#mask020x050=mask_generator(50,20)
#masknyxnx2=mask_generator(100,40)
#masknyxnx3=mask_generator(250,100)
mask250x100=mask_generator(250,100,10)
