import numpy as np
res= 0.25
Lx = 10 #km
Ly = 25
Nx=int(Lx/res)
Ny=int(Ly/res)
cero = np.zeros((Ny,1), dtype=int)
un = np.ones((Ny,1), dtype=int) 
newArray=un
for ix in range(1,Nx-1):
        if ix < (Nx-1)-int(Nx/10) and ix > 0 + int(Nx/10)-1: ap=cero
        else: ap=un
        newArray = np.append(newArray,ap, axis = 1)
print(newArray)

mat = np.matrix(newArray)

with open('mask66.dat','wb') as f:
    for line in mat:
        np.savetxt(f, line, fmt='%.0f')
