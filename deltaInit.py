
"""
Created on Mon Sep 10 16:26:55 2018

@author: bmusammartanov
"""
print("\n")
print("#########################################################")
print("######                                            #######")
print("######      Base Mesh intial cell dimension       #######")
print("######                                            #######")
print("#########################################################")
print("\n")
xmin = float(input("- Insert Xmin [m]: "))
xmax = float(input("- Insert Xmax [m]: "))
ymin = float(input("- Insert Ymin [m]: "))
ymax = float(input("- Insert Ymax [m]: "))
zmin = float(input("- Insert Zmin [m]: "))
zmax = float(input("- Insert Zmax [m]: "))
Nx= int(input("- Number of cells along x direction : "))
Ny= int(input("- Number of cells along y direction : "))
Nz= int(input("- Number of cells along z direction : "))

def deltaDist():
    dx = (xmax - xmin)/Nx 
    dy = (ymax - ymin)/Ny
    dz = (zmax - zmin)/Nz
    dAvg = (dx+dy+dz)/3  
    print("\n--> The average initial cell lenght is {:5.5f} m".format(dAvg))
    
deltaDist()

ans = input("Press any key to close!")
