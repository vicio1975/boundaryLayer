#-*- coding: utf-8 -*-
"""
Created on 13/05/2015
Updated on 16/08/2018

Changes: 1. Removed the moving results file (WIN users)
         2. Added the selection of geometry duct, three geometries kind
         3. Added the skin friction factor pipe flows
         4. Limit the Espansion number and the yplus_max number as "yplus_min + 50"
         Small changes(3): mesh details in text, format text, unregular geoms, Region Name, versions
                        turbulent length scale.
         Small changes(4):   air viscosity estimation "Sutherland Equation",
         Vogel water viscosity estimation.
         5. Termal boundary layer implementation (work in progress): 
            Prandtl Number = kinematic_viscosity/thermal_diffusivity = [μ/ρ] / [Kc/(ρ*cp)] = 
            Prandtl Number = (μ*cp)/Kc 
            ν : momentum diffusivity (kinematic viscosity), ν = μ/ρ (SI units: m2/s)
            α : thermal diffusivity, α = Kc/(ρ*cp) (SI units: m2/s)
            μ : dynamic viscosity, (SI units: Pa s = N s/m2)
            Kc: thermal conductivity, (SI units: W/m-K)
            cp: specific heat, (SI units: J/kg-K)
            ρ : density, (SI units: kg/m3).
            
            1. boundary layer improvemnt estimation;
            2. saving the selected y+ value; 
            
@author: Vincenzo Sammartano
email: v.sammartano@gmail.com
"""

import numpy as num
#from shutil import move
#import os

#### Version of the tool
V = 5
Sv = 2

###################################################Classes declarations
class physics:
    fluids = ['air','water']
    nameflu = "Not assigned"
    #Based on a mean molar mass for dry air of 28.9645 g/mol.
    Rf = 287.058 # Universal Constant of Gases [J/(Kg K)]
    pAtm = 101325 # [Pa] atmospheric pressure
    g = 9.806   # [m/s2] gravitational accelaration

    #initialization of nameflu
    def __init__(self, nameflu):  # nameflu = flu
        self.nameflu = nameflu    # instance variable unique to each instance

    def prop(self):
        t = self.T + 273.17 # Kelvin
        rep = True
        while rep == True:
            #Air
            if self.nameflu == self.fluids[0]:
                rot = self.pAtm/(self.Rf*t) #Density as function of temperature [Kg/mc]
                gamma_t = rot * self.g  #specific weight at t°C
                  #Sutherland Equation
                ba = 1.458*10**(-6) 
                sa = 110.4 #costant in Kelvin
                mi = ba * (t**1.5)/(t+sa) #Dinamic Viscosity  Pa s = kg m^-1 s^-1
                ni = mi/rot         #Cinematic Viscosity  m2·s-1
                
                  #Thermal diffusivity, α = Kc/(ρ*cp)
                # -183 < T < 218 C
                #Thermal conductivity
                Kc = 5.75e-5 * ( 1 + 0.00317 * (self.T) - 0.0000021 * (self.T**2))
                #Specific heat
                cp = (-10**-10*(self.T)**3) + ( 3*10**-7*(self.T)**2 ) - (5*10**-5*(self.T)) + 0.9917
                #cv = (-10**-10*(self.T)**3) + ( 3*10**-7*(self.T)**2 ) - (5*10**-5*(self.T)) + 0.7047
                #thermal diffusivity
                alpha = Kc/(rot*cp) 
                PrN = ni/alpha
                rep = False
                return [rot,gamma_t,mi,ni,alpha,PrN]
            #Water
            elif self.nameflu == self.fluids[1]:
                #Kell formulation
                ro = 999.83952
                At = 1 + (16.879850*10**-3)*self.T #constant to be used in the Kell formula

                #rot is the water density as function of temperature [Kg/mc]
                rot =  (ro + (16.945176*self.T) - (7.987040*10**-3)*(self.T**2.0) -
                     +(46.17046*10**-6)*(self.T**3.0) + (105.56302*10**-9)*(self.T**4.0)-
                            +(280.54253*10**-12)*(self.T**5.0))/At
                gamma_t = rot*self.g #specific weight at t°C
                #Vogel Formula
                # 0<T<100 Celsius
                c1 = -3.7188
                c2 = 578.919
                c3 = -137.546                 
                mi = (num.e**((c1+(c2/(t+c3)))))/1000 #Dinamic Viscosity  Pa s = kg m^-1 s^-1
                ni = mi/rot         #Cinematic Viscosity  m2·s-1

                #Thermal conductivity
                Kc = 5.7109 * (10**-1) + 1.7625 * (10**-3) * self.T - 6.7036 * 10**(-6) * (self.T **2)
                #Thermal Diffusivity
                alpha = 1.3168 * (10**-7) + 6.2477 * (10**-10) * self.T - 2.4022 * (10**-12) * (self.T**2)
                #Specific heat
                cp = Kc/(rot*alpha)
                PrN = ni/alpha
                rep = False
                return [rot,gamma_t,mi,ni,alpha,PrN]
            else:
                print(" ... please select the correct fluid from the list!")
                rep = True

################################################### Functions here
def spec(V,Sv):
    print("#########################################################")
    print("######                                            #######")
    print("######          BOUNDARY LAYER ESTIMATION         #######")
    print("######            - Fluid Specifications          #######")
    print("######            - Flow Velocity Field           #######")
    print("######            - Turbulence characteristics    #######")
    print("######            - Mesh design                   #######")
    print("######                                    ver.{}.{} #######".format(V,Sv))
    print("#########################################################")

    regionName = input("\n>>> Set the domain region = ") 
    print("\n---------------- Fluid Specifications -----------------")

    fluSel = {'1': 'air', '2': 'water'} #dictionary
    print("\n* Fluid to be modelled:\n",
          "  1. Air\n",
          "  2. Water")
    flu = 'None'
    while flu == "None":
        flu = fluSel.get(input("* Set the fluid [1-2]: "),"None")
        if flu == "None": print("... fluid not yet implemented!")

    fluid = physics(flu) # create the object fluid
    #Fluid vars out --> [0-rot, 1-gamma_t, 2-mi, 3-ni, 4-alpha, 5-PrN]

    fluid.T = float(input("* Set the operating temperature (-30 C < T < 150 C): "))
    print("\n--> The {one} has a temperature T = {two} [C]".format(one=fluid.nameflu,two=fluid.T))
    print("--> The density of {one} is {two:1.4f} Kg/m^3".format(one=fluid.nameflu,two = fluid.prop()[0]))
    print("--> The specific  volume of {one} is {two:1.2f} N/m^3".format(one=(fluid.nameflu),two=fluid.prop()[1]))
    print("--> The dinamic   viscosity of {one} is {two:1.4e} Pa*s".format(one=fluid.nameflu,two=fluid.prop()[2]))
    print("--> The kinematic viscosity of {one} is {two:1.4e} m^2/s".format(one=fluid.nameflu,two=fluid.prop()[3]))
    #Thermal 
    print("--> The thermal diffusivity of {one} is {two:1.4e} m^2/s".format(one=fluid.nameflu,two=fluid.prop()[4]))
    print("--> The Prandtl Number of {one} is {two:1.4e} m^2/s".format(one=fluid.nameflu,two=fluid.prop()[5]))
    prn = fluid.prop()[5]
    if prn > 1:
        print("--> Prandtl Number > 1 ==> Velocity BL > Thermal BL")
    else:
        print("--> Prandtl Number < 1 ==> Velocity BL < Thermal BL")
    return fluid,regionName

def geom(V0,fluid):
    print("* Flow geometry conditions:")
    ans = "None"
    while ans == "None":
        ans = int(input("    1. External Flow\n    2. Confined Flow\n* Set flow geometry [1-2]: "))

        #### Free stream
        if ans == 1:
            FLK = "ext"
            dc = float(input("* Set the lenght scale: "))
            RE = V0 * (dc/fluid.prop()[3])
            #Hansen(1928) approach
            if RE < 5*10**5:
                print("--> Laminar Regime")
                d1 = 4.91
                d2 = 1/2
                delta = d1 * dc/(RE**d2) ##Boundary Layer Thickness
                print("--> Boundary Layer Thickness: {:8.6f} m".format(delta))
            elif RE >= 5*10**5:
                print("--> Turbulent Regime")
                d1 = 0.37
                d2 = 1/5
                delta = d1 * dc/(RE**d2) ##Boundary Layer Thickness
                print("--> Boundary Layer Thickness: {:8.6f} m".format(delta))
            ### delta is estimated
            
        #### Confined Flows
        elif ans == 2:
            FLK = "int"
            ansS = int(input("    1. Circular Duct\n    2. Rectangular Duct\n    3. Unregular shape\n* Set the duct geometry [1-2-3]: "))
            ## case 1. Circular section
            if ansS == 1:
                dc = float(input("* Set duct diameter (m): "))
                L  = float(input("* Set the duct lenght (m): "))
                RE_0 = V0 * (dc/fluid.prop()[3])
                if RE_0 < 3000:
                    Lh  = 0.05 * RE_0 * dc
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")
                        
                elif RE_0 > 3000 and RE_0 < 10000:
                    Lh  = 4 * dc * RE_0**(1/6)
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")

                elif RE_0 > 10000:
                    Lh  = 40 * dc #Nikuradse
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")
                                        
            ## case 2. Rectangular section
            if ansS == 2:
                sA = float(input("* Set the width (m): "))
                sB = float(input("* Set the height (m): "))
                Ar = sA * sB   # area of the section
                Pr = 2 * (sA + sB) # perimeter of the section
                dc = 4 * (Ar/Pr) #hydraulic diameter
                L  = float(input("* Set the duct lenght (m): "))
                RE_0 = V0 * (dc/fluid.prop()[3])
                #Check of the boundary layer development
                if RE_0 < 3000:
                    Lh  = 0.05 * RE_0 * dc
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")
                        
                elif RE_0 > 3000 and RE_0 < 10000:
                    Lh  = 4 * dc * RE_0**(1/6)
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")

                elif RE_0 > 10000:
                    Lh  = 40 * dc #Nikuradse
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")

           ## case 3. Unregular shape
            if ansS == 3:
                Ar = float(input("* Set the area (sm): "))
                Pr = float(input("* Set the perimeter (m): "))
                dc = 4 * (Ar/Pr) #hydraulic diameter
                L  = float(input("* Set the duct lenght (m): "))
                RE_0 = V0 * (dc/fluid.prop()[3])
   
                #Check of the boundary layer development
                if RE_0 <= 2300:
                    Lh  = 0.05 * RE_0 * dc #Kays and Crawford (1993) and Shah and Bhatti (1987)
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")
                        
                elif RE_0 > 2300 and RE_0 < 10000:
                    Lh  = 1.359 * dc * RE_0**(1/4) #Bhatti and Shah (1987) and Zhi-qing (1982)
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")

                elif RE_0 > 10000:
                    Lh  = 10 * dc 
                    if L < Lh:
                        print("--> The boundary layer is not fully developed")
                        dc = L
                    else:
                        print("--> The boundary layer is developed")
            #ni = fluid.prop()[3]            
            RE = V0 * (dc/fluid.prop()[3])
            delta = 5 * (fluid.prop()[3] * dc/V0)**0.5
            #delta will be updated once the wall shear stress is estimated
        else:
            ans = "None"
            print("... wrong selection!")
    print("\n--> The characteristic lenght scale is {:5.6f} m".format(dc) )
    print("--> The Reynolds number is {:1.5e}".format(RE))
    return dc,RE,FLK,delta
    
def meshH(Vo,ymin,ymax):
    # Levels of refinement
    print("\n---------------- Mesh Design -----------------")
    CN = float(input("* Set the Courant Number (CNF <= 1) = "))
    Lo = float(input("* Set the base mesh level dimension (m) = "))
    dt = Lo*CN / Vo
    print("--> The time discretization dt = {:1.3e} sec".format(dt))
    print("* Set the rate of refinement: ")
    lev = int(input("  1. High;\n  2. Medium;\n  3. Low;\n  [1-2-3]: "))
    
    if lev == 1:
        yref = ymin
        print("... refinement length yref = ymin = {}".format(yref))
    elif lev == 3:
        yref = ymax
        print("... refinement length yref = ymax = {}".format(yref))
    elif lev == 2:
        yref = 0.5*(ymax+ymin)
        print("... refinement length yref = (ymax+ymin)/2 = {}".format(yref))
    else:
        print("... wrong selection!")
        print("... default refinement lenght = ymax")
        yref = ymax
    Lref = round((num.log2(Lo) - num.log2(yref)))
    if (Lref < 0): Lref = 0

    print("--> The refinement level Lref = {}".format(int(Lref)))

     ## Layering of the wall
    Nlrs = 0
    y1 = 999
    Er = float(input(" * Set the espansion ratio Er: [1.1-1.5] "))
    #limit the espantion ratio
    if Er<1.1: Er=1
    elif Er>1.5: Er=1.5
    while y1 > ymin:
        Nlrs = Nlrs + 1
        y1 = yref/((Er)**(Nlrs-1))
        Rs = y1/yref #relative size
    print("--> Number of layers close to the wall: {}".format(Nlrs))
    print("--> The Relative size is {}".format(Rs))
    print("--> First cell size is {}".format(y1))
    print("\n---------------------------------------------")
    print("---------------------------------------------")
    return CN,Lo,dt,Lref,yref,y1,Er,Rs,Nlrs

def calc(fluid,Name,V,Sv): #fluid is an object of the fluid class 
    Cnu = 0.09 #turbulence model constant
    print("\n---------------- Flow Velocity Field -----------------\n")
    V0 = float(input("* Set the free stream velocity [m/s]: "))

    CL,Re,FLK,delta = geom(V0,fluid) #Geometry function calling
    
    if (FLK == "ext") :
        if  Re < 5 * num.power(10, 5):
            print(">>> Streamline Flow [Re < 5*10^5]: Blausius Skin friction coefficient")
            Cf = 0.664 * num.power(Re, -0.5) #Skin-friction coefficient - laminar
            print("--> Cf = 0.664*(Re)^-0.5 = {:4.5f}\n".format(Cf))
        elif (Re >= 5 * num.power(10, 5)):
            print(">>> Turbulent Flow [5x10^5 < Re < 10^9]: Schlichting skin-friction coefficient")
            Cf = num.power(((2 * num.log10(Re))-0.65), -2.3)
            print("--> Cf = (2*log10(Re)-0.65)^(-2.3) = {:4.5f}\n".format(Cf))

    if FLK == "int":
        eps  = float(input("* Set the wall absolute roughness (mm):\n   - Stainless steel [0.0015]\n   - Steel commercial pipe    [0.045 - 0.09]\n   - New cast iron    [0.25 - 0.8]\n   - PVC and Plastic [0.0015 - 0.007]\n   ... Select a value : "))
        eps = eps * 1e-3 # absolute wall roughness
        #Estimation of the friction factor by Colebrook-White formula
        print(">>> Colebrook-White friction coefficient")
        # Iterative procedure
        # First step - Hyp. fully turbulent flow
        turbTerm =  eps/(3.7*CL) #turbulent term
        lambInf = 0.25 * (num.log10(turbTerm)**2)**-1
        lamI = lambInf #First value for the friction coefficient

        errLam = 999
        tol  = 1e-13 #admitted tolerance

        while (errLam > tol):
            lamTerm = 2.51/(Re*(lamI**0.5))
            lamII = 0.25 * (num.log10(turbTerm + lamTerm)**2)**-1
            errLam = num.abs((lamII - lamI)/lamII)
            lamI = lamII
        Cf = lamI
        print("--> Cf = {:4.5f}".format(Cf))

    #Shear stress and shear velocity
    tw = 0.5 * fluid.prop()[0] * Cf * num.power(V0, 2) #Wall shear stress
    Uw = num.power( (tw/fluid.prop()[0]) , 0.5) #Shear velocity
    print("\n>>> Wall mesh treatment")
    yplus_min = float(input("* Set the smallest y+_min =  "))
    yplus_max = yplus_min + 50
    y_min = (yplus_min*fluid.prop()[3])/Uw
    y_max = (yplus_max*fluid.prop()[3])/Uw
    I = 0.16 * num.power(Re,(-1.0/8.0))   #Turbulent intensity (The common choice is I = 0.05)
    K = (3.0/2.0) * num.power((I*V0),2.0) #Turbulent kinetic energy
    u = (2.0/3.0) * num.power(K,0.5)         #Turbulent fluctuation

    #turbulent scale estimation
    #large energy-containing eddies in a turbulent flow.
    if FLK == "ext":
        tls = 0.4*delta
        #delta estimated in the geom function
        #Thermal boundary layer Thickness
        if Re <= 5*10**5:
            delta_t = delta * (fluid.prop()[5])**(1/3) 
            print("--> The Thermal BL height is deltaT = {} m".format(delta_t))
        elif Re > 5*10**5:
            delta_t = delta
            print("--> The Thermal BL height is deltaT = {} m".format(delta_t))

    elif FLK == "int":
        #delta estimated using the shear stress
        delta = 11.6 * (fluid.prop()[3]/Uw) #Turbulent BL thickness
        #Thermal boundary layer Thickness
        if Re <= 5*10**5:
            delta_t = delta * (fluid.prop()[5])**(1/3) 
            print("--> The Thermal BL height is deltaT = {} m".format(delta_t))
        elif Re > 5*10**5:
            delta_t = delta
            print("--> The Thermal BL height is deltaT = {} m".format(delta_t))
        tls = 0.038*CL
        
    print("\n>>> Turbulence free-stream boundary conditions")
    print("--> ymin = {:8.3e} m - wall minimun cell height".format(y_min))
    print("--> ymax = {:8.3e} m - wall maximum cell height".format(y_max))
    print("--> d = {} m - Viscous BL thickness".format(delta))
    print("--> tw = {:8.3e} Pa*m^-2 - Wall shear stress".format(tw))
    print("--> Uw = {:8.3e} m*s^-1 - Shear Velocity".format(Uw))
    print("--> I = {:4.4f} (-) - Turbulent intensity ".format(I))
    print("--> K = {:4.4f} m^2*s^-2 - Turbulent Kinetic Energy".format(K))
    print("--> u' = {:4.4f} m*s^-1 - Turbulent velocity fluctuation".format(u))
    print("\n>>> Epsilon and Omega - Turbulence dissipation")

    epsilon_min = Cnu*num.power(K,1.5)/tls # min turbulence dissipation
    epsilon_max = num.power(Cnu,(3.0/4.0))*num.power(K,1.5)/tls #max turbulence dissipation
    E = (epsilon_min + epsilon_max)/2     #averaged dissipation range
    #Omega is the specific rate of dissipation, of the K
    omega = num.power(Cnu,-(1.0/4.0))*num.power(K,0.5)/tls  #Eq.7.4 in the "OpenFOAM Thecnology Primer" pag. 223
    #Estimated values of turbulence levels
    lver = Cnu*num.power(K,1.5)/E                              #updated value for turbulence length scale
    epsilon  = num.power(Cnu,(3.0/4.0))*num.power(K,1.5)/lver  #updated value for epsilon
    omega    = num.power(Cnu,-(1.0/4.0))*num.power(K,0.5)/lver #updated value for omega
    nut      = fluid.prop()[0]*K/omega #turbulence viscosity
    mi_ratio = nut/fluid.prop()[2]     #viscosity ratio away from inlet patches;

    print("--> Turbulent length scale lver = {:8.5f} m".format(lver))
    print("--> epsilon = {:8.5f} m^2*s^-3 - Turbulent Energy dissipation".format(epsilon))
    print("--> omega = {:8.5f} s^-1 - Specific rate of dissipation".format(omega))
    print("--> nut = {:8.3e} Pa*s - Turbulent viscosity".format(nut))
    print("--> nut/mi = {:4.4f} (-) - Viscosity ratio".format(mi_ratio))

    #Limitation of turbulence levels at inlet section
    if  (mi_ratio > 10):
        mi_ratio_inlet = 10
        epsilon_inlet = (1/mi_ratio_inlet) * Cnu * K * (K/fluid.prop()[3])
        omega_inlet   = (1/mi_ratio_inlet) * (K/fluid.prop()[3])
    elif (mi_ratio <= 10):
        epsilon_inlet = epsilon
        omega_inlet   = omega

    print("--> Inlet Turbulent Energy dissipation epsilon_in = {:5.5f} m^2*s^-3".format(epsilon_inlet))
    print("--> Inlet specific rate of dissipation omega_in = {:5.5f} s^-1".format(omega_inlet))
    print("\n* Wall treatment options:")
    print("    1. 100 < y+ < 300  - Wall Function")
    print("    2.  50 < y+ < 200  - Scalable Wall Function")
    print("    3.       y+ < 11   - Solve the boundary layer")
    WT = True
    while WT:
        WT = int(input("* Set the wall treatment [1,2,3]: "))
        if WT == 1:
            k_wall = K
            omega_wall = 10.0 * (6.0*fluid.prop()[3])/(0.075 * num.power(y_min,2))
            nutw = fluid.prop()[0]*k_wall/omega_wall
            print("\n>>> Wall Function treatment ")
            print("--> K_wall = K = {:8.3f} m^2*s^-2".format(K))
            print("--> Omega_wall = 10*(6*ni)/(0.075*y_min^2) = {:8.3f} s^-1".format(omega_wall))
            print("--> nutw = {:2.3e} Pa*s".format(nutw))
            wt = "Wall treatment: Wall Function"
            WT = False
        elif WT == 2:
            k_wall = 0
            omega_wall_lR = 10*(6.0*fluid.prop()[3])/(0.075*num.power(y_min,2.0))
            omega_wall_log = num.power(K,0.5)/num.power(Cnu,1/4)/0.4/y_min
            omega_wall = num.power((num.power(omega_wall_lR,2)+num.power(omega_wall_log,2)),0.5)
            nutw = fluid.prop()[0]*k_wall/omega_wall
            print("\n>>> Scalable Wall Function treatment")
            print("--> K_wall = 0 m^2*s^-2")
            print("--> Omega_wall = (K^0.5)/(0.4*y_min*Cnu^1/4) = {:8.3f} s^-1".format(omega_wall))
            print("--> nutw = {:2.3e} Pa*s".format(nutw))
            wt = "Wall treatment: Scalable Wall Function"
            WT = False

        elif WT == 3:
            k_wall     = num.power(10,-10)
            omega_wall = num.power(10,-10)
            nutw = num.power(10,-10)
            print("\n>>> Solving the boundary layer")
            print("--> K_wall     = {:8.3f} m^2*s^-2".format(k_wall))
            print("--> Omega_wall = {:8.3f} s^-1".format(omega_wall))
            print("--> nutw = {:2.3e} Pa*s".format(nutw))
            wt = "Wall treatment: Solving the Wall"
            WT = False
        else:
            WT = True
            print("\n ... you did a wrong selection!")

    #Mesh design
    CN,Lo,dt,Lref,yref,y1,Er,Rs,Nlrs = meshH(V0,y_min,y_max)    #Call to meshH function!

    #Writing on BC file
    BCfile = "BL_"+Name+".txt" #Change name of the file
    data_BC = open(BCfile,'w')

    data_BC.write("###########################################################\n")
    data_BC.write("###########################################################\n")
    data_BC.write("#########                                         #########\n")
    data_BC.write("#########       BOUNDARY LAYER ESTIMATION         #########\n")
    data_BC.write("#########                                         #########\n")
    data_BC.write("#########         - Fluid Specifications          #########\n")
    data_BC.write("#########         - Flow Velocity Field           #########\n")
    data_BC.write("#########         - Turbulence characteristics    #########\n")
    data_BC.write("#########         - Mesh design                   #########\n")
    data_BC.write("#########                               ver.{}.{}   #########\n".format(V,Sv))
    data_BC.write("###########################################################\n")
    data_BC.write("###########################################################\n")
    data_BC.write("\n------------------------------------------------------------\n")
    data_BC.write("\n--------- Domain/SubDomain Name: {}\n".format(Name))
    data_BC.write("\n------------------------------------------------------------\n")
    data_BC.write("--------------------- Fluid Specifications -----------------\n")
    data_BC.write("Fluid = {} \n".format(fluid.nameflu))
    data_BC.write("Temperature = {} C \n".format(fluid.T))
    data_BC.write("The density is {:6.3f} Kg/m^3 \n".format(fluid.prop()[0]))
    data_BC.write("The specific volume is {:8.2f} N/m^3 \n".format(fluid.prop()[1]))
    data_BC.write("The dinamic viscosity is {:2.4e} Pa*s \n".format(fluid.prop()[2]))
    data_BC.write("The Kinematic viscosity is {:2.4e} m^2/s \n".format(fluid.prop()[3]))
    data_BC.write("The Thermal diffusivity {:2.4e} m^2/s \n".format(fluid.prop()[4]))
    data_BC.write("Prandtl number is {:2.4e}\n".format(fluid.prop()[5]))
    data_BC.write("------------------------------------------------------------\n\n")

    data_BC.write("------------------------------------------------------------\n")
    data_BC.write("------------------- Flow Velocity Field  -------------------\n")
    data_BC.write("Free stream velocity Vfree is {:5.3f} m/s \n".format(V0))
    data_BC.write("The characteristic lenght scale is {:1.5e} m \n".format(CL))
    data_BC.write("The Reynolds number is {:10.8e}  \n".format(Re))
    data_BC.write("Skin-friction coefficient Cf is {:1.10e}  \n".format(Cf))
    data_BC.write("------------------------------------------------------------\n\n")

    data_BC.write("------------------------------------------------------------\n")
    data_BC.write("------------------- Turbulent characteristics --------------\n")
    data_BC.write("ymin = {:1.5e} m - First cell in boundary layer ymin \n".format(y_min))
    data_BC.write("ymax = {:1.5e} m - Last cell in boundary layer ymax \n".format(y_max))
    data_BC.write("delta = {:1.5e} m - Velocity BL thickness \n".format(delta))
    data_BC.write("deltaT = {:1.5e} m - Thermal BL thickness \n".format(delta_t))

    data_BC.write("tw = {:5.5e} Pa*m^-2 - Wall shear stress\n".format(tw))
    data_BC.write("Uw = {:5.5f} m*s^-1 - Shear Velocity\n".format(Uw))
    data_BC.write("u = {:5.5e} m*s^-1 - Turbulent velocity fluctuation\n".format(u))
    data_BC.write("lver = {:5.5e} m - Turbulent length scale\n\n".format(lver))

    data_BC.write("I = {:1.6f} - Turbulent intensity\n".format(I))
    data_BC.write("K = {:1.5e} m^2*s^-2 - Turbulent Kinetic Energy\n".format(K))
    data_BC.write("epsilon = {:1.5e} m^2*s^-3 - Turbulent Energy dissipation \n".format(epsilon))
    data_BC.write("omega = {:5.5e} s^-1 - Specific rate of dissipation \n\n".format(omega))

    data_BC.write(">>> {st}\n".format(st=wt))
    data_BC.write("nut = {:5.6e} Pa*s - Turbulent viscosity \n".format(nut))
    data_BC.write("nut/mi = {:1.4f} - Viscosity ratio \n".format(mi_ratio))
    data_BC.write("epsilon_in = {:5.6e} m^2*s^-3 - Inlet Turbulent Energy dissipation \n".format(epsilon_inlet))
    data_BC.write("omega_in = {:5.6e} sec^-1 - Inlet specific rate of dissipation \n".format(omega_inlet))
    data_BC.write("Kwall = {:5.6e} m^2*s^-2 - Wall Turbulent Kinetic Energy \n".format(k_wall))
    data_BC.write("Omega_wall = {:5.6e} sec^-1 - Wall Specific rate of dissipation \n".format(omega_wall))
    data_BC.write("nutWall = {:5.6e} Pa*s - Wall Turbulent viscosity \n".format(nutw))
    data_BC.write("------------------------------------------------------------\n\n")

    data_BC.write("------------------------------------------------------------\n")
    data_BC.write("------------------------ Mesh Design  ----------------------\n")
    data_BC.write("CN   = {:1.2f} - Courant Number\n".format(CN))
    data_BC.write("Lo   = {:5.4e} m - first level dimension\n".format(Lo))
    data_BC.write("dt   = {:1.5e} sec - time discretization\n".format(dt))
    data_BC.write("Y+   = {:5.0f} - Selected value of y+\n".format(yplus_min))
    data_BC.write("Lref = {} - Levels of refinement\n".format(int(Lref)))
    data_BC.write("Yref = {:5.6e} m - Last cell height after refinement\n".format(yref))
    data_BC.write("Y1   = {:5.6e} m - First cell over the wall with layers\n".format(y1))
    data_BC.write("Er   = {} - Espantion Ratio\n".format(Er))
    data_BC.write("Rs   = {:1.4f} - Relative Size\n".format(Rs))
    data_BC.write("Nlrs = {} - Number of layers\n".format(int(Nlrs)))
    data_BC.write("------------------------------------------------------------\n")
    data_BC.write("------------------------------------------------------------\n")
    data_BC.close()

    print("\n ... Estimation completed!\n ... Please find a report file in the code directory\n ... Enjoy your simulation!")
####################################################################  END of Functions


######################################################################## MAIN
F,Name = spec(V,Sv)      #fluid specification function -  F is an object of the class physics
calc(F,Name,V,Sv)    #Flow field and Turb. estimation function

input("\n\n >>> Press a key to exit!")
##################################################################### END MAIN

