# -*- coding: utf-8 -*-
"""
Created by Vincenzo Sammartano
email:  v.sammartano@gmail.com

"""
from tkinter import ttk
import tkinter as tk
from tkinter import messagebox, filedialog, font
import os
from numpy import log as ln

# Functions

def print_results(x, y):
    """
    This function prints results on the output frame
    - x is the reynolds number
    - y is the BL thickness
    """
    result_text = "Reynolds number = {one}".format(one=x) + "\nBL thickness = {two}".format(two=y)
    print(result_text)
    # this deletes the current text
    results_values.delete('1.0', tk.END)
    # this updates the results
    results_values.insert('1.0', result_text)

def ACTF():
    if modForm.get() == 1:
        for1.config(state='normal')
        for1.config(fg='black')
        for2.config(fg='gray')
    
    if modForm.get() == 2:
        for2.config(state='normal') 
        for1.config(fg='gray')
        for2.config(fg='black')

def inputs():
    """
    This starts the calculation
    The inputs is made by: name of fluid, temperature, velocity, and Length scale
    """
    return [flu_set.get(), float(t1.get()), float(v1.get()), float(ls1.get()), modForm.get()]

def fluid(x, y):
    """
    estimate the fluid properties:
    it is getting the fluid name, x
    and the temperature of the fluid in K, y
    """
    Fp = []
    # Temperature in Celsius
    t = y - 273.15

    if x == "Air":
        # Air
        rot = pAtm / (Rf * y)  # Density as function of temperature in Kelvin [Kg/mc]
        gamma_t = rot * g  # specific weight at t°C
        # Sutherland Equation
        ba = 1.458 * 10**(-6)
        sa = 110.4  # Kelvin
        mi = ba * (y**1.5) / (y + sa)  # Dinamic Viscosity  Pa s = kg m^-1 s^-1
        ni = mi / rot  # Cinematic Viscosity  2·s-1
        FpA = [rot, gamma_t, mi, ni]
        Fp = FpA

    if x == "Water":
        # Water - Kell formulation
        ro = 999.83952
        At = 1 + (16.879850 * 10**-3) * t  # constant to be used in the Kell formula
        # rot is the water density as function of temperature [Kg/e]
        rot = (ro + (16.945176 * t) - (7.987040 * 10**-3) * (t**2.0) -
               (46.17046 * 10**-6) * (t**3.0) + (105.56302 * 10**-9) * (t**4.0) -
               (280.54253 * 10**-12) * (t**5.0)) / At
        gamma_t = rot * g  # specific weight at t°C
        # 0<t<370 Celsius
        mi = 2.414 * (10**-5) * 10**(247.8 / (y - 140))  # Dinamic Viscosity  Pa s = kg m^-1 s^-1
        ni = mi / rot  # Kinetic Viscosity  m2·s-1
        FpW = [rot, gamma_t, mi, ni]
        Fp = FpW
    return Fp

def BL(x, y, z, w):
    """
    This function works on the BL thickness estimation. Here x and y are:
    - x is the velocity m·s-1
    - y is the kinematic viscosity m2·s-1
    - z is the length scale (mm)
    - w is the flow boundary type.
    """
    # z is scaled in meter
    z = z / 1000

    #### Free stream
    # Boundary layer thickness estimation - delta
    if w == 1:
        # Reynolds number
        RE = x * (z / y)
        if RE < 5 * 10**5:
            # Hansen(1928) approach https://ntrs.nasa.gov/api/citations/19930094831/downloads/19930094831.pdf
            # https://ntrs.nasa.gov/api/citations/20150019333/downloads/20150019333.pdf
            d1 = 4.91
            d2 = 1 / 2
            delta = d1 * z / (RE**d2)
        elif RE >= 5 * 10**5:
            d1 = 0.37
            d2 = 1 / 5
            delta = d1 * z / (RE**d2)

    #### Confined Flows
    elif w == 2:
        # Reynolds number
        RE = x * (z / y)
        # https://elmoukrie.com/wp-content/uploads/2022/04/hermann-schlichting-deceased-klaus-gersten-boundary-layer-theory-springer-verlag-berlin-heidelberg-2017.pdf
        if RE <= 2e3:
            G = 1.47
        elif RE > 2e3 and RE <= 2e4:
            G = 1.38
        elif RE > 2e4:
            G = 1.33
        
        delta = 122 * ln(RE) / (RE * G)

    return RE, delta

def calculation():
    """
    Define the main calculation workflow
    """
    # Reading the inputs from the UI entry boxes:
    name_fluid = inputs()[0]
    temp_K = inputs()[1]
    velocity = inputs()[2]
    length = inputs()[3]
    flow_boundary = inputs()[4]

    # Estimating the fluid properties:
    density_t = fluid(name_fluid, temp_K)[0]
    gamma_t = fluid(name_fluid, temp_K)[1]
    mi_t = fluid(name_fluid, temp_K)[2]
    ni_t = fluid(name_fluid, temp_K)[3]

    # Boundary Layer estimation
    Reynolds, bl_thickness = BL(velocity, ni_t, length, flow_boundary)
    # Plot results in the box
    print_results(Reynolds, bl_thickness)

def ex():
    """
    This destroys the UI
    """
    root.destroy()
    
def resize(event):
    new_width = event.width
    new_height = event.height




if __name__ == "__main__":
    ###### UI Block    
    ## Tkinter Window
    root = tk.Tk()
    root.geometry("390x450")
    root.title("Boundary Layer Thickness")
    root.resizable(width=False, height=False)
    root.bind("<Configure>", resize)
    # Fonts
    f_H12B = font.Font(family='Helvetica', size=12, weight='bold')
    f_H12 = font.Font(family='Helvetica', size=12, weight='normal')
    f_H10 = font.Font(family='Helvetica', size=10, weight='bold')

    # Constants
    Rf = 287.058  # Universal Constant of Gases [J/(Kg K)]
    pAtm = 101325  # [Pa] atmospheric pressure
    g = 9.806  # [m/s2] gravitational accelaration

    # List of fluids
    flu_opt = ['Air', 'Water']

    ## UI Objects
    Labels_Names = ["Fluid = ", "Temperature (K) = ", "FreeStream Velocity (m/s) =", "Length Scale (mm) = "]

    ### Frames
    Frames_Names = ["Input", "Output", ""]
    Frames = []  # preallocation

    # input frame 
    fr1 = tk.LabelFrame(root, text=Frames_Names[0])
    fr1.grid(row=0, column=0, padx=12, pady=12, sticky='ew')
    fr1.config(borderwidth=2, font=f_H12B)

    # output frame 
    fr2 = tk.LabelFrame(root, text=Frames_Names[1])
    fr2.grid(row=1, column=0, padx=12, pady=12, sticky='ew')
    fr2.config(borderwidth=2, font=f_H12B)

    # Buttons frame 
    fr3 = tk.Frame(root)
    fr3.grid(row=2, column=0, padx=12, pady=12, sticky='ew')

    Frames.append(fr1)
    Frames.append(fr2)
    Frames.append(fr3)

    # Buttons
    Buttons_Names = ["Estimate", "Exit"]
    Buttons = []  # preallocation

    for i, B in enumerate(Buttons_Names):
        Buttons.append(tk.Button(Frames[2], text=B, font=f_H12B))
        Buttons[i].grid(row=0, column=i, sticky='ew', padx=20)

    Buttons[0].config(command=calculation)
    Buttons[1].config(command=ex)

    ### Labels
    Labels = []
    for i, L in enumerate(Labels_Names):
        Labels.append(tk.Label(Frames[0], text=L, font=f_H12))
        Labels[i].grid(row=i, column=0, sticky='E', ipadx=10, pady=5)

    # Entries
    # Fluid selection
    flu_set = tk.StringVar()
    f1 = tk.OptionMenu(Frames[0], flu_set, *flu_opt)
    f1.config(width=8, font=f_H10)
    f1.grid(row=0, column=1, sticky='e', padx=20)
    flu_set.set(flu_opt[0])

    # temperature selection
    t1 = tk.Entry(Frames[0], justify="center", font=f_H12)
    t1.config(width=10)
    t1.grid(row=1, column=1, pady=5)
    t1.insert("end", 290)

    # Velocity selection
    v1 = tk.Entry(Frames[0], justify="center", font=f_H12)
    v1.config(width=10)
    v1.grid(row=2, column=1, pady=5)
    v1.insert("end", 1)

    # Length scale selection
    ls1 = tk.Entry(Frames[0], justify="center", font=f_H12)
    ls1.config(width=10)
    ls1.grid(row=3, column=1, pady=5)
    ls1.insert("end", 1)

    # Boundary inputs
    modForm = tk.IntVar()
    for1 = tk.Radiobutton(Frames[0], text="free stream", padx=10, variable=modForm, value=1, font=f_H12)
    for1.grid(row=4, column=0, sticky='ew')
    for1.configure(command=ACTF, indicatoron=1)
    for1.select()
    for1.configure(state='normal')

    for2 = tk.Radiobutton(Frames[0], text="confined flow", padx=10, variable=modForm, value=2, font=f_H12)
    for2.grid(row=4, column=1, sticky='ew')
    for2.configure(command=ACTF)
    for2.config(fg='gray')

    # Add a label to display the results
    results_values = tk.Text(fr2, font=f_H12B, width=30, height=5, padx=2, pady=2, bd=2, wrap='word')
    results_values.grid(row=0, column=0, padx=2, pady=20, sticky="ew")
    s = "Reynolds number = \n" + "BL thickness = "
    results_values.insert('1.0', s)

    # Make the widgets expandable
    root.grid_columnconfigure(0, weight=1)
    fr1.grid_columnconfigure(1, weight=1)
    fr2.grid_columnconfigure(0, weight=1)
    fr3.grid_columnconfigure(0, weight=1)
    fr3.grid_columnconfigure(1, weight=1)

    root.mainloop()
