# -*- coding: utf-8 -*-
"""
Created by Vincenzo Sammartano
email:  v.sammartano@gmail.com

"""
from tkinter import ttk
import tkinter as tk
from tkinter import messagebox, font
import math

# Functions

def print_results(x, y,z,w):
    """
    This function prints results on the output frame
    - x is the reynolds number
    - y is the BL thickness
    """
    result_text = "Reynolds number (-) = {:4.0f} \n\nBL thickness (mm) = {:2.2f}\n\nFirst cell height (mm) = {:2.2f}\n\nCells numbers = {:2.0f}".format(x, y,z,w)
    # this deletes the current text
    results_values.delete('1.0', tk.END)
    # this updates the results
    results_values.insert('1.0', result_text)

def ACTF(for1, for2):
    """
    Toggles between free stream and confined flow inputs.
    """
    if modForm.get() == 1:  # Free Stream selected
        for1.config(state='normal', fg='black')
        for2.config(state='normal', fg='gray')
    elif modForm.get() == 2:  # Confined Flow selected
        for2.config(state='normal', fg='black')
        for1.config(state='normal', fg='gray')

def fluid(x, y):
    """
    Estimates the fluid properties based on fluid name and temperature.
    - x: fluid name (string)
    - y: temperature in Kelvin
    Returns: List of fluid properties
    """
    Fp = []
    t = y - 273.15  # Convert temperature to Celsius

    if x == "Air":
        # Air properties
        rot = pAtm / (Rf * y)  # Density
        gamma_t = rot * g  # Specific weight
        ba = 1.458 * 10**(-6)
        sa = 110.4  # Sutherland constant
        mi = ba * (y**1.5) / (y + sa)  # Dynamic viscosity
        ni = mi / rot  # Kinematic viscosity
        Fp = [rot, gamma_t, mi, ni]

    elif x == "Water":
        # Water properties
        ro = 999.83952
        At = 1 + (16.879850 * 10**-3) * t
        rot = (ro + (16.945176 * t) - (7.987040 * 10**-3) * (t**2.0) -
               (46.17046 * 10**-6) * (t**3.0) + (105.56302 * 10**-9) * (t**4.0) -
               (280.54253 * 10**-12) * (t**5.0)) / At
        gamma_t = rot * g
        mi = 2.414 * (10**-5) * 10**(247.8 / (y - 140))  # Dynamic viscosity
        ni = mi / rot  # Kinematic viscosity
        Fp = [rot, gamma_t, mi, ni]

    return Fp

def BL(velocity, length_mm, kinematic_viscosity, flow_type):
    """
    Calculates the boundary layer thickness (delta) based on Reynolds number and flow type.
    """
    length_m = length_mm / 1000.0  # Convert length to meters
    reynolds_number = (velocity * length_m) / kinematic_viscosity
    
    if flow_type == 1:  # Free-stream boundary layer
        if reynolds_number < 3e5:  # Laminar flow
            delta = length_m * ((0.375 / (reynolds_number ** 0.5)) ** (3 / 4))
        else:  # Turbulent flow
            delta = 0.16 * length_m * (reynolds_number ** (-1 / 7))

    elif flow_type == 2:  # Confined flows
        if reynolds_number <= 2e3:
            G = 1.47
        elif 2e3 < reynolds_number <= 2e4:
            G = 1.38
        elif reynolds_number > 2e4:
            G = 1.33
        delta = (122 * math.log(reynolds_number)) / (reynolds_number * G)
    else:
        raise ValueError("Unsupported flow type. Use 1 for free-stream or 2 for confined flows.")

    return delta, reynolds_number

def calculation():
    """
    Define the main calculation workflow
    """
    try:
        # Retrieve inputs
        name_fluid = flu_set.get()
        temp_K = float(t1.get())
        velocity = float(v1.get())
        length = float(ls1.get())
        flow_boundary = modForm.get()
        yplus = int(yplus1.get())
        layer_ratio = float(layer_ratio1.get())
        
        if flow_boundary not in [1, 2]:
            raise ValueError("Flow type must be 1 (free stream) or 2 (confined flow).")

        # Estimate fluid properties
        fluid_properties = fluid(name_fluid, temp_K)
        if not fluid_properties:
            messagebox.showerror("Error", "Unsupported fluid type.")
            return

        density_t, gamma_t, mi_t, ni_t = fluid_properties

        # Calculate boundary layer thickness
        bl_thickness, reynolds_number = BL(velocity, length, ni_t, flow_boundary)

        #Mesh cell size over walls
        #friction factor
        Cf = 0.058*reynolds_number**-0.2
        #wall shear stress
        t_star = 0.5 * density_t * velocity**2 * Cf
        #shear velocity
        u_star = (t_star/density_t)**0.5
        #minimum distance from wall y+  = u_star * y / ni
        y = (yplus/u_star)*ni_t
        #Number of layers needed
        Nl = bl_thickness/y
        #layer ratio
        ytot = y
        layers = 1

        while ytot<=bl_thickness:
            layers = layers + 1
            ytot = ytot*(layer_ratio +1)
            
        #bl_thickness and y in mm
        bl_thickness = bl_thickness*1000
        y = y *1000
        # Display results
        print_results(reynolds_number, bl_thickness, y, layers)

    except ValueError as e:
        messagebox.showerror("Input Error", str(e))
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {e}")

def ex():
    root.destroy()

if __name__ == "__main__":
    # Constants
    Rf = 287.058
    pAtm = 101325
    g = 9.806

    # Tkinter UI Setup
    root = tk.Tk()
    root.title("Boundary Layer Thickness Calculator")
    root.geometry("400x500")
    root.resizable(width=False, height=False)
    
    flu_opt = ["Air", "Water"]
    modForm = tk.IntVar()

    fr1 = tk.LabelFrame(root, text="Inputs")
    fr1.pack(fill="both", padx=10, pady=10)

    fr2 = tk.LabelFrame(root, text="Results")
    fr2.pack(fill="both", padx=10, pady=10)

    tk.Label(fr1, text="Fluid:").grid(row=0, column=0, padx=10, pady=4)
    flu_set = ttk.Combobox(fr1, values=flu_opt,width=10)
    flu_set.set("Air")
    flu_set.grid(row=0, column=1, padx=10)

    tk.Label(fr1, text="Temperature (K):").grid(row=1, column=0, padx=10, pady=5)
    t1 = tk.Entry(fr1)
    t1.insert(0, "300")
    t1.grid(row=1, column=1, padx=10)

    tk.Label(fr1, text="Velocity (m/s):").grid(row=2, column=0, padx=10, pady=5)
    v1 = tk.Entry(fr1)
    v1.insert(0, "1.0")
    v1.grid(row=2, column=1, padx=10)

    tk.Label(fr1, text="Length Scale (mm):").grid(row=3, column=0, padx=10, pady=5)
    ls1 = tk.Entry(fr1)
    ls1.insert(0, "100.0")
    ls1.grid(row=3, column=1, padx=10)

    tk.Label(fr1, text="y+ (-):").grid(row=4, column=0, padx=10, pady=5)
    yplus1 = tk.Entry(fr1)
    yplus1.insert(0, "1")
    yplus1.grid(row=4, column=1, padx=10)

    tk.Label(fr1, text="layer ratio:").grid(row=5, column=0, padx=10, pady=5)
    layer_ratio1 = tk.Entry(fr1)
    layer_ratio1.insert(0, "1.2")
    layer_ratio1.grid(row=5, column=1, padx=10)
    
    for1 = tk.Radiobutton(fr1, text="Free Stream", variable=modForm, value=1, 
                          command=lambda: ACTF(for1, for2))
    for1.grid(row=6, column=0, pady=5)

    for2 = tk.Radiobutton(fr1, text="Confined Flow", variable=modForm, value=2, 
                          command=lambda: ACTF(for1, for2))
    for2.grid(row=6, column=1, pady=5)

    results_values = tk.Text(fr2, height=7)
    results_values.pack(padx=10, pady=10)

    btn_calculate = tk.Button(root, text="Calculate", command=calculation)
    btn_calculate.pack(pady=2)

    btn_exit = tk.Button(root, text="Exit", command=ex)
    btn_exit.pack(pady=2)

    root.mainloop()
