import matplotlib.pyplot as plt
import numpy as np
 
# Given constants
ni = 15e-6  # kinematic viscosity (m^2/s)
L = 1  # characteristic length (m)
U = np.arange(50, 200, 1)  # velocity values (m/s)
RE_values = []  # list to store Reynolds numbers
d_sqrt_values = []  # list to store delta*sqrt(U/(nu*L))

# Loop through velocities
for u in U:
    # Calculate Reynolds number
    RE = u * (L / ni)
    RE_values.append(RE)
    
    if RE < 3 * 10**5:
        # Laminar flow - boundary layer thickness approximation
        delta = 5 * L / np.sqrt(RE)
    else:
        # Turbulent flow - boundary layer thickness approximation
        delta = 0.16 * L*RE**(1/7)
    # Calculate delta * sqrt(U / (nu * L))
    d_sqrt = delta * np.sqrt(u / (ni * L))
    d_sqrt_values.append(d_sqrt)

# Plotting the results
plt.plot(RE_values, d_sqrt_values)
plt.xlabel('Reynolds Number (RE)')
plt.ylabel('delta * sqrt(U / (nu * L))')
plt.title('Boundary Layer Thickness Scaling vs Reynolds Number')
plt.xscale('log')  # Set x-axis to logarithmic scale
plt.show()
