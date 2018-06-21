#Project object 
This is about a code written in Python for the estimation of the boundary layer for a confined or free turbulent flow and a mesh design according to a set wall treatment.

Here the main steps:
- Fluid Specifications (air and water)
- Flow Velocity Field
- Turbulence characteristics
- Mesh design 

Once the kind of flow and the temperature are defined, the flow characteristics are estimated: density, dynamic and kinematic viscosity (Sutherland formula for air dynamic viscosity, Kell formula for water density and Vogel equation for water dynamic density). 
The free stream velocity and the geometry settings (confined or external flow, type of geometry and/or the characteristic lenght) allow to estimate the number of Reynolds.
The function "calc" estimates the skin friction factor dependending on the Reynolds number and the kind of flow (external or internal flow): for external flow the Blausius formaula is used for laminar flow and Schlichting for turbulent flows; while for external flows the Colebrook-White formula is used for all flow regimes.
Once the Shear stress and shear velocity are estimated, the Turbulence free-stream boundary conditions are calculated: Viscous BL thickness, Turbulent intensity, Turbulent Kinetic Energy, Turbulent velocity fluctuation, Epsilon and Omega - Turbulence dissipation and dissipation frequency, turbulence viscosity and the viscosity ratio. 
A Mesh design is proposed according to the selected wall treatment (wall function, scalable wall function and fully resolved wall).
The code write the results on a text file.


Bibliography

- Menter F. R. 1993, Zonal 2-equation k-omega turbulence models for aerodynamic flows. AIAA journal.
- Lanunder and Spalding 1972, Mathematical Models of Turbulence. Academic Press. 
- Ferziger and Peric, 2002. Computational methods for fluid dynamics. Ed Berlin Springer. 
- Pope 2000, Turbulent Flows. Cambridge University Press.
- Schlichting, hermann and Klauss Gersten , 2001. Boundary Layer Theory. Ed. Berlin: Springer.
- Wilcox 1988, Re-assessment of the scale determining equation for advanced turbulence models. AIAA Journal. Turbulence modeling for CFD. DCW Industries.
- Maric, Hopken and Mooney, 2014. The OpenFOAM Technology Primer. Ed Sourceflux UG.

- About Thermal Conductivity
Title: The Temperature Dependence of the Thermal Conductivity of Air
Authors: Kannuluik, W. G. & Carman, E. H.
Journal: Australian Journal of Scientific Research, Series A: Physical Sciences, vol. 4, p.305          

- Heat Capacity at constant volume
https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
http://thermopedia.com/content/553/ 
The values up to 1000 K were originally published in "Tables of Thermal Properties of Gases", NBS Circular 564,1955. The last five rows were calculated from a formula by B G Kyle "Chemical and Process Thermodynamics", Englewood Cliffs / Prentice Hall, 1984, and have <1% error.
