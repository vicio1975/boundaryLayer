#Project object 
This is about a code written in Python for the estimation of the boundary layer for a confined or free turbulent flow and a mesh design according to a set wall treatment.

Here the main steps:
- Fluid Specifications (air or water)
- Flow Velocity Field
- Turbulence characteristics
- Mesh design 

Here some steps detail:
a) Once the kind of flow and the temperature are defined, the flow characteristics are estimated: density, dynamic and kinematic viscosity (Sutherland formula for air dynamic viscosity, Kell formula for water density and Vogel equation for water dynamic density). 
The free stream velocity and the geometry settings (internal or external flow, type of geometry and/or the characteristic length of the flow-Dc) allow to estimate the number of Reynolds (V*Dc/ni).
b) The function "calc" estimates the skin friction factor dependending on the Reynolds number and the kind of flow (external or internal flow): the Blausius formula is used for laminar flow and Schlichting for turbulent flows; for internal flows the Colebrook-White formula is used for all flow regimes. The Shear stress and the shear velocity are estimated using the skin friction factor.
c) The shear stress and shear velocity allow to estimate the Turbulence free-stream boundary conditions, such as: Viscous BL thickness, Turbulent intensity, Turbulent Kinetic Energy, Turbulent velocity fluctuation, Epsilon and Omega - Turbulence dissipation and dissipation frequency, turbulence viscosity and the viscosity ratio. 
d) A Mesh design is proposed according to the selected wall treatment (wall function, scalable wall function and fully resolved wall).

The code writes the results on a text file.

~~~~

Termal boundary layer implementation:
	Prandtl Number = kinematic_viscosity/thermal_diffusivity = [mi/ro] / [Kc/(ro*cp)] = 
	Prandtl Number = (mi*cp)/Kc 
	ni : momentum diffusivity (kinematic viscosity), ni = mi/ro (SI units: m2/s)
	alfa : thermal diffusivity, alfa = Kc/(ro*cp) (SI units: m2/s)
	mi : dynamic viscosity, (SI units: Pa s = N s/m2)
	Kc: thermal conductivity, (SI units: W/m-K)
	cp: specific heat, (SI units: J/kg-K)
	ro : density, (SI units: kg/m3).

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
The values up to 1000 K were originally published in "Tables of Thermal Properties of Gases", NBS Circular 564,1955.
The last five rows were calculated from a formula by B G Kyle "Chemical and Process Thermodynamics", Englewood Cliffs / Prentice Hall, 1984, and have <1% error.

- 2006 ASHRAE Handbook Refrigeration (SI) - Table 2 - Thermal Property Models for Water and Ice (t range: 40 < t < 150°C) - Choi and Okos (1986)

