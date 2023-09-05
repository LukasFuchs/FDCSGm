**The models were solved with the following solvers and initial conditions:**

Diffusion: explicit
Advection: semi-lag
Viscosity: const
Initial temperature field: block
Resolution (nx x nz): 51 x 51
Reference viscosity [Pa s]: 1.00e+23

## **Constants**
Gravitional acceleartion, **g** [m/s]: 10

Model height, **H** [km]: 1000

Model width, **L** [km]: 1000

Temperature at the top, **T<sub>top</sub>** [K]: 273

Temperautre at the bottom, **T<sub>bottom</sub>** [K]: 1273

Thermal conductivity, **k** [W/m/K]: 5

Specific heat capacity, **c<sub>p</sub>** [J/kg]: 1250

Density, **ρ** [kg/m<sup>3</sup>]: 4000

Thermal expansion coefficient, **alpha** [1/K]:	2.5e-5

Für die dynamischen Randbedingungen nehmen wir überall free-slip an. Die thermischen 
Randbedingungen sind konstante Temperatur oben und unten und keinen Waermefluss an den 
lateralen Rändern. 

## Isoviscous Convection

## Temperature Dependent Convection
