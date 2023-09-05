### **The models were solved with the following solver routines and initial conditions:**

Diffusion: explicit

Advection: semi-lag

Viscosity: const

Initial temperature field: block

Resolution (nx x nz): 51 x 51

Reference viscosity [Pa s]: 1.00e+23

### **Constants**
Gravitional acceleartion, **g** [m/s]: 10

Model height, **H** [km]: 1000

Model width, **L** [km]: 1000

Temperature at the top, **T<sub>top</sub>** [K]: 273

Temperautre at the bottom, **T<sub>bottom</sub>** [K]: 1273

Thermal conductivity, **k** [W/m/K]: 5

Specific heat capacity, **c<sub>p</sub>** [J/kg]: 1250

Density, **ρ** [kg/m<sup>3</sup>]: 4000

Thermal expansion coefficient, **alpha** [1/K]:	2.5e-5

The kinematic boundary conditions are free-slip on all boundaries and the thermal boundary conditions are constant temperature at the top and bottomt and no lateral heat flux at the sides. 

## Isoviscous Convection

Ra = 1e4

![Field_SS](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a307c72b-d33e-411a-9f8d-4d7849c6b55a)

*Time series*

![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7e7863a1-9360-41af-94a9-58f26065bb37)

Ra = 1e6



## Temperature Dependent Convection
