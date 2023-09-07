The Blankenbach benchmark is an effective benchmark to test and compare a ode for isoviscous and temperature-dependent viscosity thermal convection. Blankenbach et al. (1989) tested different two dimensional thermal convection codes with a broad variety of numerical techniques and reported steady state values for a number of model parameters. 

Convection is studied in a rectangular box of height *H* and width *L*. The kinematic boundary conditions are free slip along all boundaries, a specified temperature on the top (T<sub>top</sub>) and at the bottom (T<sub>bottom</sub>) and thermal insulation along the lateral boundaries. The difference between T<sub>top</sub> and T<sub>bottom</sub> in all experiments is 1000 K. The following formulation for temperature- and depth-dependent viscosity of the mantle is used: 

$\eta=\eta_0 exp(-b \frac{T-T_{top}}{T_{bottom}-T{top}} + c\frac{y}{H})$

where $\eta_0$ is the viscosity at the top of the model and *b* and *c* are coefficients establishing the dependences of viscosity with temperature and depth, respectively. The model constants are listed in the table below. The density depends linearly on temperature (see *equation of state* formulation). 

--------------------------------------------------------------

### **Constants**
Gravitional acceleartion, **g** [m/s]: 10

Model height, **H** [km]: 1000

Model width, **L** [km]: 1000

Temperature at the top, **T<sub>top</sub>** [K]: 273

Temperautre at the bottom, **T<sub>bottom</sub>** [K]: 1273

Thermal conductivity, **k** [W/m/K]: 5

Specific heat capacity, **c<sub>p</sub>** [J/kg]: 1250

Density, **ρ** [kg/m<sup>3</sup>]: 4000

Thermal expansion coefficient, **α** [1/K]:	2.5∙10<sup>-5</sup>

--------------------------------------------------------------

### **The models were solved with the following solver routines and initial conditions:**

Diffusion: explicit<br>
Advection: semi-lag<br>
Viscosity: const<br>
Initial temperature field: block<br>
Resolution (nx x nz): 51 x 51<br>

--------------------------------------------------------------

## Isoviscous Convection

***Ra* = 10<sup>4</sup>**<br>
Reference viscosity [Pa s]: 10<sup>23</sup>
![Field_SS](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a307c72b-d33e-411a-9f8d-4d7849c6b55a)

**Time series**
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7e7863a1-9360-41af-94a9-58f26065bb37)

**Resolution Test**
![ResTest_eta_const](https://github.com/LukasFuchs/FDCSGm/assets/25866942/b1837f13-d1b2-4a8b-882a-1005013cc6bf)

***Ra* = 10<sup>6</sup>**<br>
Reference viscosity [Pa s]: 10<sup>21</sup>
![Field_SS](https://github.com/LukasFuchs/FDCSGm/assets/25866942/24db94dd-0c7d-4566-aac3-49a995cd3ff0)
**Time series**
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/0d87a365-3347-40a6-83d1-150ea76f34cd)

--------------------------------------------------------------

## Temperature Dependent Convection<br>
nx = nz = 101

***Ra* = 10<sup>4</sup>**<br>
Reference viscosity [Pa s]: 10<sup>23</sup><br>
![Field43800](https://github.com/LukasFuchs/FDCSGm/assets/25866942/d0c64608-e208-4ac3-912b-890e939a1644)
**Time Series**<br>
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a8d16cfe-739b-4233-be88-e06bb98a753f)


