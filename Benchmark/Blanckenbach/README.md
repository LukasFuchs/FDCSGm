... Additional info missing ...

--------------------------------------------------------------

### **The models were solved with the following solver routines and initial conditions:**

Diffusion: explicit

Advection: semi-lag

Viscosity: const

Initial temperature field: block

Resolution (nx x nz): 51 x 51

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

The kinematic boundary conditions are free-slip on all boundaries and the thermal boundary conditions are constant temperature at the top and bottomt and no lateral heat flux at the sides. 

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
The temperature dependent viscosity is given by: 

$\eta=\eta_0 exp(-b \frac{T-T_{top}}{T_{bottom}-T{top}} + c\frac{y}{H})$

where $\eta_0$ is the viscosity at the top of the model and *b* and *c* are coefficients establishing the dependences of viscosity with temperature and depth, respectively. 

***Ra* = 10<sup>4</sup>**<br>
Reference viscosity [Pa s]: 10<sup>23</sup>

