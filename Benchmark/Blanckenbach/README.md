The Blankenbach benchmark is an effective benchmark to test and compare a code for isoviscous and temperature-dependent viscosity thermal convection. Blankenbach et al. (1989) tested different two dimensional thermal convection codes with a broad variety of numerical techniques and reported steady-state values for a number of model parameters. 

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

The temperature and velocity patterns are resolved well and the steady state values for the Nusselt number and the root mean square velocity are nearly reached, too. Variations from the benchmark values might be due to a rather coarse resolution (51x51) and the given solving method, that is the explicit diffusion solver, the semi-lagrangian advection method, and the operator splitting method. 

This becomes even more obvious for higher Rayleigh number (*Ra* = 10<sup>6</sup>) calculations for which the deviations of the Nusselt number are quite significant as well as for a temperature-dependend thermal convection. For the latter, a higher resolution (101x101) is necessary to obtain a stable solution of the temperature field. A higher resolution significantly increases the computation time (at least to reach the same final time as in the previous calculations; steady state seems to be reached rather early within the computations), which should be compansated using a variable grid size, which needs to be implemented first. 

... Resolution test for higher Ra ...

***Ra* = 10<sup>4</sup>**<br>
Reference viscosity [Pa s]: 10<sup>23</sup>
![Field3150](https://github.com/LukasFuchs/FDCSGm/assets/25866942/115d5bb6-3c3d-44a6-992a-b6606d7d5144)
**Figure 1.** Steady-state temperautre (top) and velcoity field (bottom). The arrows whos the velocity vectors at certain grid points.

**Time series**
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7e7863a1-9360-41af-94a9-58f26065bb37)
**Figure 2.** Time series of Nusselt number and root mean square velocity as well as the vertical temperautre profile in the middle of the model domain and the temperautre differences at the corners of the model (in the following order: top left (1), top right (2), bottom right (3), and bottom left (4)).

***Evolutuion***<br>
![Evolution_small](https://github.com/LukasFuchs/FDCSGm/assets/25866942/fb4f6e36-29d7-4f2a-affa-3bca238ae59d)

**Resolution Test**<br>
**Figure 3.** Resolution test for high Rayleigh number calculations. 

***Ra* = 10<sup>6</sup>**<br>
Reference viscosity [Pa s]: 10<sup>21</sup>
![Field_SS](https://github.com/LukasFuchs/FDCSGm/assets/25866942/24db94dd-0c7d-4566-aac3-49a995cd3ff0)
**Figure 4.** Steady-state temperautre and velcoity field. For more details see captions Figure 1. 

**Time series**
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/0d87a365-3347-40a6-83d1-150ea76f34cd)
**Figure 5.** Time series, temperature profile, and corner heat fluxes. For more details see captions of Figure 2.

--------------------------------------------------------------

## Temperature Dependent Thermal Convection<br>
nx = nz = 101

***Ra* = 10<sup>4</sup>**<br>
Reference viscosity [Pa s]: 10<sup>23</sup><br>
![Field43800](https://github.com/LukasFuchs/FDCSGm/assets/25866942/d0c64608-e208-4ac3-912b-890e939a1644)
**Figure 6.**  Steady-state temperautre and velcoity field. For more details see captions Figure 1. 

**Time Series**<br>
![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a8d16cfe-739b-4233-be88-e06bb98a753f)
**Figure 7.** Time series, temperature profile, and corner heat fluxes. For more details see captions of Figure 2.


