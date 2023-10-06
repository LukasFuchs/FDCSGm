# General information 

This directory contains two files: 

1. [*BlanckenbachBenchmark.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/Blanckenbach/BlanckenbachBenchmark.m)<br>
     -> Main script to calculate one benchmark model (define *Ra* and viscosity). 
   
2. [*Blanckenbach_Res_Test.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/Blanckenbach/Blanckenbach_Res_Test.m)<br>
     -> Script to run a resolution test for a given *Ra* and viscosity. See below for more details. 

-------------------------------------------------------------

&emsp;The Blankenbach benchmark is an effective benchmark to test and compare a code for isoviscous and temperature-dependent viscosity thermal convection. Blankenbach et al. (1989) tested different two dimensional thermal convection codes with a broad variety of numerical techniques and reported steady-state values for a number of model parameters. 

&emsp;Convection is studied in a rectangular box of height *H* and width *L*. The kinematic boundary conditions are free slip along all boundaries, a specified temperature on the top (*T<sub>top</sub>*) and at the bottom (*T<sub>bottom</sub>*) and thermal insulation along the lateral boundaries. The difference between *T<sub>top</sub>* and *T<sub>bottom</sub>* in all experiments is 1000 K. The following formulation for temperature- and depth-dependent viscosity of the mantle is used: 

$\eta=\eta_0 exp(-b \frac{T-T_{top}}{T_{bottom}-T{top}} + c\frac{y}{H})$, &emsp; &emsp; &emsp; (1)

where $\eta_0$ is the viscosity at the top of the model and *b* and *c* are coefficients establishing the dependences of viscosity with temperature and depth, respectively. The model constants are listed in the table below. The density depends linearly on temperature (see [*equation of state* formulation](https://github.com/LukasFuchs/FDCSGm/tree/main/StokesProblem#equation-of-state)). 

--------------------------------------------------------------

## **Constants**
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

## **The models were solved with the following solver routines and initial conditions:**

Diffusion: explicit<br>
Advection: semi-lag<br>
Viscosity: const<br>
Initial temperature field: block<br>
Resolution (nx x nz): 51 x 51<br>

--------------------------------------------------------------

## Isoviscous Convection

### Low Rayleigh number cases

&emsp;The temperature and velocity patterns are resolved well (Figure 1) and the steady state values for the Nusselt number and the root mean square velocity are nearly reached (Figure 2), too. Variations from the benchmark values might be due to a rather coarse resolution (51x51) and the given solving methods, that is the explicit diffusion solver, the semi-lagrangian advection method, and the operator splitting method. 

![Field3150](https://github.com/LukasFuchs/FDCSGm/assets/25866942/0c34c0a8-2e05-4ead-9c6d-a85a92b1a3b3)<br>
**Figure 1.** Steady-state temperautre (top) and velcoity field (bottom). The arrows whos the velocity vectors at certain grid points. The reference viscosity is 10<sup>23</sup> [Pa s] and *Ra* = 10<sup>4</sup>.

![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7e7863a1-9360-41af-94a9-58f26065bb37)<br>
**Figure 2.** Time series of Nusselt number and root mean square velocity as well as the vertical temperautre profile in the middle of the model domain and the temperautre differences at the corners of the model (in the following order: top left (1), top right (2), bottom right (3), and bottom left (4)).

![Evolution_small](https://github.com/LukasFuchs/FDCSGm/assets/25866942/fb4f6e36-29d7-4f2a-affa-3bca238ae59d)<br>

### High Rayleigh number cases

&emsp;This becomes even more obvious for higher Rayleigh number (*Ra* = 10<sup>6</sup>) calculations for which the deviations of the Nusselt number are quite significant as well as for a temperature-dependent thermal convection (Figures 3 & 4).

![Field_SS](https://github.com/LukasFuchs/FDCSGm/assets/25866942/24db94dd-0c7d-4566-aac3-49a995cd3ff0)<br>
**Figure 3.** Steady-state temperautre and velcoity field. For more details see captions Figure 1. The reference viscosity is 10<sup>21</sup> [Pa s] and *Ra* = 10<sup>6</sup>.

![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/0d87a365-3347-40a6-83d1-150ea76f34cd)<br>
**Figure 4.** Time series, temperature profile, and corner heat fluxes. For more details see captions of Figure 2.

![Evolution_small](https://github.com/LukasFuchs/FDCSGm/assets/25866942/2cf47636-250b-4494-9f8c-c27fb24aac47)

***Resolution test for higher Ra***

Depdending on the Rayleigh number, we can estimate the number of grid points in the upper thermal boundary layer (*d*) and, thus, the total resolution of our model, with:

$(\frac{H}{d})^3 = \frac{1}{4}Ra$, &emsp; &emsp; &emsp; (2)

Assuming that we want to use *n* grid points within the upper thermal boundary layer, the total number of vertical grid points is given by : 

$nz = (n-1)\sqrt[3]{\frac{Ra}{4}}+1$, &emsp; &emsp; &emsp; (3)

![ResTest_eta_constRa1e6](https://github.com/LukasFuchs/FDCSGm/assets/25866942/29b92c83-ddb2-498b-a0fb-c3ca3da6f24a)<br>
**Figure 5.** Resolution test for high Rayleigh number calculations. Needs some further exploration, why the rms is not propergating towards the benchmark value!

--------------------------------------------------------------

## Temperature-Dependent Thermal Convection

nx = 101<br>
nz = 101<br>
a = log(1000)<br>
b = 0

&emsp;For a temperature-dependent thermal convection, a higher resolution (101x101) is necessary to even obtain a stable solution of the temperature field (Figures 6 & 7). A higher resolution significantly increases the computation time (at least to reach the same final time as in the previous calculations; steady state seems to be reached rather early within the computations), which should be compansated using a variable grid size, which needs to be implemented first. 

![Field43800](https://github.com/LukasFuchs/FDCSGm/assets/25866942/d0c64608-e208-4ac3-912b-890e939a1644)<br>
**Figure 6.**  Steady-state temperautre and velcoity field. For more details see captions Figure 1. The reference viscosity is 10<sup>23</sup><br> [Pa s] and *Ra* = 10<sup>4</sup>.

![TimeSeries](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a8d16cfe-739b-4233-be88-e06bb98a753f)<br>
**Figure 7.** Time series, temperature profile, and corner heat fluxes. For more details see captions of Figure 2.

***Evolution***<br>
![Evolution_small](https://github.com/LukasFuchs/FDCSGm/assets/25866942/6dd49bc4-258c-4334-8d90-513984750067)<br>
