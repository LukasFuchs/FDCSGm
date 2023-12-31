# General Information

&emsp; This benchmark enables to compare the numerical velocity solution of a purely dynamically driven system (i.e., here only by the density) with its analytical solution (e.g., Ramber, 1968). The model consists of two layers of certain thicknes (*h<sub>1</sub>*, *h<sub>2</sub>*), density (*ρ<sub>1</sub>*, *ρ<sub>2</sub>*) and viscosity (*η<sub>1</sub>*, *η<sub>2</sub>*), where its boundary is perturbed by a certain sinisoidal perturbation defined as: 

$\Delta A cos(2 \pi \frac{x-0.5L}{\lambda}),$&emsp;&emsp;&emsp; (1)

where $\Delta A$ is the amplitude, *x* the horizontal coordinate, *L* the length of the model, and $\lambda$ the wavelength of the perturbation. I use the marker field to define the area of the different layers and interpolate the density and viscosity from the markers back to the regular grid. 

&emsp;I follow the description of *Gerya (2009)* to calculate the vertical velocity of the centered amplitude and its analytical solution. 

-------------------------------------------------------------------------------

# Rayleigh-Taylor Instability (RTI) Growth Rate

![Setup](https://github.com/LukasFuchs/FDCSGm/assets/25866942/107fd9e6-e00a-4fc8-bc33-396b9cd85a5b)<br>
**Figure 1.** Rayleigh-Taylor Instability. Initial setup and solution of the first time step for a model with a wavelength of 500 m and an amplitude of 100 m. **LEFT:** Density field overlain by the velocity vectors on the regular grid. **RIGHT:** Viscosity field. 

![SetupZoom](https://github.com/LukasFuchs/FDCSGm/assets/25866942/8a51eb05-920a-4f92-85b9-c6aa9f059afb)<br>
**Figure 2.** Zoom in on the area of interest for calculating the vertical velocity of the central wave in the perturbation. The thin solid line is the compositional boundary. The white dots are the tracers within the lower layer. The red dots are the nodes of the vertical velcoity grid surrounding top point of the central wave, which are used to interpolate the vertical velocity from the staggered grid to the point of interest. The arrows show the velocity field of the instantaneus solution of the inital condition. 

&emsp; Below, the figure shows the results of a series of RTIs of certain viscosity and density contrast for a range of amplitudes and perturbation wavelenghts (for more details see figure caption). The code does calculate the vertical velocity for the centered amplitude reasonably well, however, there are some deviations from the analytical solution and the code struggels for very small perturbation (i.e., A = h<sub>1</sub>/1500; not shown here for visualization reasons). These are topics remaining to analyze in detail regarding the accuracy of the code. 

![GrowthRate_15_150_51](https://github.com/LukasFuchs/FDCSGm/assets/25866942/1dfa3f04-03f3-43ca-8e81-d6ba6c56c8b1)<br>
**Figure 3.** Growth rate of a rayleigh-taylor instability. The solid lines are the analytical solution of the growth rate for a certain perturbation amplitude, wavelength, density and viscosity ratio (for more details see Gerya, 2009, Figure. 16.1). The circles are the numerical solution for an amplitude of $\Delta A = h_1/15$ and the dots for an amplitude of $\Delta A = h_1/150$. The resolution of the numerical grid is defined as 101 in vertical direction and the number of grid points varies with the horizontal direction such that the gird is always squared, i.e., dx = dz. The number of marker varies as well, where I assume 25 marker per each finite difference cell.

# Directory Content 
[GrowthRate.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/RTI/GrowthRate.m)<br>
&emsp;-> Script to calculate the series of models for the growth rate shown in **Figure 2**. 

[RTI_semiLag.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/RTI/RTI_semiLag.m)<br>
&emsp; -> Script to calculate a time-dependent RTI using the semi-lagrangian advection method for composition. I simply wrote this script to check if advection of composition is possible with semi-lagrangian, too. It works to a certain extent, however, caution needs to be taken along the compositional interface. 

[RTI_tracers.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/RTI/RTI_tracers.m)<br>
&emsp; -> Script to calculate a time-dependent RTI using the tracer advection method. 

# References 

*Ramberg, H. (1968). Instability of layered systems in the field of gravity. I. Physics of the Earth and Planetary Interiors, 1(7), 427-447.*
