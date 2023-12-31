# General Information

This directory contains two files: 

1. [*GaussianDiffusion.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/GaussDiffusion/GausschenDiffusion.m)<br>
   -> Script to calculate and compare the numerical solution of the diffusion of an initial Gaussian temperature distribution with its analytical solution. The diffusive part of the energy equation can be solved by different finite difference discretization methods (explicit, implicit, Crank-Nicholson, Alternative direct implicit).<br><br>
   The initial Gaussian temperature anomaly is defined by: <br><br>
   $T=T_0 + A exp(-\frac{(x-0.5L)^2 + (z-0.5H)^2}{2\sigma^2 /\pi})$,<br><br>
   where *A* is the amplitude, *σ* its width (here defined as a percentage of *L*), and T<sub>0</sub> is the background temperature. <br>
   The analytical solution is given by (*reference*?):<br><br>
   $T_{ana} = T_0 + \frac{A}{(1 + 2\pi t \kappa / \sigma^2)} exp(-\frac{(x-0.5L)^2 + (z-0.5H)^2}{2\sigma^2 / \pi + 4t\kappa})$<br>

2. [*ResolutionTest.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/GaussDiffusion/ResolutionTest.m)<br>
   -> Script to run a resolution test for each finite difference discretization scheme (explicit, implicit, CNV, ADI).

------------------------------------------------------------------------

## Example of a gaussian diffusion 

The model is solved with an explicit solver and for a resolution of 41 x 41. <br>
The parameter for the Gaussian temperature anomaly are: <br> 
A = 200 [K], <br>
σ = 10 [% of L], <br>
T<sub>0</sub> = 1000 [K]. <br>

![Evolution](https://github.com/LukasFuchs/FDCSGm/assets/25866942/4ef2feff-8b47-41f8-83df-d8ed04f7f5a3)<br>
**Figure 1.** Transient solution. **Top Left:** Numerical solution of the temperature field (background color and black contour lines), overlain by the same contour lines of the analytical solution (yellow dashed lines) for every second iteration step. **Top Right:** Total difference between the analytical and numerical temperature solution (*T*<sub>ana</sub> - *T*<sub>num</sub>). **Bottom Left:** Temperature profile through the middle of the model domain; solid line - numerical solution; yellow dashed line - analytical solution. **Bottom Right:** Root mean square error with time. 

--------------------------------------------------------------

## Resolution Test

&emsp;For the resolution test, we compared the RMS between the analytical and numerical solution, the maximum, and mean of the temperature of the analytical and numerical solution for each finite differnece scheme after a period of 10 Million years. The resolution is 21, 41, 61, 81, 101, 121, 141, 161, respectively.

![Comparison](https://github.com/LukasFuchs/FDCSGm/assets/25866942/d298959d-8622-462f-b888-c434b57747c8)<br>
**Figure 2.** RMS, *T*<sub>max</sub>, and *T*<sub>mean</sub> for each finite difference discretization scheme. 

