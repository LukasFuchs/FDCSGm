# General Information

This directory contains three files: 

1. [*FallingBlock.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/FallingBlock/FallingBlock.m)<br>
   -> A general script to calculate the falling block problem (define model geometry, viscosity ratio, density ratio, and time). The parameters are not scaled for all scripts here. 

2. [*FallingBlockBenchmark.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/FallingBlock/FallingBlockBenchmark.m)<br>
   -> A script to calculate the instantaneous sinking velocity of a block with a certain geometry and density ratio for differnet viscosity ratios between the block and the background material.

3. [*FallingBlockBenchmark_timedependent.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/FallingBlock/FallingBlockBenchmark_timedependent.m)<br>
   -> Script to calculate the time-dependent solution of the previous script. The final time for each model is set to the time given by *Gerya* (2009). Here, the code does advect the composition with passive tracers to transport density and viscosity.

------------------------------------------------------------------

## **Constants**

Model height **H** [km]: 500 <br>
Model width **W** [km]: 500 <br>
Graviational acceleration **g** [m/s<sup>2</sup>]: 9.81 <br>
Medium viscosity $\eta_0$ [Pa s]: 10<sup>21</sup> <br>
Block viscosity $\eta_1$ [Pa s]: 10<sup>15</sup> -  10<sup>27</sup> <br>
Medium density $\rho_0$ [kg/m<sup>3</sup>]: 3200 <br>
Block density $\rho_1$ [kg/m<sup>3</sup>]: 3300 <br>

nx = 51<br>
nz = 51

Number of tracers per cell: nmx*nmz<br>
nmx = 5 <br>
nmz = 5

------------------------------------------------------------------

## Instantaneous Solution

![Field1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/5a53c246-d05e-45ef-babe-4fcd9dfab735)<br>
**Figure 1.** Denstiy [kg/m<sup>3</sup>] and viscosity [Pa s] fields for a viscosity ratio between the block and the background medium of 6 orders of magnitude. The arrows show the instantaneous solution for the velocity field. 

![drho_100_vz_eta_r_nx_51](https://github.com/LukasFuchs/FDCSGm/assets/25866942/ab71d51e-21da-4185-9d03-54e8e57c0b70)<br>
**Figure 2.** Maximum sinking velocity [m/s] over the logarithm of the viscosity ratio between the block and the background medium. The absolut viscosity deviates a little bit from the solution of *Gerya* (2009), but the trend of a decreasing sinking velocity with increasing viscosity ratio is the same. 

--------------------------

## Time Dependent

&emsp;The properties from the tracers are interpolated onto the regular finite differnece grid using a bilinear interpolation scheme. Thereby, the interpolation of the viscosity is limited to a cell of dx*dz around a grid node, where the density is interpolated within a cell of 2 * dx 2 * dz. The material defroms and sinks similar to the ones presented in *Gerya* (2009).  

![TracerCompTPRT05](https://github.com/LukasFuchs/FDCSGm/assets/25866942/b63c1fc0-47cc-4e52-b9fc-f62b981d4e92)<br>
**Figure 3.** Density fields [kg/m<sup>3</sup>] after a certain time for a sinking block model with different viscosity ratios (from top left to bottom right: -1,1,2,3,4,5, respectively). Shown is the density field overlaid by the tracers distribution. 
