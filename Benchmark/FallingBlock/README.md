This directory contains three scripts: 
1. FallingBlock.m<br>
   -> A general script to calculate the falling block problem (define model geometry, viscosity ratio, density ratio, and time). The parameters are not scaled for all scripts here. 

2. FallongBlockBenchmark.m<br>
   -> A script to calculate the instantaneous sinking velocity of a block with a certain geometry and density ratio for differnet viscosity ratios between the block and the background material.

3. FallingBlockBenchmark_timedependent.m
   -> Script to calculate the time-dependent solution of the previous script. The final time for each model is set to the time given by *Gerya* (2009). Here, the code does advect the composition with passive tracers to transport the density and the viscosity.

------------------------------------------------------------------
## **Constants**<br>
Model height **H** [km]: 500 
Model width **W** [km]: 500
Graviational acceleration **g** [m/s<sup>2</sup>]: 9.81
Medium viscosity $\eta_0$ [Pa s]: 10<sup>21</sup>
Medium density $\rho_0$ [kg/m<sup>3</sup>]: 3200
Block density $\rho_1$ [kg/m<sup>3</sup>]: 3300


nx = 51
nz = 51

Number of tracers per cell: nmx*nmz
nmx = 5; 
nmz = 5

------------------------------------------------------------------

## Instantaneous Solution <br>

![Field](https://github.com/LukasFuchs/FDCSGm/assets/25866942/95af2887-d065-4be7-ab7f-a0b09d7bf73c)<br>
**Figure 1.** Denstiy [kg/m<sup>3</sup>] and viscosity [Pa s] fields. The arrows show the instantaneous solution for the velocity field for a viscosity ratio between the block and the background medium of 6 orders of magnitude. 

![vz_eta_r_nx_51](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a5fa6fc1-a989-4e27-bd4f-3a67562ecad9)<br>
**Figure 2.** Maximum sinking velocity [m/s] over the logarithm of the viscosity ratio between the block and the background medium. The absolut viscosity deviates a little bit from the solution of *Gerya* (2009), but the trend of a decreasing sinking velocity with increasing viscosity ratio is the same. 

## Time Dependent <br>
![TracerComparison_nx_51](https://github.com/LukasFuchs/FDCSGm/assets/25866942/320bdef2-dcab-478c-b46b-cea2e712c5bf)<br>
**Figure 3.** Density fields [kg/m<sup>3</sup>] for a sinking block model with different viscosity ratios. 
