# General Information

This directory contains two files: 

1. [*RigidBodyRotation.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/RigidBodyRotation/RigidBodyRotation.m)<br>
    -> Script to test the efficiency of a certain advection scheme for a specific background velocity field.
   
2. [*RigidBodyRotationComp.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/RigidBodyRotation/RigidBodyRotationComp.m)<br>
    -> Similar as the script above, but it executes the different advection schemes conecutively to compare the effect of different numerical parameters on the advection scheme, e.g., the grid resolution. The final result is plotted in one figure.

For more details on the advection see [/FDCSGm/AdvectionProblem/](https://github.com/LukasFuchs/FDCSGm/tree/main/AdvectionProblem)

-----------------------------------------------------------------

## Perturbation setting

&emsp;A rigid body rotation test is an effective benchmark to check the efficiency of a specific advection scheme. Ideally, an initial perturbation (depending of which property is advected; here it is temperature) should end at the same position after one full rotation and should not lose its initial shape. 

One can choose one of three different initial perturbations: 

1. **A rectangular, constant perturbation**<br>
    -> Choose **'block'** as *B.Tini*; a rectangular block of certain width and height is positioned on the left side of the model domain, with a perturbation amplitude of *B.TAmpl* and a background (temperature) of *B.T0*.<br>

2. **A gaussian perturbation**<br>
    -> Choose **'gaussianRBR'** as *B.Tini*; a gaussian temperature distribution is positioned on the left side of the model domain (for more details on how the perturbation is defined, see [/FDCSGm/Benchmark/GaussDiffusion/](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/GaussDiffusion)).<br>

4. **A circular, constant perturbation**<br>
    -> Choose **'circle'** as *B.Tini*; a circular perturbation of certain radius (*B.Tsigma* * *M.L*) and amplitude (*B.TAmpl*) relative to the background (*B.T0*) is positioned on the left side of the model domain.<br>

------------------

## Advection scheme

&emsp;As simple as this problem sounds, it is rather difficult to preserve the initial shape, mainly due to numerical diffusion or due to inaccuracies of an interpolation. Here, I use *four* different advection schemes to advect the total temperature, where the courant criteria is equal to one (for more details, see [/FDCSGm/AdvectionProblem/](https://github.com/LukasFuchs/FDCSGm/tree/main/AdvectionProblem)): 

1. *The upwind scheme*<br>
   
2. *The Staggered Leap Frog (SLF) scheme*<br>
   
3. *The semi-lagragian scheme*<br>
   
5. *Passive tracers*<br>
  
To avoid the effect of boundary conditions as best as possible, I set the velocity outside of the circular rigid body rotation field to zero. 

----------------------------------------------------

## Examples of the advection schemes

![Field1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7ce5423c-4aa8-45f6-a5b7-d1b7710f85b7)<br>
**Figure 1.** Initial setup for a circular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is equal to one inside. 

![2D_Advection_RigidBody](https://github.com/LukasFuchs/FDCSGm/assets/25866942/81e45839-b93f-4631-ab39-70be21370ecb)<br>
**Figure 2.** Example of an evolution of a full rotation. The advection scheme is the semi-lagragian scheme for a grid resolution of 301x301, shown at every 50th iteration step. The shape of the initial perturbation is preserved rather well and only minor *diffusion* due to the interpolation of the temperature is observable. 

![Comparison_301_301_circle](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a3671efb-a046-4d17-8395-5093c4b111f4)<br>
**Figure 3.** Comparison of the final step for each advection method. For more details see the titles of each subplot. The initial condition is plotted on top.

![Comparison_101_101_Circular](https://github.com/LukasFuchs/FDCSGm/assets/25866942/6d9d0c90-cb81-4d3f-9bc6-a08ab3478903)<br>
**Figure 4.** Same as *Figure 3*, but with a smaller resolution of 101x101. For more details see the titles of each subplot. The initial condition is plotted on top.

![Comparison_101_101_gaussian_ini](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7f0b01f6-a4b4-4a98-98f7-18e1aea797ee)<br>
**Figure 5.** Initial setup for a gaussian (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus the maximum is equal to one.

![Comparison_101_101_Gaussian](https://github.com/LukasFuchs/FDCSGm/assets/25866942/63341da1-1c2f-4d45-91b9-84bf959e9f88)<br>
**Figure 6.** Same as *Figure 4*, but with a gaussian (temperature) perturbation. For more details see the titles of each subplot. The initial condition is plotted on top.

![Comparison_101_101_block_ini](https://github.com/LukasFuchs/FDCSGm/assets/25866942/0483945e-69ad-4d1d-9c09-e8b5b6f3cd7e)<br>
**Figure 7.** Initial setup for a rectangular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is equal to one inside.

![Comparison_101_101_block](https://github.com/LukasFuchs/FDCSGm/assets/25866942/bc566701-9a00-4509-bf4e-9ef691af59b4)<br>
**Figure 8.** Same as *Figure 4*, but with a rectangular, constant (temperature) perturbation. For more details see the titles of each subplot. The initial condition is plotted on top.

