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

![Field1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f632c6e3-2051-45d8-ae48-c62a40ac2242)<br>
**Figure 1.** Initial setup for a circular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is equal to one inside. 

![2D_Advection_RigidBody](https://github.com/LukasFuchs/FDCSGm/assets/25866942/66ea8ad6-a277-4cd5-a91e-627f1b3f41fd)<br>
**Figure 2.** Example of an evolution of a full rotation. The advection scheme is the semi-lagragian scheme for a grid resolution of 301x301, shown at every 50th iteration step. The shape of the initial perturbation is preserved rather well and only minor *diffusion* due to the interpolation of the temperature is observable. 

![Comparison_301_301_circular](https://github.com/LukasFuchs/FDCSGm/assets/25866942/e3bea260-2b1c-4f4e-9da7-79cc429069f4)<br>
**Figure 3.** Comparison of the final step for each advection method. For more details see the titles of each subplot. 

![Comparison_101_101_circular](https://github.com/LukasFuchs/FDCSGm/assets/25866942/ca206459-0f94-4500-8982-61b030ae71b8)<br>
**Figure 4.** Same as *Figure 3*, but with a smaller resolution of 101x101. For more details see the titles of each subplot. 

![Initial_101_101_gaussian](https://github.com/LukasFuchs/FDCSGm/assets/25866942/c6b09e02-5ce4-4cd7-b9a6-fbd78cebf3a4)<br>
**Figure 5.** Initial setup for a gaussian (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus the maximum is equal to one.

![Comparison_101_101_gaussian](https://github.com/LukasFuchs/FDCSGm/assets/25866942/c516f57b-4860-4a68-9e25-41e50a271d08)<br>
**Figure 6.** Same as *Figure 4*, but with a gaussian (temperature) perturbation. For more details see the titles of each subplot. 

![Initial_101_101_block](https://github.com/LukasFuchs/FDCSGm/assets/25866942/3fe2ccf2-a341-4446-bca1-c274bf2b1c75) <br>
**Figure 7.** Initial setup for a rectangular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is equal to one inside.

![Comparison_101_101_block](https://github.com/LukasFuchs/FDCSGm/assets/25866942/2c1c431c-7b63-4f82-8a6b-34894afc8261)<br>
**Figure 8.** Same as *Figure 4*, but with a rectangular, constant (temperature) perturbation. For more details see the titles of each subplot.
