# General Information

This directory contains two files: 

1. *RigidBodyRotation.m*<br>
    -> Script to test the efficiency of a certain advection scheme for a specific background velocity field.
   
2. *RigidBodyRotationComp.m*<br>
    -> Similar as the script above, but it executes the different advection schemes conecutively to compare the effect of different numerical parameters on the advection scheme, e.g., the grid resolution. The final result is plotted in one figure.

-----------------------------------------------------------------

A rigid body rotation test is an effective benchmark to check the efficiency of a specific advection scheme. Ideally, an initial perturbation (depending of which property is advected; here it is temperature) should end at the same position after one full rotation and should not lose its initial shape. 

One can choose one of three different initial perturbations: 

1. A rectangular, constant perturbation<br>
    -> Choose 'block' as B.Tini; a rectangular block of certain width and height is positioned on the left side of the model domain, with a perturbation amplitude of B.TAmpl and a background (temperature) of B.T0. <br>

2. A gaussian perturbation<br>
    -> Choose 'gaussianRBR' as B.Tini; a gaussian temperature distribution is positioned on the left side of the model domain (for more details see */FDCSGm/Benchmark/GaussDiffusion/README.md*). <br>

3. A circular, constant perturbation<br>
    -> Choose 'circle' as B.Tini; a circular perturbation of certain radius (*B.Tsigma* * *M.L*) and amplitude (*B.TAmpl*) relative to the background (*B.T0*) is positioned on the left side of the model domain.<br>

As simple as this problem sounds, it is rather difficult to preserve the initial shape, mainly due to numerical diffusion effects or due to inaccuracies of an interpolation. Here, I use four different advection schemes to advect the temperature: 

1. *The upwind scheme*<br>
    -> An stable and effective way, however, with a certain amount of numerical diffusion if the *courant criteria* is not fulfilled. The courant criteria implies that the time step is smaller than the minimum grid spacing divided by the maximum velocity, that is, a property should not be advected over a distance larger than the grid spacing, or:<br>
$\Delta t \le \frac{\Delta x}{max(|v|)}$<br>
The upwind scheme is similar to a forward in time and centered in space discretization, however, one needs to consider the advection velocity as well to ensure that the discretization in space is always upstream.<br>
   
2. *The Staggered Leap Frog (SLF) scheme*<br>
    -> This method considers a centered in time and centered in space discretization of the partial differentials, thus it has a higher order of error and is suppose to not have any numerical diffusion. As promissing as this scheme sounds it is not properly working here yet (I believe)!<br>
   
3. *The semi-lagragian scheme*<br>
    -> This method assumes that an *imaginary tracer* is located at a certain postion and lands directly at a finite difference grid node after advection within one time step. Thus, one needs to calculate the *origin point* from each grid node back in time with the given velocity field (using a central point iteration method) and then interpolates the property from the regular grid points to the determined *origin points*. The method does not have any numerical diffusion but shows inaccuracies due to the interpolation method.<br>
   
5. *Passive tracers*<br>
    -> Here, one assumes that the model domain is completely filled with so-called tracers or markers. These tracers are then advected by a certain method (e.g., Euler or Runge Kutta) and transport any property stored on them. However, care needs to be taken when interpolating those properties from the regular grid onto the tracers and back. This is even more complex if the property advected does have an effect on parameters controlling the governing equations (e.g., the viscosity in continuum euqation). Here, I advect the tracers using Runge-Kutta fourth order; the tracers do transport the absolute temperature, which is only interpolated to the regular grid points every time step and not back to the tracers again (since it is not suppose to change here).<br>

To avoid the effect of boundary conditions, I set the velocity outside of the circular rigid body rotation field to zero. 

----------------------------------------------------

![Field1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f632c6e3-2051-45d8-ae48-c62a40ac2242)
**Figure 1.** Initial setup for a circular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is one inside. 

![2D_Advection_RigidBody](https://github.com/LukasFuchs/FDCSGm/assets/25866942/66ea8ad6-a277-4cd5-a91e-627f1b3f41fd)
**Figure 2.** Example of an evolution of a full rotation. The advection scheme is the semi-lagragian scheme for a grid resolution of 301x301, shown at every 50th iteration step. The shape of the initial perturbation is preserved rather well and only minor *diffusion* due to the interpolation of the temperature is obervable. 

![Comparison_301_301_circular](https://github.com/LukasFuchs/FDCSGm/assets/25866942/e3bea260-2b1c-4f4e-9da7-79cc429069f4)
**Figure 3.** Comparison of the final step for each advection method. For more details see the titles of each subplot. 

![Comparison_101_101_circular](https://github.com/LukasFuchs/FDCSGm/assets/25866942/ca206459-0f94-4500-8982-61b030ae71b8)
**Figure 4.** Same as *Figure 3*, but with a smaller resolution of 101x101. For more details see the titles of each subplot. 

![Initial_101_101_gaussian](https://github.com/LukasFuchs/FDCSGm/assets/25866942/c6b09e02-5ce4-4cd7-b9a6-fbd78cebf3a4)
**Figure 5.** Initial setup for a gaussian (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus the maximum is one.

![Comparison_101_101_gaussian](https://github.com/LukasFuchs/FDCSGm/assets/25866942/c516f57b-4860-4a68-9e25-41e50a271d08)
**Figure 6.** Same as *Figure 4*, but with a gaussian (temperature) perturbation. For more details see the titles of each subplot. 

![Initial_101_101_block](https://github.com/LukasFuchs/FDCSGm/assets/25866942/3fe2ccf2-a341-4446-bca1-c274bf2b1c75)
**Figure 7.** Initial setup for a rectangular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is one inside.

![Comparison_101_101_block](https://github.com/LukasFuchs/FDCSGm/assets/25866942/2c1c431c-7b63-4f82-8a6b-34894afc8261)
**Figure 8.** Same as *Figure 4*, but with a rectangular, constant (temperature) perturbation. For more details see the titles of each subplot.






