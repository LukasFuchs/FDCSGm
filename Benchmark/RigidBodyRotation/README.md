# General Information

This directory contains two files: 

1. *RigidBodyRotation.m*<br>
    -> Script to test the efficiency of a certain advection scheme for a specific background velocity field.
   
2. *RigidBodyRotationComp.m*<br>
    -> Similar as the script above, but it executes the different advection schemes conecutively to compare the effect of different numerical parameters on the advection scheme, e.g., the grid resolution. The final result is plotted in one figure.

-----------------------------------------------------------------

A rigid body rotation test is an effective benchmark to check the efficiency of a specific advection scheme. Ideally, an initial perturbation (depending of which property is advected; here it is temperature) should end at the same position after one full rotation and should not lose its initial shape. 

One can choose three different initial perturbation to be advected: 

1. A rectangular, constant perturbation<br>
    -> Choose 'block' as B.Tini; a rectangular block of certain widht and height is positioned on the left side of the model domain, with a perturbation amplitude of B.TAmpl and a background (temperature) of B.T0. <br>

2. A gaussian perturbation<br>
    -> Choose 'gaussian' as B.Tini; a gaussian temperature distribution is position on the left side of the model domain (for more details see */FDCSGm/Benchmark/GaussDiffusion/*). <br>

3. A circular, constant perturbation<br>
    -> Choose 'circle' as B.Tini; a circular perturbation of certain radius (*B.Tsigma* * *M.L*) and amplitude (*B.TAmpl*) relative to the background (*B.T0*) is positioned on the left side of the model domain.(<br>

As simple as this problem sounds, it is rather difficult to preserve the initial shape, mainly due to numerical diffusion effects or due to inaccuracies due to an interpolation. Here, I use four different advection schemes to advect the temperature: 

1. The upwind scheme<br>
    -> <br>
   
2. The Staggered Leaped Frog scheme<br>
    -><br>
   
3. The semi-lagragian scheme<br>
    -><br>
   
4. Passive tracers<br>
    -><br>

To avoid the effect of boundary conditions, I set the velocity outside of the circular rigid body rotation field to zero. 

----------------------------------------------------

![Field1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f632c6e3-2051-45d8-ae48-c62a40ac2242)
**Figure 1.** Initial setup for a circular, constant (temperature) perturbation. The temperature field is scaled by the maximum of the temperature, thus it is one inside. 

![2D_Advection_RigidBody](https://github.com/LukasFuchs/FDCSGm/assets/25866942/66ea8ad6-a277-4cd5-a91e-627f1b3f41fd)
**Figure 2.** Example of an evolution of a full rotation. The advection scheme is a semi-lagragian scheme for a grid resolution of 301x301, shown at every 50th iteration step. The shape of the initial perturbation is preserved rather well and only minor *diffusion* due to the interpolation of the temperature is obervable. 
