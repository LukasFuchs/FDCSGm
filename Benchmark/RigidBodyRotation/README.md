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
    ->

2. A gaussian perturbation<br>
    ->

3. A circular, constant perturbation<br>
    -> 

As simple as this problem sounds, it is rather difficult to preserve the initial shape, mainly due to numerical diffusion effects or due to inaccuracies due to an interpolation. Here, I use four different advection schemes to advect the temperature: 

1. The upwind scheme<br>
    -> 
   
2. The Staggered Leaped Frog scheme<br>
    ->
   
3. The semi-lagragian scheme<br>
    ->
   
4. Passive tracers<br>
    ->

