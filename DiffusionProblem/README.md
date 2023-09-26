# General Information

&emsp;This directory contains several rountines to solve the diffusive part of the *temperature conservation equation* (1- and 2-D, stationary and time-dependent) using different numerical discretization methods. The routines are avaible for a dimensional or non-dimensional (files ending with *Sc.m) form of the equation (so far the 1-D routines are only availabe for a dimensional version!). 

-------------

The general *temperature equation* describes the variation of temperature due to a *conductive* and *convective* process. However, for certain situations, one can consider those terms in a separate manner and the equation simplifies (in 2-D) to (neglecting the convective term):

$\rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} -\frac{\partial q_z}{\partial z} + \rho H$,&emsp;&emsp;&emsp;(1)

or including Fourier’s law (assuming variable thermal parameters):

$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial z} k \frac{\partial T}{\partial z} + \rho H$.&emsp;&emsp;&emsp;(2) 

Assuming that the thermal parameters are constant, equation (2) simplifies to: 

$\frac{\partial T}{\partial t} = \kappa (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2}) + \frac{Q}{\rho c_p}$,&emsp;&emsp;&emsp;(3)
  
where *κ* is the thermal diffusivity [m<sup>2</sup>/s] and $Q=\rho H$ is the heat production rate per volume [W/m<sup>3</sup>]. 

&emps;Equation (3) is a *parabolic partial differential equation* which can be solve numerically in different manners, assuming initial and boundary conditions are defined. A detailed description on how to solve equation (3) using an *explicit* finite difference scheme and how to implement the most common boundary conditions (*Dirichlet* and *Neumann*) is given in the [introduction](https://github.com/LukasFuchs/FDCSGm/tree/main#energy-equation) of this code. In the following I would like to show some additional ways to discretize the diffusive part of the energy equation and discuss their advantages and disadvantages a little bit. All discretization methods can be used in the [thermal convection code](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems) and the [Blankenbach Benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach). A more detailed analysis on the accuracy of each discretization scheme and the effect of the grid resolution is given in the [Gaussian Diffusion Benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/GaussDiffusion). So far, however, variable thermal parameters are only included in the 1-D solutions and the 2-D steady state solution (i.e., $\frac{\partial T}{\partial t}=0$).  

---------------------

## Example: Geotherms

### 1-D Geotherms

&emsp;The 1-D temperature profile is calculated by solving the diffusive part of the 1-D temperature conservation equation (so far only with a radiogenic heat source) for variable thermal parameters with a proper conserving finite difference scheme. That is, the heat flow is calculated on the centered and the remaining parameters on the regular grid points, respectively. The discretization scheme for variable thermal parameters is picked to solve for a temperature profile of a continental lithosphere with upper, lower crust, and mantle. 
The 1-D temperature equation is given by: 

$\rho c_{p} \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}(k \frac{\partial T}{\partial z}) + \rho H$, &emsp; &emsp; &emsp; (1)

where $\rho, c_{p}, T, t, k, H, z$ are the density [kg/m<sup>3</sup>], the specific heat capacity [J/kg/K], the temperature [K], the time [s], the thermal conductivity [W/m/K], the heat generation rate per mass [W/kg], and the depth [m] respectively. 

&emsp;Here, a proper conservative finite difference scheme means that the heat flux is calculated on the centered grid points (A, B, etc.). The 1-D vertical heat flux is given by the Fourier’s law:

$q_{z} = -k \frac{\partial T}{\partial z}$. &emsp; &emsp; &emsp; (2)

***Solving the equation***

&emsp;Following the discretization shown in Figure 1 one needs to solve the following equation (in an [implicit finite difference formulation](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dimplicit_vary.m)):

$\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = -\frac{q_{z,j+1/2}^{n+1} - q_{z,j-1/2}^{n+1} }{\Delta z} + \rho_j H_j$, &emsp;&emsp;&emsp; (3)

$\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = \frac{ k_{j+1/2} \frac{T_{j+1}^{n+1} - T_{j}^{n+1}}{\Delta z} - k_{j-1/2} \frac{T_{j}^{n+1} - T_{j-1}^{n+1}}{\Delta z} }{\Delta z} + \rho_j H_j$. &emsp;&emsp;&emsp; (4)

Sorting the variables (known variables on the right-hand side, unknown on the left-hand side): 

$\frac{k_{j+1/2}}{\Delta z^2} T_{j+1}^{n+1} - \frac{k_{j+1/2}}{\Delta z^2} T_{j}^{n+1} - \frac{k_{j-1/2}}{\Delta z^2} T_{j}^{n+1} + \frac{k_{j-1/2}}{\Delta z^2} T_{j-1}^{n+1} = \frac{\rho_j c_{p,j}}{\Delta t} T_{j}^{n+1} - \frac{\rho_j c_{p,j}}{\Delta t} T_{j}^{n} - \rho_j H_j$, &emsp;&emsp;&emsp;(5)

$\frac{k_{j+1/2}}{\Delta z^2} \frac{\Delta t}{\rho_j c_{p,j}} T_{j+1}^{n+1} - \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{j+1/2} + k_{j-1/2}}{\Delta z^2}) T_{j}^{n+1} - T_{j}^{n+1} + \frac{k_{j-1/2}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}} T_{j-1}^{n+1} = -T_{j}^{n} - \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(6)

$-\frac{k_{j+1/2}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}T_{j+1}^{n+1} + (1 + \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{j+1/2} + k_{j-1/2}}{\Delta z^2}))T_{j}^{n+1} - \frac{k_{j-1/2}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}} T_{j-1}^{n+1} = T_j^n + \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(7)

$aT_{j-1}^{n+1} + bT_{j}^{n+1} + cT_{j+1}^{n+1} = T_j^n + \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(8)

with

$a = -\frac{k_{j-1/2}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}$

$b = (1 + \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{j+1/2} + k_{j-1/2}}{\Delta z^2}))$, &emsp;&emsp;&emsp; (9) 

$c = - \frac{k_{j+1/2}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}$

and

$k_{j+1/2} = \frac{k_j + k_{j+1}}{2}$

$k_{j-1/2} = \frac{k_{j-1} + k_{j}}{2}$. &emsp;&emsp;&emsp; (10)

An [*explicit*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dexplicit_vary.m) solver for a 1-D thermal profile with variable thermal parameters and a radiogenic heat source is also available.

***Thermal boundary conditions***

The thermal boundary conditions are defined as: 

1. Constant temperature (*Dirichlet*)
The temperature at the top or bottom can just be set as constant to *T<sub>top</sub>* or *T<sub>bot</sub>*, respectively.
      
2. Constant temperature gradient (*Neumann*)
The gradient of temperature (and thus the vertical heat flux) can be defined using so called ghost nodes at the top and the bottom of the profile. Therefore we define the condition at the top and bottom as:

   $\frac{\partial T}{\partial z} \vert_{j=1} = c_{top} = \frac{T_2-T_0}{2\Delta z}$

   $\frac{\partial T}{\partial z} \vert_{j=nz} = c_{bottom} = \frac{T_{nz+1}-T_{nz-1}}{2\Delta z}$, &emsp;&emsp;&emsp; (11)

   Where T<sub>0</sub> and T<sub>nz+1</sub> are the ghost nodes for temperature at the top and bottom, respectively. The constants *c<sub>top</sub>* and *c<sub>bottom</sub>* are defined as:

   $c_{top,bottom} = -\frac{q_{top,bottom}}{2\Delta z}$
   
   Using these conditions, we can define formulations for the temperature at the ghost nodes as:

   $T_0 = T_2 - 2\Delta z c_{top}$

   $T_{nz+1} = T_{nz-1} + 2\Delta z c_{bottom}$

   Now one can solve equation (3) for the top and the bottom using the formulations of the temperature at the ghost nodes with equation (13), which results in:

   $\rho_{j=1}c_{p,j=1}\frac{\partial T}{\partial t}\vert_{j=1} = \frac{\partial}{\partial z}(k\frac{\partial T}{\partial z})$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{ k_A\frac{dT}{dz}\vert_A - k_B\frac{dT}{dz}\vert_B }{\Delta z}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{ k_A\frac{T_2^{n+1} - T_1^{n+1}}{\Delta z} - k_B\frac{T_1^{n+1} - T_0^{n+1}}{\Delta z} }{\Delta z}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{k_A}{\Delta z^2}T_2^{n+1} - \frac{k_A}{\Delta z^2}T_1^{n+1} - \frac{k_B}{\Delta z^2}T_1^{n+1} + \frac{k_B}{\Delta z^2}T_0^{n+1}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{k_A}{\Delta z^2}T_2^{n+1} - \frac{k_A+k_B}{\Delta z^2}T_1^{n+1} + \frac{k_B}{\Delta z^2}(T_2^{n+1} - 2\Delta z c_{top})$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = - \frac{k_A + k_B}{\Delta z^2}T_1^{n+1} + \frac{k_A+k_B}{\Delta z^2}T_2^{n+1} - \frac{k_B}{\Delta z}2c_{top}$

   $(1 + \frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2} (k_A+k_B) ) T_1^{n+1} - \frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2} (k_A+k_B) T_2^{n+1} = T_1^n - \frac{2\Delta t c_{top}}{ \rho_1 c_{p,1} \Delta z} k_B$

   $aT_1^{n+1}+bT_2^{n+1}=T_1^n+Q_{top}$, &emsp;&emsp;&emsp; (14)

   with

   $a = 1+\frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2}(k_A+k_B)$

   $b = -\frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2}(k_A+k_B)$

   $Q_{top} = -\frac{2\Delta t c_{top}}{\rho_1 c_{p,1} \delta z}k_A$. &emsp;&emsp;&emsp; (15)

   Similar for the bottom boundary with:

   $aT_{nz}^{n+1}+bT_{nz-1}^{n+1}=T_{nz}^{n} + Q_{bottom}$

   $a = 1+\frac{\Delta t}{\rho_{nz} c_{p,nz} \Delta z^2}(k_A+k_B)$

   $b = -\frac{\Delta t}{\rho_{nz} c_{p,nz} \Delta z^2}(k_A+k_B)$

   $Q_{top} = -\frac{2\Delta t c_{bottom}}{\rho_{nz} c_{p,nz} \delta z}k_A$. &emsp;&emsp;&emsp; (16)

#### Oceanic Geotherms
![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/247cb1ca-b0d1-40b0-bfdb-aaac00a66222)
***Figure 2. Oceanic Lithosphere.** LEFT: Temperature profile [K]  for an oceanic lithosphere of 60 Ma of age and constant thermal boundary conditions at the top and bottom. The blue line shows the initial temperature profile. The yellow dashed line shows the solution for a half-space cooling model. RIGHT: Heat flux [mW/m<sup>2</sup>] with depth. The parameters of this model are defined as the default values in the routine OceanicGeotherm.m.*

![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/6b13a316-ec6f-4ebc-9720-3290881a9b49)
***Figure 3. Oceanic Lithosphere II**. Same as Figure 1 but with constant heat flux boundary conditions qbottom =10 mW/m<sup>2</sup> and qtop = 90 mW/m<sup>2</sup>.*

#### Continental Geotherms
![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/3e7bdc48-24dc-48d8-a058-49b01cc819da)
***Figure 4. Continental Lithosphere.** LEFT: Temperature profile for a continental lithosphere of 1000 Ma of age with constant upper and lower thermal boundary conditions. The blue line shows the initial condition, the red line shows the solution of equation (1), the yellow dashed line shows the solution of the time-independent heat equation (1-D poisson equation), and the magenta dashed line shows the solution of a 2D, staggered finite difference code. MIDDLE: Heat flux with depth. RIGHT: Thermal parameter for the lithosphere setup: thermal conductivity [k], specific heat [c<sub>p</sub>], density [ρ], and volumetric heat generation rate [Q].*

![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/98e9da70-f343-4ed9-be71-a44299116c72)
***Figure 5. Continental Lithosphere II**. Same as Figure 3 but with constant upper and lower heat flux boundary conditions, q<sub>top</sub> = 40 mW/m<sup>2</sup> and q<sub>bottom</sub> = 10 mW/m<sup>2</sup>.*

## Steady State Solution

## Discretization Methods
...


Needs, in detail:
- discretization of     
     - implicit
     - ADI
     - CNV
     - Poisson solution, constant and variable,
     - explicit (variable thermal parameters, at some point!)

Needs briefly: 
- boundary conditions

