# Routines to Solve the Diffusion Problem

   This directory contains all rountines to solve the diffusive part of the *temperature conservation equation* (1- and 2-D) using different numerical discretization methods. For more details on the discretization methods see the Numerical Methods section in *FDCSGm/README.md*. The routines are avaible for a dimensional or non-dimensional (files ending with *Sc.m) form of the equation (so far the 1-D routines are only availabe for a dimensional version!). 

## 1-D Geotherms
   The 1-D temperature profile is calculated by solving diffusive parte of the 1-D temperature conservation equation (so far only with a radiogenic heat source) for variable thermal parameters with a proper conserving finite difference scheme. That is, the heat flow is calculated on the centered and the remaining parameters on the regular grid points, respectively. The discretization scheme for variable thermal parameters is picked to solve for a temperature profile of a continental lithosphere with upper, lower crust, and mantle. 
The 1-D heat equation is given by: 

$\rho c_{p} \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}(k \frac{\partial T}{\partial z}) + \rho H$, &emsp; &emsp; &emsp; (1)

where $\rho, c_{p}, T, t, k, H, z$ are the density [kg/m<sup>3</sup>], the specific heat capacity [J/kg/K], the temperature [K], the time [s], the thermal conductivity [W/m/K], the heat generation rate per mass [W/kg], and the depth [m] respectively. 

   Here, a proper conservative finite difference scheme means that the heat flux is calculated on the centered grid points (A, B, etc.). The 1-D vertical heat flux is given by the Fourier’s law:

$q_{z} = k \frac{\partial T}{\partial z}$. &emsp; &emsp; &emsp; (2)

*Solving the equation* 
   Following the discretization shown in Figure 1 we need to solve the following equation (in an implicit finite difference formulation):

$\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = \frac{q_{z,j+1/2}^{n+1} - q_{z,j-1/2}^{n+1} }{\Delta z} + \rho_j H_j$, &emsp;&emsp;&emsp; (3)

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

*Thermal boundary conditions*

The thermal boundary conditions are defined as: 
   a) Constant temperature (*Dirichlet*)
      The temperature at the top or bottom can just be set as constant to *T<sub>top</sub>* or *T<sub>bot</sub>*, respectively.
   b) b)	Constant temperature gradient (Neumann)
The gradient of temperature (and thus the vertical heat flux) can be defined using so called ghost nodes at the top and the bottom of the profile. Therefore we define the condition at the top and bottom as:


### Oceanic Geotherms

### Continental Geotherms

## Steady State Solution
