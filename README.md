# FDCSGm

   This is a **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) for **m**ATLAB (**FDCSGm**). 
   
   I use the code mainly for teaching purposes. The students learn how to discretize the equations for the conservation of energy, momentum, and mass using different discretization forms (explicit, implicit, ADI, etc.), how to connect the momentum equation with the energy equation (a simple operator splitting method), and finally create a simple two-dimensional thermal convection code for an isoviscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is a very effective way to teach numerical methods in combination with geodynamical issues as well as to check simple geodynamic problems. 
   
   The staggered grid method and how to discretize and solve the stokes equation and how their variables are indexed are based on the methods described by *Gerya (2010)*, as well as the coupling between the momentum and energy equation. The energy equation can be discretized in multiple different ways (explicit, implicit, Crank Nicolson, Alternating-Direction implicit) and also be used in the thermal convection code. The advection equation can be discretized in multiple ways (upwind, semi-lagrangian, and passive tracers). For the thermal convection code, temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far only the advection of the absolute temperature is implemented, though) if necessary (e.g., for the falling block or rigid body rotation benchmark, a growth rate of a Rayleigh-Taylor instability test, and the viscous inclusion benchmark). 
   
   The *MasterFile.m* in the main directory (*FDCSGm/*) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple benchmarks included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

# Energy equation

# Stokes equation 

# Advection equation 

# Benchmarks 
## Blanckenbach 
## Channel Flow
## Falling Block
## Gaussian Diffusion 
## Rayleigh Taylor Instability 
## Spherical Viscous Inclusion 

# References
*Gerya, T. (2010). Introduction to numerical geodynamic modelling.*
