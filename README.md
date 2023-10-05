# FDCSGm

&emsp;This is a two-dimensional **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) in **m**ATLAB (**FDCSGm**) I wrote, mainly for teaching purposes. The students learn how to discretize the conservation equations of energy, momentum, and mass using different discretization forms (*explicit*, *implicit*, *CNA*, *ADI*, etc.), how to couple the momentum equation with the energy equation (here, a simple *operator splitting* method), and finally build a simple two-dimensional thermal convection code for isoviscous and temperature-dependent viscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is an effective way to teach numerical methods in combination with geodynamical problems as well as to check simple geodynamic settings. 
   
&emsp;The staggered grid method, how to discretize and solve the stokes equation, and how to index these variables are based on the methods described by *Gerya (2009)*. The energy equation can be discretized in multiple different ways (e.g., *explicit*, *implicit*, *Crank Nicolson*, *Alternating-Direction implicit*) which can also be used in the thermal convection code. The advection equation can be discretized in multiple ways, too (e.g., *upwind*, *semi-lagrangian*, and *passive tracers*). For the thermal convection code (e.g., the [Blankenback benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach) or for [mixed heated systems](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems), temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far, only the advection of the **absolute** temperature is implemented, though) if necessary (e.g., for the [falling block](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/FallingBlock) or [rigid body rotation benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RigidBodyRotation), a growth rate of a [Rayleigh-Taylor instability](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RTI) test, and the [viscous inclusion benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/ViscousInclusion)). 
   
&emsp;The [*MasterFile.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/MasterFile.m) in the main directory (*FDCSGm/*) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple [benchmarks](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark) included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

-----------------

&emsp;For the given scaling parameters, the non-dimensional governing equations are given as (assuming a constant viscosity; scaling for a variable viscosity is applicable in the same way):

$\frac{\partial v_x}{\partial x} + \frac{\partial v_z}{\partial z} = 0$,&emsp;&emsp;&emsp;(32)

$\frac{\partial T}{\partial t} + v_x \frac{\partial T}{\partial x} + v_z \frac{\partial T}{\partial z} = (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2} + H)$,&emsp;&emsp;&emsp;(33)

$-\frac{\partial P}{\partial x} + \eta \frac{\partial^2 v_x}{\partial x^2} + \eta \frac{\partial^2 v_x}{\partial z^2} = 0$,&emsp;&emsp;&emsp;(34)

$-\frac{\partial P}{\partial z} + \eta \frac{\partial^2 v_z}{\partial z^2} + \eta \frac{\partial^2 v_z}{\partial x^2} - RaT = 0$,&emsp;&emsp;&emsp;(35)

where *Ra* is the so-called thermal *Rayleigh number* and *P* the *dynamic pressure*. In case of a basally heated thermal convection, the convective vigor is defined by the Rayleigh number, which describes a relationship between heat transported by buoyancy and conduction, and the effect of the layers thickness and bulk viscosity.

# References

*Gerya, T. (2009). Introduction to numerical geodynamic modelling.*

*Becker, T. W. and Kaus, B. J. P (2020 update): Numerical Modeling of Earth Systems: An introduction to computational methods with focus on solid Earth applications of continuum mechanics, The University of Texas at Austin, v. 1.2.2, 2020.*

