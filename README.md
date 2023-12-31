# FDCSGm

&emsp;This is a two-dimensional **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) I wrote in **m**ATLAB (**FDCSGm**), mainly for [teaching purposes](https://lukasfuchs.wordpress.com/teaching-2/). The students learn how to discretize the conservation equations of energy, momentum, and mass using different discretization forms. Here, for the sake of a better teaching experience (at least in my opinion), I do consider the energy equation in a seperate manner by first discussing the *conductive* part, second the *convective* part, and finally how to couple them (here, a simple *operator splitting* method, e.g., *Becker and Kaus* (2020)). The final goal of the course for the students is to build a simple two-dimensional thermal convection code for isoviscous and temperature-dependent viscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is an effective way to teach numerical methods in combination with geodynamical problems as well as to check simple geodynamic settings. 
   
&emsp;The staggered grid method, how to discretize and solve the [stokes equation](https://github.com/LukasFuchs/FDCSGm/tree/main/StokesProblem), and how to index these variables are based on the methods described by *Gerya (2009)*. The [conductive part of the energy equation](https://github.com/LukasFuchs/FDCSGm/tree/main/DiffusionProblem) can be discretized in multiple different ways (e.g., *explicit*, *implicit*, *Crank Nicolson*, *Alternating-Direction implicit*) which can also be used in the thermal convection code. The [advection equation](https://github.com/LukasFuchs/FDCSGm/tree/main/AdvectionProblem) can be discretized in multiple ways, too (e.g., *upwind*, *semi-lagrangian*, and *passive tracers*). 

&emsp;For the thermal convection code (e.g., the [Blankenback benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach) or for [mixed heated systems](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems)), temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far, only the advection of the **absolute** temperature is implemented, though) if necessary (e.g., for the [falling block](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/FallingBlock) or [rigid body rotation benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RigidBodyRotation), a growth rate of a [Rayleigh-Taylor instability](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RTI) test, and the [viscous inclusion benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/ViscousInclusion)). 
   
&emsp;The [*MasterFile.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/MasterFile.m) in the main directory (**FDCSGm/**) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple [benchmarks](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark) included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

For the colormaps I use the scientific colormaps by *Crameri (2021)*.

-----------------

# References

*Crameri, Fabio. (2021). Scientific colour maps (7.0.1). Zenodo. https://doi.org/10.5281/zenodo.5501399*

*Gerya, T. (2009). Introduction to numerical geodynamic modelling.*

*Becker, T. W. and Kaus, B. J. P (2020 update): Numerical Modeling of Earth Systems: An introduction to computational methods with focus on solid Earth applications of continuum mechanics, The University of Texas at Austin, v. 1.2.2, 2020.*

*...additional references will be added from time to time...*
