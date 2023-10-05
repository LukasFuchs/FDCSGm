# FDCSGm

&emsp;This is a two-dimensional **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) in **m**ATLAB (**FDCSGm**) I wrote, mainly for teaching purposes. The students learn how to discretize the conservation equations of energy, momentum, and mass using different discretization forms (*explicit*, *implicit*, *CNA*, *ADI*, etc.), how to couple the momentum equation with the energy equation (here, a simple *operator splitting* method), and finally build a simple two-dimensional thermal convection code for isoviscous and temperature-dependent viscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is an effective way to teach numerical methods in combination with geodynamical problems as well as to check simple geodynamic settings. 
   
&emsp;The staggered grid method, how to discretize and solve the stokes equation, and how to index these variables are based on the methods described by *Gerya (2009)*. The energy equation can be discretized in multiple different ways (e.g., *explicit*, *implicit*, *Crank Nicolson*, *Alternating-Direction implicit*) which can also be used in the thermal convection code. The advection equation can be discretized in multiple ways, too (e.g., *upwind*, *semi-lagrangian*, and *passive tracers*). For the thermal convection code (e.g., the [Blankenback benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach) or for [mixed heated systems](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems), temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far, only the advection of the **absolute** temperature is implemented, though) if necessary (e.g., for the [falling block](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/FallingBlock) or [rigid body rotation benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RigidBodyRotation), a growth rate of a [Rayleigh-Taylor instability](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RTI) test, and the [viscous inclusion benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/ViscousInclusion)). 
   
&emsp;The [*MasterFile.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/MasterFile.m) in the main directory (*FDCSGm/*) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple [benchmarks](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark) included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

-----------------

# Stokes equation 
&emsp;On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one might consider, i.e., *inertia*, *surface*, and *volumetric* forces. A common equation to describes such motion is given by: 

$\rho \frac{D \overrightarrow{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \boldsymbol{\rho}$,&emsp;&emsp;&emsp;(14)

where *ρ* is the density [kg/m<sup>3</sup>], $\overrightarrow{v}$ is the velocity vector [m/s], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [Pa], $\boldsymbol{g}$ is the gravitational acceleration [m/s<sup>2</sup>], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial}{\partial t} + \overrightarrow{v} \cdot \nabla$. The *Cauchy stress tensor* is given by: 

$\boldsymbol{\sigma} = -\nabla P + \nabla \cdot \boldsymbol{\tau}$,&emsp;&emsp;&emsp;(15)

where *P* is the total pressure (*P = P<sub>dynamic</sub> + P<sub>hydrostatic</sub>*) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. In Eulerian form, equation (14) is given by (Navier-Stokes equation):

$\rho (\frac{\partial v_{i}}{\partial x} + v_{j}\frac{v_{i}}{\partial x_{j}}) = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(16)

where summation over repeated indices is implied. 

&emsp;To solve equation (16), one needs to define a rheology which, for a purely viscous medium, can be given by a constitutive relationship between stress and strain rate in the form of, e.g.:

$\tau_{ij} = 2 \eta \dot{\varepsilon}_{ij}$,&emsp;&emsp;&emsp;(17)

where $\eta$ is the dynamic viscosity [Pa s] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* [1/s] and given by: 

$\dot{\varepsilon}_{ij} = \frac{1}{2} (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i})$.&emsp;&emsp;&emsp;(18)

&emsp;Assuming that the inertia forces are negligible in comparison to the gravitational forces, one can further simplify equation (16) to:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(19)

or in the form of the unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\partial}{\partial x_j} \eta (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i}) + \rho g_{i}$,&emsp;&emsp;&emsp;(20)

Assuming constant viscosity equation (20) simplifies further to (Stokes equation): 

$0 = -\frac{\partial P}{\partial x_{i}} + \eta (\frac{\partial^2 v_i}{\partial x_j^2} + \frac{\partial^2 v_j}{\partial x_i^2}) + \rho g_{i}$.&emsp;&emsp;&emsp;(21)

&emsp;Equation (20) provides us two equations for our three unknowns. Thus, one needs to also consider the mass conservation equation, where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$\frac{\partial v_i}{\partial x_i} = 0$.&emsp;&emsp;&emsp;(22)

Equations (20) and (22) enable us to solve for the three unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*. 

&emsp;The conservation equations of *momentum* and *mass* are solved properly in two dimensions (*x* and *z*, *z* pointing downwards in negative direction) using a staggered finite difference grid, based on the description by *Gerya* (2009), where the viscosity and the density are defined on the regular grid points, the horizontal and vertical velocity in between the regular grid points, and the pressure within a finite difference cell. A staggered grid enables the conservation of the stress between adjacent grid points and one can solve equations (20)-(22) for the unknows. For more details on the discretization and the solving equations see *Gerya* (2009). 

# Scaling and equation of state
&emsp;To better compare different kinds of thermal convection, one can scale the governing equations. The equation can be scaled by the following reference parameters: 

$h_{sc} = h$,&emsp;&emsp;&emsp;(24)

$\eta_{sc} = \eta_0$,&emsp;&emsp;&emsp;(25)

$t_{sc} = \frac{h^2}{\kappa}$,&emsp;&emsp;&emsp;(26)

$v_{sc} = \frac{\kappa}{h}$,&emsp;&emsp;&emsp;(27)

$\tau_{sc} = \frac{\eta_0 \kappa}{h^2}$,&emsp;&emsp;&emsp;(28)

$T_{sc} = \Delta T$,&emsp;&emsp;&emsp;(29)

$H_{sc} = \frac{c_p \Delta T \kappa}{h^2}$,&emsp;&emsp;&emsp;(30)

where the subscript *sc* stands for the scaling parameters, and *h*, *η<sub>0</sub>*, *t*, *κ*, *v*, *τ*, *T*, *H*, *c<sub>p</sub>*, are the height, the reference viscosity, the time, the thermal diffusivity, the velocity, the stress, the temperature, the heat generation source, and the specific heat capacit, respectively. 

&emsp;The buoyance term on the right-hand side of equation (19), that is the density term which is temperature dependent (and pressure, but I do neglect this effect here so far), can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$\rho = \rho_0 (1-\alpha T)$,&emsp;&emsp;&emsp;(31)

where *ρ<sub>0</sub>* is the reference density and *α* the thermal expansion coefficient [1/K]. 

&emsp;For the given scaling parameters, the non-dimensional governing equations are given as (assuming a constant viscosity; scaling for a variable viscosity is applicable in the same way):

$\frac{\partial v_x}{\partial x} + \frac{\partial v_z}{\partial z} = 0$,&emsp;&emsp;&emsp;(32)

$\frac{\partial T}{\partial t} + v_x \frac{\partial T}{\partial x} + v_z \frac{\partial T}{\partial z} = (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2} + H)$,&emsp;&emsp;&emsp;(33)

$-\frac{\partial P}{\partial x} + \eta \frac{\partial^2 v_x}{\partial x^2} + \eta \frac{\partial^2 v_x}{\partial z^2} = 0$,&emsp;&emsp;&emsp;(34)

$-\frac{\partial P}{\partial z} + \eta \frac{\partial^2 v_z}{\partial z^2} + \eta \frac{\partial^2 v_z}{\partial x^2} - RaT = 0$,&emsp;&emsp;&emsp;(35)

where *Ra* is the so-called thermal *Rayleigh number* and *P* the *dynamic pressure*. In case of a basally heated thermal convection, the convective vigor is defined by the Rayleigh number, which describes a relationship between heat transported by buoyancy and conduction, and the effect of the layers thickness and bulk viscosity.

# References

*Gerya, T. (2009). Introduction to numerical geodynamic modelling.*

*Becker, T. W. and Kaus, B. J. P (2020 update): Numerical Modeling of Earth Systems: An introduction to computational methods with focus on solid Earth applications of continuum mechanics, The University of Texas at Austin, v. 1.2.2, 2020.*

