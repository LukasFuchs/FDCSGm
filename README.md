
# FDCSGm

   This is a two-dimensional **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) in **m**ATLAB (**FDCSGm**) I wrote, mainly for teaching purposes. The students learn how to discretize the conservation equations of energy, momentum, and mass using different discretization forms (explicit, implicit, CNV, ADI, etc.), how to couple the momentum equation with the energy equation (here, a simple operator splitting method), and finally build a simple two-dimensional thermal convection code for isoviscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is an effective way to teach numerical methods in combination with geodynamical problems as well as to check simple geodynamic settings. 
   
   The staggered grid method, how to discretize and solve the stokes equation, and how to index those variables are based on the methods described by *Gerya (2010)*, as well as the coupling between the momentum and energy equation. The energy equation can be discretized in multiple different ways (explicit, implicit, Crank Nicolson, Alternating-Direction implicit) and also be used in the thermal convection code. The advection equation can be discretized in multiple ways (upwind, semi-lagrangian, and passive tracers). For the thermal convection code, temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far only the advection of the absolute temperature is implemented, though) if necessary (e.g., for the falling block or rigid body rotation benchmark, a growth rate of a Rayleigh-Taylor instability test, and the viscous inclusion benchmark). 
   
   The *MasterFile.m* in the main directory (*FDCSGm/*) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple benchmarks included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

# Energy equation
   The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal within a closed system. In terms of a geodynamical problem, the energy can be described by the temperature, which is transported mainly through *conductive* and *convective* processes, such that a simple general energy equation is defined as follows (assuming only radioactive heat sources):

$(\frac{\partial E}{\partial t} + \overrightarrow{v} \cdot \nabla E) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H$,&emsp;&emsp;&emsp;(1)

where the energy can be described as $E=c_{p} \rho T$, *c<sub>p</sub>* is the specific heat capacity [J/kg/K], *ρ* is a reference density [kg/m<sup>3</sup>], *T* is the temperature [K], *t* is the time, $\overrightarrow{v}$ is the velocity vector, *q<sub>i</sub>* is the heat flux in direction *i*, *∂/∂xi* is a directional derivative in direction *i*, and *H* the heat production rate per mass [W/kg]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as follows: 

$\overrightarrow{q} = - k \nabla T$,&emsp;&emsp;&emsp;(2)

where *k* is the thermal conductivity [W/m/K]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation* in an Eulerian form can then be written as: 

$\rho c_p (\frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T) = -\frac{\partial q_i}{\partial x_i} + \rho H$,&emsp;&emsp;&emsp;(3) 
 
   This form of the temperature equation describes the variation of temperature due to a *conductive* (right hand side of the equation) and *convective* (left hand side of the equation) processes. For a matter of simplicity, one can consider those terms in a separate manner  to solve the energy equation. Thus, we first focus on the conductive part of equation (3) and the equation simplifies (in 2-D) to (neglecting the convective term):

$\rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} -\frac{\partial q_z}{\partial z} + \rho H$,&emsp;&emsp;&emsp;(4)

or including Fourier’s law (assuming variable thermal parameters):

$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial z} k \frac{\partial T}{\partial z} + \rho H$.&emsp;&emsp;&emsp;(5) 

Assuming that the thermal parameters are constant, equation (5) simplifies to: 

$\frac{\partial T}{\partial t} = \kappa (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2}) + \frac{Q}{\rho c_p}$,&emsp;&emsp;&emsp;(6)
  
where κ is the thermal diffusivity [m<sup>2</sup>/s] and Q is the heat production rate per volume [W/m<sup>3</sup>]. Equation (6) is a *parabolic partial differential equation* which can be solve numerically in different manners, assuming initial and boundary conditions are defined. 

   Here we would like to discuss a simple, but effective finite difference method to discretize and solve the equation, that is the forward in time and centered in space (FTCS) method in an explicit manner. The advantage of an explicit description is that it is simple to derive and rather fast computationally, however, it is only numerically stable as long as the *heat diffusive stability criterion* is fulfilled. The stability criterion can be determined by a simple Von Neumann stability analysis, which analysis the growth of an eigenmode perturbation for a certain finite difference approach. In case of an explicit 2-D finite difference approach, the heat diffusive stability criterion is defined as $\Delta t < \frac{\Delta x^2}{3 \kappa}$, that is the time step is limited by the model’s resolution. Within the code different discretization methods can be chosen to solve the diffusive part of the temperature conservation equation (i.e., implicit, CNV, ADI). In the following, we would like to focus on the explicit method as an example of how to discretize the equation. For the other discretization methods check the numerical methods section below. 

   Using a general finite difference discretization method to approximate the partial derivatives from equation (6), the temperature equation can be written as:

$\frac{T_{i,j}^{n+1} - T_{i,j}^{n} }{\Delta t} = \kappa (\frac{T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}}{(\Delta x)^2} + \frac{T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}}{(\Delta x)^2}) + \frac{Q_{i,j}^n}{\rho c_p}$&emsp;&emsp;&emsp;(7)

where *i*,*j* are the vertical and horizontal indices of the numerical finite difference grid, *n* is the time step index, Δ*t* is the time step, and Δ*x* and Δ*z* are the widths of the grid in horizontal and vertical direction. This equation contains know and unknow parameters and we can rearrange them to solve the equation for the unknowns as:

$T_{i,j}^{n+1} = T_{i,j}^{n} + s_x(T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}) + s_z(T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}) + \frac{Q_{i,j}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(8)

where $s_x = \frac{\kappa \Delta t}{(\Delta x)^2}$ and $s_z = \frac{\kappa \Delta t}{(\Delta z)^2}$. Equation (8) can be solved *iteratively* for every inner grid point assuming an initial condition is defined (multiple initial conditions can be set in the code). 

   Different conditions can be set for the boundaries of our modell domain. Here, we focus on two fundamental conditions, the *Dirichlet* and *Neumann* boundary conditions. The Dirichlet boundary condition simply define a constant temperature along the boundary, such that the temperature, for example, along the *left* boundary can be defined as:

$T_{i,j=1} = T_{left}$, for all i. &emsp;&emsp;&emsp;(9)

   The same applies to all other boundaries (*right*, *top*, and *bottom*). The Neumann boundary condition defines that the variation of a certain parameter does not change across the boundary, that is the temperature across the boundary or thermal heat flux *q* through the boundary. A sophisticated method to describe the heat flux across the boundary using finite differences is assuming the existence of so-called ghost nodes outside of our modell domain, that are nodes which we can use to describe the flux across a boundary, but which do not actually exist. Therefore, one first needs to define the variation of the temperature across the boundary (e.g., the *left*, for *i* = 2,*nz*-1) as:

 $\frac{\partial T}{\partial x}|_{i,j} = c_l$,&emsp;&emsp;&emsp;(10)

or using finite differences: 

$\frac{T_{i,2} - T_{i,0}}{2 \Delta x} = c_l$,&emsp;&emsp;&emsp;(11)
 
where *c<sub>l</sub>* is a constant defining the flux across the boundary and *T<sub>i,0</sub>* are the ghost nodes outside the left boundary. Now, we can solve for an expression of the temperature at the ghost nodes which fulfills the condition of equation (10) as:

$T_{i,0} + T_{i,2} - 2 \Delta x c_l$.&emsp;&emsp;&emsp;(12)

Assuming equation (8) is also valid along the left boundary nodes and using the condition for the ghost nodes outside the numerical domain, that is we assume Neumann boundary conditions, one can rewrite equation (8) for the left boundary nodes as follows:

$T_{i,1}^{n+1} = T_{i,1}^{n} + s_x(2T_{i,2}^{n} - 2(T_{i,1}^{n} + \Delta x c_l)) + s_z(T_{i+1,1}^{n} - 2T_{i,1}^{n} + T_{i-1,1}^{n}) + \frac{Q_{i,1}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(13)
 
The same applies to the other boundaries. Caution needs to be taken at the corners of the model. 

   The routines for the explicit, implicit, and ADI discretization methods are available in a dimensional and non-dimensional way. For more details on the scaling of the parameters see *FDCSGm/ScaleParam*. Besides the routines to solve the diffusive part of the 2-D temperature conservation equation using different discretization methods there are additional  routines available within the code. There are two routines to calculate 1-D temperature profiles for a oceanic and continental lithosphere based on an 1-D implicit solver for constant and variable thermal parameters. A 1- and 2-D Poisson solver for a constant and variable *k* is available, to solve the steady-state diffusive equation (i.e. $\frac{\partial T}{\partial t} = 0$). For more details see *FDCSGm/DiffusionProblem/README.md*.

# Stokes equation 

# Advection equation 

# Scaling

# Numerical Methods

# Benchmarks 
## Blanckenbach 
## Channel Flow
## Falling Block
## Gaussian Diffusion 
## Rayleigh Taylor Instability 
## Spherical Viscous Inclusion 

# References
*Gerya, T. (2010). Introduction to numerical geodynamic modelling.*
