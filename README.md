# FDCSGm

   This is a two-dimensional **F**inite **D**ifference **C**ode with a **S**taggered **G**rid for the stokes variables (velocity and pressure) in **m**ATLAB (**FDCSGm**) I wrote, mainly for teaching purposes. The students learn how to discretize the conservation equations of energy, momentum, and mass using different discretization forms (explicit, implicit, CNV, ADI, etc.), how to couple the momentum equation with the energy equation (here, a simple operator splitting method), and finally build a simple two-dimensional thermal convection code for isoviscous convection with application (that is, the Blanckenbach benchmark). The style of the code is very flexible since I do optimize and vary routines from time to time, but it is an effective way to teach numerical methods in combination with geodynamical problems as well as to check simple geodynamic settings. 
   
   The staggered grid method, how to discretize and solve the stokes equation, and how to index those variables are based on the methods described by *Gerya (2009)*. The energy equation can be discretized in multiple different ways (explicit, implicit, Crank Nicolson, Alternating-Direction implicit) which can also be used in the thermal convection code. The advection equation can be discretized in multiple ways, too (upwind, semi-lagrangian, and passive tracers). For the thermal convection code, temperature is advected using the semi-lagrangian method, however, passive tracers can also be used to advect composition (density and viscosity) or temperature (so far, only the advection of the absolute temperature is implemented, though) if necessary (e.g., for the falling block or rigid body rotation benchmark, a growth rate of a Rayleigh-Taylor instability test, and the viscous inclusion benchmark). 
   
   The *MasterFile.m* in the main directory (*FDCSGm/*) contains a list and description of all constants, parameters, and variables, as well as a general structure of the code to solve the equations. Some variables can be removed if not needed, otherwise they need to be defined as *‘none’*. Within the code there are multiple benchmarks included to test the accuracy and different discretization methods for the energy, advection, and stokes equations. The thermal convection code is based on the routines used in those benchmarks. 

# Energy equation
   The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal within a closed system. In terms of a geodynamical problem, energy can be described by the temperature, which is transported mainly through *conductive* and *convective* processes, such that a general energy equation is defined as follows (assuming only radioactive heat sources):

$(\frac{\partial E}{\partial t} + \overrightarrow{v} \cdot \nabla E) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H$,&emsp;&emsp;&emsp;(1)

where the energy is described as $E=c_{p} \rho T$, *c<sub>p</sub>* is the specific heat capacity [J/kg/K], *ρ* is a reference density [kg/m<sup>3</sup>], *T* is the temperature [K], *t* is the time [s], $\overrightarrow{v}$ is the velocity vector [m/s], *q<sub>i</sub>* is the heat flux in direction of *i*  [W/m<sup>2</sup>], *∂/∂xi* is a directional derivative in direction of *i*, and *H* the heat production rate per mass [W/kg]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as follows: 

$\overrightarrow{q} = - k \nabla T$,&emsp;&emsp;&emsp;(2)

where *k* is the thermal conductivity [W/m/K]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation* in an Eulerian form can then be written as: 

$\rho c_p (\frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T) = -\frac{\partial q_i}{\partial x_i} + \rho H$,&emsp;&emsp;&emsp;(3) 
 
   This form of the temperature equation describes the variation of temperature due to a *conductive* (right hand side of the equation) and *convective* (left hand side of the equation) process. For a matter of simplicity, one can consider those terms in a separate manner  to solve the energy equation. Thus, I first focus on the conductive part of equation (3) and the equation simplifies (in 2-D) to (neglecting the convective term):

$\rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} -\frac{\partial q_z}{\partial z} + \rho H$,&emsp;&emsp;&emsp;(4)

or including Fourier’s law (assuming variable thermal parameters):

$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial z} k \frac{\partial T}{\partial z} + \rho H$.&emsp;&emsp;&emsp;(5) 

Assuming that the thermal parameters are constant, equation (5) simplifies to: 

$\frac{\partial T}{\partial t} = \kappa (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2}) + \frac{Q}{\rho c_p}$,&emsp;&emsp;&emsp;(6)
  
where κ is the thermal diffusivity [m<sup>2</sup>/s] and $Q=\rho H$ is the heat production rate per volume [W/m<sup>3</sup>]. Equation (6) is a *parabolic partial differential equation* which can be solve numerically in different manners, assuming initial and boundary conditions are defined. 

   Here I would like to discuss a simple, but effective finite difference method to discretize and solve the equation, that is the forward in time and centered in space (FTCS) method in an explicit manner. The advantage of an explicit description is that it is simple to derive and rather fast computationally, however, it is only numerically stable as long as the *heat diffusive stability criterion* is fulfilled. The stability criterion can be determined by a simple *Von Neumann* stability analysis, which analysis the growth of an eigenmode perturbation for a certain finite difference approach. In case of an explicit 2-D finite difference approach, the heat diffusive stability criterion is defined as $\Delta t < \frac{\Delta x^2}{3 \kappa}$, that is the time step is limited by the model’s resolution. Within the code different discretization methods can be chosen to solve the diffusive part of the temperature conservation equation (i.e., implicit, CNV, ADI). In the following, I would like to focus on the explicit method as an example of how to discretize the equation. For the other discretization methods check the *DiffusionProblem* directory. 

   Using a general finite difference discretization method to approximate the partial derivatives from equation (6), the temperature equation can be written as:

$\frac{T_{i,j}^{n+1} - T_{i,j}^{n} }{\Delta t} = \kappa (\frac{T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}}{(\Delta x)^2} + \frac{T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}}{(\Delta x)^2}) + \frac{Q_{i,j}^n}{\rho c_p}$&emsp;&emsp;&emsp;(7)

where *i*,*j* are the vertical and horizontal indices of the numerical finite difference grid, *n* is the time step index, Δ*t* is the time step, and Δ*x* and Δ*z* are the widths of the grid in horizontal and vertical direction. This equation contains know and unknow parameters and one can rearrange them to solve the equation for the unknowns as:

$T_{i,j}^{n+1} = T_{i,j}^{n} + s_x(T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}) + s_z(T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}) + \frac{Q_{i,j}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(8)

where $s_x = \frac{\kappa \Delta t}{(\Delta x)^2}$ and $s_z = \frac{\kappa \Delta t}{(\Delta z)^2}$. Equation (8) can be solved *iteratively* for every inner grid point assuming an initial condition is defined (multiple initial conditions can be set in the code). 

   Different conditions can be set for the boundaries of our model domain. Here, I focus on two fundamental conditions, the *Dirichlet* and *Neumann* boundary conditions. The Dirichlet boundary condition simply defines a constant temperature along the boundary, such that the temperature, for example, along the *left* boundary can be defined as:

$T_{i,j=1} = T_{left}$, for all i. &emsp;&emsp;&emsp;(9)

   The same applies to all other boundaries (*right*, *top*, and *bottom*). The Neumann boundary condition defines that the variation of a certain parameter does not change across the boundary, that is the temperature across the boundary or thermal heat flux *q* through the boundary. A sophisticated method to describe the heat flux across the boundary using finite differences is assuming the existence of so-called ghost nodes outside of our model domain, that are nodes which one can use to describe the flux across a boundary, but which do not actually exist. Therefore, one first needs to define the variation of the temperature across the boundary (e.g., the *left*, for *i* = 2,*nz*-1) as:

 $\frac{\partial T}{\partial x}|_{i,j} = c_l$,&emsp;&emsp;&emsp;(10)

or using finite differences: 

$\frac{T_{i,2} - T_{i,0}}{2 \Delta x} = c_l$,&emsp;&emsp;&emsp;(11)
 
where *c<sub>l</sub>* is a constant defining the flux across the boundary and *T<sub>i,0</sub>* are the ghost nodes outside the left boundary. Now, one can solve for an expression of the temperature at the ghost nodes which fulfills the condition of equation (10) as:

$T_{i,0} + T_{i,2} - 2 \Delta x c_l$.&emsp;&emsp;&emsp;(12)

Assuming equation (8) is also valid along the left boundary nodes and using the condition for the ghost nodes outside the numerical domain, that is one assumes Neumann boundary conditions, one can rewrite equation (8) for the left boundary nodes as follows:

$T_{i,1}^{n+1} = T_{i,1}^{n} + s_x(2T_{i,2}^{n} - 2(T_{i,1}^{n} + \Delta x c_l)) + s_z(T_{i+1,1}^{n} - 2T_{i,1}^{n} + T_{i-1,1}^{n}) + \frac{Q_{i,1}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(13)
 
The same applies to the other boundaries. Caution needs to be taken at the corners of the model. 

   The routines for the explicit, implicit, and ADI discretization methods are available in a dimensional and non-dimensional form in the *DiffusionProblem* directory. For more details on the scaling of the parameters see the scaling [section below](https://github.com/LukasFuchs/FDCSGm/blob/main/README.md#scaling-and-equation-of-state) and [*/FDCSGm/ScaleParam*](https://github.com/LukasFuchs/FDCSGm/tree/main/ScaleParam). Besides the routines to solve the diffusive part of the 2-D temperature conservation equation using different discretization methods, there are additional routines available within the code. There are two routines to calculate 1-D temperature profiles for a oceanic and continental lithosphere based on an 1-D implicit solver for constant and variable thermal parameters. A 1- and 2-D Poisson solver for a constant and variable *k* is available, to solve the steady-state diffusive equation (i.e. $\frac{\partial T}{\partial t} = 0$). For more details, also on the different discretization schemes, see [*/FDCSGm/DiffusionProblem/*](https://github.com/LukasFuchs/FDCSGm/tree/main/DiffusionProblem).

# Stokes equation 
   On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one might consider, i.e., inertia, surface, and volumetric forces. A common equation to describes such motion is given by: 

$\rho \frac{D \overrightarrow{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \boldsymbol{\rho}$,&emsp;&emsp;&emsp;(14)

where *ρ* is the density [kg/m<sup>3</sup>], $\overrightarrow{v}$ is the velocity vector [m/s], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [Pa], $\boldsymbol{g}$ is the gravitational acceleration [m/s<sup>2</sup>], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial}{\partial t} + \overrightarrow{v} \cdot \nabla$. The *Cauchy stress tensor* is given by: 

$\boldsymbol{\sigma} = -\nabla P + \nabla \cdot \boldsymbol{\tau}$,&emsp;&emsp;&emsp;(15)

where *P* is the total pressure (*P = P<sub>dynamic</sub> + P<sub>hydrostatic</sub>*) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. In Eulerian form, equation (14) is given by (Navier-Stokes equation):

$\rho (\frac{\partial v_{i}}{\partial x} + v_{j}\frac{v_{i}}{\partial x_{j}}) = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(16)

where summation over repeated indices is intended. 

   To solve equation (16), one needs to define a rheology which, for a purely viscous medium, can be given by a constitutive relationship between stress and strain rate in the form of, e.g.:

$\tau_{ij} = 2 \eta \dot{\varepsilon}_{ij}$,&emsp;&emsp;&emsp;(17)

where $\eta$ is the dynamic viscosity [Pa s] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* [1/s] and given by: 

$\dot{\varepsilon}_{ij} = \frac{1}{2} (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i})$.&emsp;&emsp;&emsp;(18)

   Assuming that the inertia forces are negligible in comparison to the gravitational forces, one can further simplify equation (16) to:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(19)

or in the form of the unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\partial}{\partial x_j} \eta (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i}) + \rho g_{i}$,&emsp;&emsp;&emsp;(20)

Assuming  constant viscosity equation (20) simplifies further to (Stokes equation): 

$0 = -\frac{\partial P}{\partial x_{i}} + \eta (\frac{\partial^2 v_i}{\partial x_j^2} + \frac{\partial^2 v_j}{\partial x_i^2}) + \rho g_{i}$.&emsp;&emsp;&emsp;(21)

   Equation (20) provides us two equations for our three unknowns. Thus, one needs to also consider the mass conservation equation, where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$\frac{\partial v_i}{\partial x_i} = 0$.&emsp;&emsp;&emsp;(22)

Equations (20) and (22) enable us to solve for the three unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*. 

   The conservation equations of momentum and mass are solved properly in two dimensions (*x* and *z*, *z* pointing downwards in negative direction) using a staggered finite difference grid, based on the description by *Gerya* (2009), where the viscosity and the density are defined on the regular grid points, the horizontal and vertical velocity in between the regular grid points, and the pressure within a finite difference cell. A staggered grid enables the conservation of the stress between adjacent grid points and one can solve equations (20)-(22) for the unknows. For more details on the discretization and the solving equations see *Gerya* (2009). 

# Advection equation 
   In  case the material is not moving, one can solve the energy equation only for the diffusive part (e.g., an intrusion problem or a non-deformation lithosphere). Generally, however, the material is moving and certain properties need to be advected with the flow (e.g., the temperature). Thermal mantle convection is a perfect example on how to transport heat with both diffusion (especiall in the thermal boundary layers) and advection (mainly within the interior). 
   
   The energy equation can be solved simultaneously with the diffusive and convective part using different discretization methods (interestingly, *FTCS* is stable with some numerical diffusion). However, for the sake of simplicity and a more conveniant way to teach both mechanisms (at least in my opinion), I do prefer, so far, the operator-splitting method, that is, I first solve for the convective part of the energy equation, followed by the conductive part. The conducitve part can be solved by the different discretization methods as described above and the convective part by the e.g., upwind, semi-lagrangian, or passive tracer method. Advection of a certain property *f* can then generally be described in an Eulerian reference frame as followed: 

$\frac{\partial f}{\partial t} = - \overrightarrow{v} \cdot \nabla f$.&emsp;&emsp;&emsp;(23)

For more details on how to discretize and sovle the advection equation see [*/FDCSGm/AdvectionProblem/*](https://github.com/LukasFuchs/FDCSGm/tree/main/AdvectionProblem). 

   For the thermal convection code I do prefer, so far, the semi-lagrangian method. However, I do only advect the total temperature and not the increments, so far. In case of a temperature independent material (e.g., a Rayleigh-Taylor Instability), advection with passive tracers, wich advect different properties (so far, density and viscosity), should be used (this partly works with a semi-lagrangian advection method, too). The tracers are advected using fourth order Runge-Kutta. For more details and examples see the [benchmark directory](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark).

# Scaling and equation of state
   To better compare different kinds of thermal convection, one can scale the governing equations. The equation can be scaled by the following reference parameters: 

$h_{sc} = h$,&emsp;&emsp;&emsp;(24)

$\eta_{sc} = \eta_0$,&emsp;&emsp;&emsp;(25)

$t_{sc} = \frac{h^2}{\kappa}$,&emsp;&emsp;&emsp;(26)

$v_{sc} = \frac{\kappa}{h}$,&emsp;&emsp;&emsp;(27)

$\tau_{sc} = \frac{\eta_0 \kappa}{h^2}$,&emsp;&emsp;&emsp;(28)

$T_{sc} = \Delta T$,&emsp;&emsp;&emsp;(29)

$H_{sc} = \frac{c_p \Delta T \kappa}{h^2}$,&emsp;&emsp;&emsp;(30)

where the subscript *sc* stands for the scaling parameters, and *h*, *η<sub>0</sub>*, *t*, *κ*, *v*, *τ*, *T*, *H*, *c<sub>p</sub>*, are the height, the reference viscosity, the time, the thermal diffusivity, the velocity, the stress, the temperature, the heat generation source, and the specific heat capacit, respectively. 

   The buoyance term on the right-hand side of equation (19), that is the density term which is temperature dependent (and pressure, but I do neglect this effect here so far), can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$\rho = \rho_0 (1-\alpha T)$,&emsp;&emsp;&emsp;(31)

where *ρ<sub>0</sub>* is the reference density and *α* the thermal expansion coefficient [1/K]. 

For the given scaling parameters, the non-dimensional governing equations are given as (assuming a constant viscosity; scaling for a variable viscosity is applicable in the same way):

$\frac{\partial v_x}{\partial x} + \frac{\partial v_z}{\partial z} = 0$,&emsp;&emsp;&emsp;(32)

$\frac{\partial T}{\partial t} + v_x \frac{\partial T}{\partial x} + v_z \frac{\partial T}{\partial z} = (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2} + H)$,&emsp;&emsp;&emsp;(33)

$-\frac{\partial P}{\partial x} + \eta \frac{\partial^2 v_x}{\partial x^2} + \eta \frac{\partial^2 v_x}{\partial z^2} = 0$,&emsp;&emsp;&emsp;(34)

$-\frac{\partial P}{\partial z} + \eta \frac{\partial^2 v_z}{\partial z^2} + \eta \frac{\partial^2 v_z}{\partial x^2} - RaT = 0$,&emsp;&emsp;&emsp;(35)

where *Ra* is the so-called thermal *Rayleigh number* and *P* the *dynamic pressure*. In case of a basally heated thermal convection, the convective vigor is defined by the Rayleigh number, which describes a relationship between heat transported by buoyancy and conduction, and the effect of the layers thickness and bulk viscosity.

# References
*Gerya, T. (2009). Introduction to numerical geodynamic modelling.*
