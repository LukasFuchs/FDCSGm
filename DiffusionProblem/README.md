# General Information

&emsp;This directory contains several rountines to solve the diffusive part of the *temperature conservation equation* (1- and 2-D, steady state and time-dependent) using different numerical discretization methods. The routines are avaible in a dimensional or non-dimensional form (files ending with *Sc.m; see [*FDCSGm/ScaleParam*](https://github.com/LukasFuchs/FDCSGm/tree/main/ScaleParam) regarding the scaling; so far the 1-D routines are only availabe in a dimensional version!). 

-------------

## Energy equation
&emsp;The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal. In terms of a geodynamical problem, energy can be described by temperature, which is transported mainly through *conductive* and *convective* processes, such that a general energy equation is defined as followed (assuming only radioactive heat sources):

$(\frac{\partial E}{\partial t} + \overrightarrow{v} \cdot \nabla E) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H$,&emsp;&emsp;&emsp;(1)

where the energy is described as $E=c_{p} \rho T$, and *c<sub>p</sub>* is the specific heat capacity [J/kg/K], *ρ* is a reference density [kg/m<sup>3</sup>], *T* is the temperature [K], *t* is the time [s], $\overrightarrow{v}$ is the velocity vector [m/s], *q<sub>i</sub>* is the heat flux in direction of *i*  [W/m<sup>2</sup>], *∂/∂xi* is a directional derivative in direction of *i*, and *H* the heat production rate per mass [W/kg]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as followed: 

$\overrightarrow{q} = - k \nabla T$,&emsp;&emsp;&emsp;(2)

where *k* is the thermal conductivity [W/m/K]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation* in an Eulerian form can then be written as: 

$\rho c_p (\frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T) = -\frac{\partial q_i}{\partial x_i} + \rho H$,&emsp;&emsp;&emsp;(3) 
 
&emsp;This form of the temperature equation describes the variation of temperature due to a *conductive* (right hand side of the equation) and *convective* (left hand side of the equation) process. For a matter of simplicity, one can consider those terms in a separate manner  to solve the energy equation. Thus, I first focus on the conductive part of equation (3) and the equation simplifies (in 2-D) to (neglecting the convective term):

$\rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} -\frac{\partial q_z}{\partial z} + \rho H$,&emsp;&emsp;&emsp;(4)

or including Fourier’s law (assuming variable thermal parameters):

$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial z} k \frac{\partial T}{\partial z} + \rho H$.&emsp;&emsp;&emsp;(5) 

Assuming that the thermal parameters are constant, equation (5) simplifies to: 

$\frac{\partial T}{\partial t} = \kappa (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2}) + \frac{Q}{\rho c_p}$,&emsp;&emsp;&emsp;(6)
  
where $\kappa = k/\rho /c_p$ is the thermal diffusivity [m<sup>2</sup>/s] and $Q=\rho H$ is the heat production rate per volume [W/m<sup>3</sup>]. Equation (6) is a *parabolic partial differential equation* which can be solved numerically in different manners, assuming initial and boundary conditions are defined. 

&emsp;Here, I first would like to discuss a simple, but effective, finite difference method to discretize and solve the equation, that is the forward in time and centered in space (FTCS) method in an *explicit* manner. This finite difference scheme will converge to the exact solution for small $\Delta x$ and $\Delta t$. The advantage of an explicit description is that it is **simple** to derive and rather **fast** computationally. However, it is only numerically stable as long as the *heat diffusive stability criterion* is fulfilled. The stability criterion can be determined by a *Von Neumann* stability analysis, which analyzes the growth of an eigenmode perturbation for a certain finite difference approach. In case of an **explicit 2-D finite difference approach**, the heat diffusive stability criterion is defined as $\Delta t < \frac{\Delta x^2}{3 \kappa}$ (assuming equal grid spacing in *x*- and *z*-direction), that is the time step is limited by the model’s resolution. 

### Explicit, FTCS

&emsp;Using an FTCS, explicit finite difference scheme to approximate the partial derivatives from equation (6), the temperature equation can be written as:

$\frac{T_{i,j}^{n+1} - T_{i,j}^{n} }{\Delta t} = \kappa (\frac{T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}}{(\Delta x)^2} + \frac{T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}}{(\Delta z)^2}) + \frac{Q_{i,j}^n}{\rho c_p}$&emsp;&emsp;&emsp;(7)

where *i* and *j* are the vertical and horizontal indices of the numerical finite difference grid, *n* is the time step index, Δ*t* is the time step, and Δ*x* and Δ*z* are the widths of the grid in horizontal and vertical direction, respectively. This equation contains know and unknow parameters and one can rearrange them to solve the equation for the unknowns as:

$T_{i,j}^{n+1} = T_{i,j}^{n} + s_x(T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}) + s_z(T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}) + \frac{Q_{i,j}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(8)

where $s_x = \frac{\kappa \Delta t}{(\Delta x)^2}$ and $s_z = \frac{\kappa \Delta t}{(\Delta z)^2}$. Equation (8) can be solved *iteratively* for every inner grid point assuming an initial condition is defined (multiple initial conditions can be set in the code, check [*SetUpInitialConditions.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/SetUpInitialConditions.m)). For more details on how this is implemented in MATLAB see [*SolveDiff2Dimplicit.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dexplicit.m).

&emsp;For the boundaries of our model domain, different thermal conditions can be set. Here, I focus on two fundamental conditions, the *Dirichlet* and *Neumann* boundary conditions. The Dirichlet boundary condition defines a constant temperature along the boundary, such that the temperature, for example, along the *left* boundary can be defined as:

$T_{i,j=1} = T_{left}$, for all i. &emsp;&emsp;&emsp;(9)

The same applies to all other boundaries (*right*, *top*, and *bottom*). 

&emsp;The Neumann boundary condition defines that the variation of a certain parameter does not change across the boundary, that is, for example, the temperature across the boundary or thermal heat flux *q* through the boundary. A sophisticated method to describe the heat flux across the boundary using finite differences is assuming the existence of so-called ghost nodes, or fictitious nodes, outside of the model domain, that are nodes one can use to describe the flux across a boundary, but do not actually exist. Therefore, one first needs to define the variation of temperature across the boundary (e.g., the *left*, for *j* = 1 and *i* = 2,*nz*-1) as:

 $\frac{\partial T}{\partial x}\vert_{i,1} = c_{left}$,&emsp;&emsp;&emsp;(10)

or using finite differences: 

$\frac{T_{i,2} - T_{i,0}}{2 \Delta x} = c_{left}$,&emsp;&emsp;&emsp;(11)
 
where *c<sub>left</sub>* is a constant defining the flux through the boundary and *T<sub>i,0</sub>* are the ghost nodes outside the left boundary. Now, one can solve for an expression of the temperature at the ghost nodes which fulfils the condition of equation (10) as:

$T_{i,0} = T_{i,2} - 2 \Delta x c_l$.&emsp;&emsp;&emsp;(12)

&emsp;Considering that equation (8) is also valid along the left boundary nodes assuming Neuman conditions, one can rewrite equation (8) for the left boundary nodes using the condition for the ghost nodes outside the numerical domain (equation (12)) as followed:

$T_{i,1}^{n+1} = T_{i,1}^{n} + s_x(2T_{i,2}^{n} - 2(T_{i,1}^{n} + \Delta x c_l)) + s_z(T_{i+1,1}^{n} - 2T_{i,1}^{n} + T_{i-1,1}^{n}) + \frac{Q_{i,1}^n \Delta t}{\rho c_p}$,&emsp;&emsp;&emsp;(13)
 
The same applies for the other boundaries. Caution needs to be taken at the corners of the model. 

## Discretization Methods

&emsp; Within the code different discretization methods can be chosen to solve the diffusive part of the temperature conservation equation (i.e., *explicit*, *implicit*, *CNA*, *ADI*). While the explicit scheme seems to be the most accurate finite difference scheme for the time-dependent diffusion equation I implemented (see [Gaussian Diffusion Benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/GaussDiffusion)), the dependency of the time stepping on the grid resolution might become problematic in models with a high resolution and could significantly slow down the model calculations (that is, smaller time steps are necessary for higher resolutions!). In the following, I will present some well-know alternative discretization methods for the diffusive part of the temperature equation, which can help to resolve this issue, and briefly discuss their advantages and disadvantages. All discretization methods can be used in the [thermal convection code](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems) and the [Blankenbach Benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach) and are generally available to chose in the code **FDCSGm**. A more detailed analysis on the accuracy of each discretization scheme and the effect of the grid resolution is given in the [Gaussian Diffusion Benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/GaussDiffusion). 

### Implicit, FTCS

&emsp;The fully implicit finite difference scheme is unconditionally stable and one can use time steps larger than the diffusion time criterion. In 2-D, the diffusion equation is then given as: 

$\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t}=\kappa (\frac{T_{i,j+1}^{n+1}-2T_{i,j}^{n+1}+T_{i,j-1}^{n+1}}{\Delta x^2} + \frac{T_{i+1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i-1,j}^{n+1}}{\Delta z^2}) + \frac{Q_{i,j}^n}{\rho c_p}$, &emsp;&emsp;&emsp; (14)

where *n* is the current and *n+1* the next time step, $\Delta t$ is the time step length, $\Delta x$ and $\Delta z$ are the horizontal and vertical grid spacing, and *i* and *j* are the vertical and horizontal indices, respectively. Rearranging equation (14) into known and unknown variables, one obtains a linear system of equations in the form of: 

$-s_zT_{i-1,j}^{n+1}-s_xT_{i,j-1}^{n+1}+(1+2s_z+2s_x)T_{i,j}^{n+1}-s_xT_{i,j+1}^{n+1}-s_zT_{i+1,j}^{n+1}=T_{i,j}^n+\frac{Q_{i,j}^n \Delta t}{\rho c_p}$, &emsp;&emsp;&emsp; (15)

where $s_x=\frac{\kappa \Delta t}{\Delta x^2}$ and $s_z=\frac{\kappa \Delta t}{\Delta z^2}$. This is a linear system of equation in the form of $\boldsymbol{A}\cdot x = rhs$, where ***A*** is a coefficient matrix (here with five non-zero diagonals), *x* the unknown vector, and *rhs* the known right-hand side. The main advantage of the implicit method is that there are no restrictions on the time step, but this does not mean that it is accurate. Taking too large time steps may result in an inaccurate solution for features with small spatial scales. Here, I do solve the linear system of equations in MATLAB with a right-array divison (but other solver method are applicable, too). 

&emsp;As promising as the implicit method is, the band-width of the coefficient matrix could be problematic in models with a very high resoltuion. In addition, one needs to be careful in setting up the coefficient matrix and the corresponding unknown and *rhs* vectors (especially in 2-D and even more in 3-D). Here, one needs to do a proper indexing of the regular grid points by using a global index variable, which goes from 1 to *nx* x *nz*. This global index replaces the local *i* and *j* indices in equation (14), such that the spatial derivatives at a point *i*, *j* are given as: 

$\frac{\partial^2 T}{\partial x^2} \vert_{i,j} = \frac{T_{(i-1)nx+j+1}-2T_{(i-1)nx+j} + T_{(i-1)nx+j+1}}{\Delta x^2}$,&emsp;&emsp;&emsp; (16)

$\frac{\partial^2 T}{\partial z^2} \vert_{i,j} = \frac{T_{i \cdot nx+j}-2T_{(i-1)nx+j} + T_{(i-2)nx+j}}{\Delta z^2}$.&emsp;&emsp;&emsp; (17)

&emsp;Equations (16) and (17) in combination with equation (14) enables to setup the linear system of equations to solve for the temperature at the new time step. For more details on how this is implemented in MATLAB see [*SolveDiff2Dimplicit.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dimplicit.m).

&emsp;Similar to the *explicit* discretization, one can use fictitious grid points outside the model domain to define Neumann boundary conditions, such that equation (15) results in (e.g., for the left boundary, i.e., *j* = 1): 

$-s_zT_{i-1,1}^{n+1}+(1+2s_z+2s_x)T_{i,1}^{n+1}-2s_xT_{i,2}^{n+1}-s_zT_{i+1,1}^{n+1}=T_{i,1}^n+2s_x\Delta xc_{left}+\frac{Q_{i,1}^n \Delta t}{\rho c_p}$, &emsp;&emsp;&emsp; (18)

where *c<sub>left</sub>* needs to fulfil the following condition at the left boundary: 

$\frac{\partial T}{\partial x} = c_{left} = \frac{T_{i,2}-T_{i,0}}{2\Delta x}$. &emsp;&emsp;&emsp; (19) 

The same applies for the remaining boundaries. 

### Cranck-Nicolson approach (CNA)

&emsp;The fully implicit method works well, but is only first order accurate in time. A way to modify this is to employ a Crank-Nicolson time step discretization, which is implicit and thus second order accurate in time. In 2-D, equation (6) is then discritized as: 

$\frac{T_{i,j}^{n+1} - T_{i,j}^{n}}{\Delta t} = \frac{\kappa}{2}\frac{(T_{i,j+1}^{n+1}-2T_{i,j}^{n+1}+T_{i,j-1}^{n+1})+(T_{i,j+1}^{n}-2T_{i,j}^{n}+T_{i,j-1}^{n})}{\Delta x^2} + \frac{\kappa}{2}\frac{(T_{i+1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i-1,j}^{n+1})+(T_{i+1,j}^{n}-2T_{i,j}^{n}+T_{i-1,j}^{n})}{\Delta z^2}$. &emsp; &emsp; &emsp;  (20)

&emsp;However, the band-width of the coefficient matrix increases as in the fully implicit case. Thus, the method becomes memory intensiv for models with a high resoltuion. For more details on how this is implemented in MATLAB, see [*SolveDiff2DCNV.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DCNV.m).


### Alternating Direction Implicit (ADI)

&emsp; Within the ADI method, one basically decomposes the calculation of one time step into two half-steps. For the first step $(n \ \text{to}\ n+1/2)$, the *x*-direction is solved explicitly and the *z*-direction implicitly and, for the second step $(n+1/2 \ \text{to}\ n+1)$, the *x*-direction is solved implicitly and the *z*-direction explicitly. The advantage of the ADI method is that the equation in each step has a simpler structure and can be solved more efficiently (e.g., with the tridiagonal matrix algorithm). Equation (6) for each half-step is then given as: 

$\frac{T_{i,j}^{n+1/2}-T_{i,j}^n}{\Delta t/2}=\kappa (\frac{T_{i+1,j}^n-2T_{i,j}^n+T_{i-1,j}^n}{\Delta x^2} + \frac{T_{i,j+1}^{n+1/2}-2T_{i,j}^{n+1/2}+T_{i,j-1}^{n+1/2}}{\Delta z^2})$; &emsp; &emsp; &emsp; (21)

$\frac{T_{i,j}^{n+1}-T_{i,j}^{n+1/2}}{\Delta t/2}=\kappa (\frac{T_{i+1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i-1,j}^{n+1}}{\Delta x^2} + \frac{T_{i,j+1}^{n+1/2}-2T_{i,j}^{n+1/2}+T_{i,j-1}^{n+1/2}}{\Delta z^2})$; &emsp; &emsp; &emsp; (22)

For more details on how this is implemented in MATLAB, see [*SolveDiff2DADI.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DADI.m).

------------------------------------------------------------------------------------------------------------------------------------

&emsp;The routines for the *explicit*, *implicit*, and *ADI* discretization methods are available in a **dimensional** and **non-dimensional** form. For more details on the scaling of the parameters see [*/FDCSGm/ScaleParam*](https://github.com/LukasFuchs/FDCSGm/tree/main/ScaleParam). 

## Steady State Solution

&emsp;So far, I only implemented variable thermal parameters in the 1-D and the 2-D steady state solutions. In steady state, one assumes that the temperature does not vary over time (i.e., $\frac{\partial T}{\partial t}=0$) and the temperature equation simplifies to an *elliptic partial differential* equation (i.e., the *Poission equation*).  

### Poisson solution, constant *k*

&emsp;For constant thermal parameters the diffusive temperature equation is given by (in 2-D): 

$0=(\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial z^2}) + \frac{Q}{k}$.&emsp;&emsp;&emsp; (23)

&emsp;For the approximation of the spatial partial derivatives with finite difference expressions, I chose a central finite difference and equation (23) is then given as: 

$0=(\frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{\Delta x^2} + \frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{\Delta z^2}) + \frac{Q}{k}$,&emsp;&emsp;&emsp; (24)

where *i* and *j* are the indices for the *z*- and *x*-direction, respectively. Now, one can rearrange the equation by known (*Q*, *k*) and unknown (*T*) variables, wich results in a linear system of equations in the form of: 

$s_zT_{i-1,j}+s_xT_{i,j-1}-2(s_x+s_z)T_{i,j}+s_xT_{i,j+1}+s_zT_{i+1,j}=-\frac{Q}{k}$, &emsp;&emsp;&emsp; (25)

where $s_x = \frac{1}{\Delta x^2}$ and $s_z = \frac{1}{\Delta z^2}$. For more details on how this is implemented in MATLAB see [*SolvePoisson2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolvePoisson2D.m).

&emsp;Here, one can again assume Dirichlet or Neumann boundary conditions, where for Dirichlet the temperature along the boundary is kept constant and for Neumann equation (25) is given by (e.g., along the left boundary): 

$s_zT_{i-1,j}-2(s_x+s_z)T_{i,j}+2s_xT_{i,j+1}+s_zT_{i+1,j} = -\frac{Q}{k}+\frac{2c_{left}}{\Delta x}$.&emsp;&emsp;&emsp; (26)

&emsp;Here, one again uses imaginary points outside of the model domain to define the flux boundary conditions. The same applies for the remaining boundaries. 

### Poisson solution, variable *k*

&emsp; For variable thermal parameters the steady-state temperature equation is given by (in 2-D): 

$0=\frac{\partial}{\partial x}(k\frac{\partial T}{\partial x})+\frac{\partial}{\partial z}(k\frac{\partial T}{\partial z})+Q$.&emsp;&emsp;&emsp; (27)

&emsp;To properly solve equation (27), one needs to apply a conservative finite difference scheme, such that the heat flux $(q_i = k\frac{\partial T}{\partial i})$ is defined between the regular grid nodes (e.g., points A,B,C, and D) and the temperature on the regular grid nodes. Therefore, one needs to average the conductivity and equation (27) results in a linear system of equations in the form of: 

$\frac{k_D}{\Delta z^2}T_{i-1,j}+\frac{k_A}{\Delta x^2}T_{i,j-1}+(-\frac{1}{\Delta x^2}(k_A+k_B)-\frac{1}{\Delta z^2}(k_c+k_D))T_{i,j}+\frac{k_B}{\Delta x^2}T_{i,j+1}+\frac{k_C}{\Delta z^2}T_{i+1,j} = -Q_{i,j}$, &emsp;&emsp;&emsp; (28)

where $k_A = \frac{k_{i,j+1}+k_{i,j}}{2}$, $k_B = \frac{k_{i,j-1}+k_{i,j}}{2}$, $k_C = \frac{k_{i+1,j}+k_{i,j}}{2}$, and $k_D = \frac{k_{i-1,j}+k_{i,j}}{2}$. For more details on how this is implemented in MATLAB see [*SolvePoisson2Dvaryk.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolvePoisson2Dvaryk.m).

---------------------

## Example: Geotherms

### 1-D Geotherms

&emsp;The 1-D temperature profile is calculated by solving the diffusive part of the 1-D temperature conservation equation (so far only with a radiogenic heat source) for variable thermal parameters with a proper conserving finite difference scheme. That is, the heat flow is calculated on the centered and the remaining parameters on the regular grid points, respectively. The discretization scheme for variable thermal parameters is picked to solve for a temperature profile of a continental lithosphere with upper, lower crust, and mantle. The 1-D temperature equation is given by: 

$\rho c_{p} \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}(k \frac{\partial T}{\partial z}) + \rho H$, &emsp; &emsp; &emsp; (29)

where $\rho, c_{p}, T, t, k, H, z$ are the density [kg/m<sup>3</sup>], the specific heat capacity [J/kg/K], the temperature [K], the time [s], the thermal conductivity [W/m/K], the heat generation rate per mass [W/kg], and the depth [m] respectively. For values and references of the given thermal parameters see [*OceanicGeotherm_1D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/OceanicGeotherm_1D.m) and [*ContinentalGeotherm_1D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/ContinentalGeotherm_1D.m).

&emsp;Here, a proper conservative finite difference scheme means that the heat flux is calculated on the centered grid points (A, B, etc.). The 1-D vertical heat flux is given by the Fourier’s law:

$q_{z} = -k \frac{\partial T}{\partial z}$. &emsp; &emsp; &emsp; (30)

#### ***Solving the equation***

&emsp;Following the discretization as described above, one needs to solve the following equation (in an [implicit finite difference formulation](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dimplicit_vary.m)):

$\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = -\frac{q_{z,B}^{n+1} - q_{z,A}^{n+1} }{\Delta z} + \rho_j H_j$, &emsp;&emsp;&emsp; (31)

$\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = \frac{ k_{B} \frac{T_{j+1}^{n+1} - T_{j}^{n+1}}{\Delta z} - k_{A} \frac{T_{j}^{n+1} - T_{j-1}^{n+1}}{\Delta z} }{\Delta z} + \rho_j H_j$. &emsp;&emsp;&emsp; (32)

Sorting the variables (known variables on the right-hand side, unknown on the left-hand side): 

$\frac{k_{B}}{\Delta z^2} T_{j+1}^{n+1} - \frac{k_{B}}{\Delta z^2} T_{j}^{n+1} - \frac{k_{A}}{\Delta z^2} T_{j}^{n+1} + \frac{k_{A}}{\Delta z^2} T_{j-1}^{n+1} = \frac{\rho_j c_{p,j}}{\Delta t} T_{j}^{n+1} - \frac{\rho_j c_{p,j}}{\Delta t} T_{j}^{n} - \rho_j H_j$, &emsp;&emsp;&emsp;(33)

$\frac{k_{B}}{\Delta z^2} \frac{\Delta t}{\rho_j c_{p,j}} T_{j+1}^{n+1} - \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{B} + k_{A}}{\Delta z^2}) T_{j}^{n+1} - T_{j}^{n+1} + \frac{k_{A}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}} T_{j-1}^{n+1} = -T_{j}^{n} - \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(34)

$-\frac{k_{B}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}T_{j+1}^{n+1} + (1 + \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{B} + k_{A}}{\Delta z^2}))T_{j}^{n+1} - \frac{k_{A}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}} T_{j-1}^{n+1} = T_j^n + \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(35)

$aT_{j-1}^{n+1} + bT_{j}^{n+1} + cT_{j+1}^{n+1} = T_j^n + \frac{H_j \Delta t}{c_{p,j}}$, &emsp;&emsp;&emsp;(36)

with

$a = -\frac{k_{A}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}$

$b = (1 + \frac{\Delta t}{\rho_j c_{p,j}} (\frac{k_{B} + k_{A}}{\Delta z^2}))$, &emsp;&emsp;&emsp; (37) 

$c = - \frac{k_{B}}{\Delta z^2}\frac{\Delta t}{\rho_j c_{p,j}}$

and

$k_{A} = \frac{k_{j-1} + k_{j}}{2}$

$k_{B} = \frac{k_j + k_{j+1}}{2}$. &emsp;&emsp;&emsp; (38)

An [*explicit*](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dexplicit_vary.m) solver for a 1-D thermal profile with variable thermal parameters and a radiogenic heat source is also available.

#### ***Thermal boundary conditions***

The thermal boundary conditions are defined as: 

1. **Constant temperature** (*Dirichlet*)<br>
The temperature at the top or bottom can just be set as constant to *T<sub>top</sub>* or *T<sub>bot</sub>*, respectively.
      
2. **Constant temperature gradient** (*Neumann*)<br>
The gradient of temperature (and thus the vertical heat flux) can be defined using so called ghost nodes at the top and the bottom of the profile. Therefore, we define the condition at the top and bottom as:

   $\frac{\partial T}{\partial z} \vert_{j=1} = c_{top} = \frac{T_2-T_0}{2\Delta z}$

   $\frac{\partial T}{\partial z} \vert_{j=nz} = c_{bottom} = \frac{T_{nz+1}-T_{nz-1}}{2\Delta z}$, &emsp;&emsp;&emsp; (39)

   Where *T<sub>0</sub>* and *T<sub>nz+1</sub>* are the ghost nodes for temperature at the top and bottom, respectively. The constants *c<sub>top</sub>* and *c<sub>bottom</sub>* are defined as:

   $c_{top,bottom} = -\frac{q_{top,bottom}}{2\Delta z}$
   
   Using these conditions, we can define formulations for the temperature at the ghost nodes as:

   $T_0 = T_2 - 2\Delta z c_{top}$

   $T_{nz+1} = T_{nz-1} + 2\Delta z c_{bottom}$

   Now one can solve equation (31) for the top and the bottom using the formulations of the temperature at the ghost nodes with equation (39), which results in:

   $\rho_{j=1}c_{p,j=1}\frac{\partial T}{\partial t}\vert_{j=1} = \frac{\partial}{\partial z}(k\frac{\partial T}{\partial z})$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{ k_A\frac{dT}{dz}\vert_A - k_B\frac{dT}{dz}\vert_B }{\Delta z}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{ k_A\frac{T_2^{n+1} - T_1^{n+1}}{\Delta z} - k_B\frac{T_1^{n+1} - T_0^{n+1}}{\Delta z} }{\Delta z}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{k_A}{\Delta z^2}T_2^{n+1} - \frac{k_A}{\Delta z^2}T_1^{n+1} - \frac{k_B}{\Delta z^2}T_1^{n+1} + \frac{k_B}{\Delta z^2}T_0^{n+1}$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = \frac{k_A}{\Delta z^2}T_2^{n+1} - \frac{k_A+k_B}{\Delta z^2}T_1^{n+1} + \frac{k_B}{\Delta z^2}(T_2^{n+1} - 2\Delta z c_{top})$

   $\rho_1c_{p,1}\frac{T_1^{n+1} - T_1^n}{\Delta t} = - \frac{k_A + k_B}{\Delta z^2}T_1^{n+1} + \frac{k_A+k_B}{\Delta z^2}T_2^{n+1} - \frac{k_B}{\Delta z}2c_{top}$

   $(1 + \frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2} (k_A+k_B) ) T_1^{n+1} - \frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2} (k_A+k_B) T_2^{n+1} = T_1^n - \frac{2\Delta t c_{top}}{ \rho_1 c_{p,1} \Delta z} k_B$

   $aT_1^{n+1}+bT_2^{n+1}=T_1^n+Q_{top}$, &emsp;&emsp;&emsp; (40)

   with

   $a = 1+\frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2}(k_A+k_B)$

   $b = -\frac{\Delta t}{\rho_1 c_{p,1} \Delta z^2}(k_A+k_B)$

   $Q_{top} = -\frac{2\Delta t c_{top}}{\rho_1 c_{p,1} \delta z}k_A$. &emsp;&emsp;&emsp; (41)

   Similar for the bottom boundary with:

   $aT_{nz}^{n+1}+bT_{nz-1}^{n+1}=T_{nz}^{n} + Q_{bottom}$

   $a = 1+\frac{\Delta t}{\rho_{nz} c_{p,nz} \Delta z^2}(k_A+k_B)$

   $b = -\frac{\Delta t}{\rho_{nz} c_{p,nz} \Delta z^2}(k_A+k_B)$

   $Q_{top} = -\frac{2\Delta t c_{bottom}}{\rho_{nz} c_{p,nz} \delta z}k_A$. &emsp;&emsp;&emsp; (42)

#### Oceanic Geotherms
![OLGT1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/fa77bb60-1314-4be4-8852-499b1a1be17c)<br>
***Figure 1. Oceanic Lithosphere.** LEFT: Temperature profile [K]  for an oceanic lithosphere of 60 Ma of age and constant thermal boundary conditions at the top and bottom. The blue line shows the initial temperature profile. The yellow dashed line shows the solution for a half-space cooling model. RIGHT: Heat flux [mW/m<sup>2</sup>] with depth. The parameters of this model are defined as the default values in the routine OceanicGeotherm.m.*

![OLGT2](https://github.com/LukasFuchs/FDCSGm/assets/25866942/638c6d66-7c96-4b6c-9924-4d00718504ec)<br>
***Figure 2. Oceanic Lithosphere II**. Same as Figure 1 but with constant heat flux boundary conditions qbottom =10 mW/m<sup>2</sup> and qtop = 90 mW/m<sup>2</sup>.*

#### Continental Geotherms
![CLGT1](https://github.com/LukasFuchs/FDCSGm/assets/25866942/a040e96d-2896-4632-ba0c-67dce616ea33)<br>
***Figure 3. Continental Lithosphere.** LEFT: Temperature profile for a continental lithosphere of 1000 Ma of age with constant upper and lower thermal boundary conditions. The blue line shows the initial condition, the red line shows the solution of equation (1), the yellow dashed line shows the solution of the time-independent heat equation (1-D poisson equation), and the magenta dashed line shows the solution of a 2D, staggered finite difference code. MIDDLE: Heat flux with depth. RIGHT: Thermal parameter for the lithosphere setup: thermal conductivity [k], specific heat [c<sub>p</sub>], density [ρ], and volumetric heat generation rate [Q].*

![CLGT2](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f6e8bc3f-dd9f-4e8b-b97a-64ff9ec1b151)<br>
***Figure 5. Continental Lithosphere II**. Same as Figure 3 but with constant upper and lower heat flux boundary conditions, q<sub>top</sub> = 40 mW/m<sup>2</sup> and q<sub>bottom</sub> = 10 mW/m<sup>2</sup>.*

# Directory Content
[ContinentalGeotherm_1D.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/ContinentalGeotherm_1D.m)<br>
&emsp;- Script to calculate the 1-D continental geotherm.
    
[Diff1D_stationary.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/Diff1D_stationary.m)<br>
&emsp;-> Script to calculate the 1-D steady state solution.

[Diffusion.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/Diffusion.m)<br>
&emsp;-> General function used in the thermal convection code to call all diffusion solver functions.

[OceanicGeotherm_1D.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/OceanicGeotherm_1D.m)<br>
&emsp;-> Script to calculate the 1-D oceanic geotherm.

[SolveDiff1Dexplicit_vary.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dexplicit_vary.m)<br>
&emsp;-> Function to solve the 1-D diffusive temperature equation with variable thermal parameters using an explicit scheme.

[SolveDiff1Dimplicit](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dimplicit.m)[(_vary).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff1Dimplicit_vary.m)<br>
&emsp;-> Functions to solve the 1-D diffusive temperature equation for variable or constant thermal parameters using an implicit scheme.

[SolveDiff2DADI](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DADI.m  )[(Sc).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DADISc.m)<br>
&emsp;-> Functions to solve the 2-D diffusive temperature equation for constant thermal parameters in a dimensional or non-dimensional (Sc) form using the ADI method.

[SolveDiff2DCNV.m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DCNV.m)<br>
&emsp;-> Function to solve the 2-D diffusive temperature equation for constant thermal parameters in a dimensional form using the Crank-Nicolson method.

[SolveDiff2Dexplicit](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dexplicit.m)[(Sc).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DexplicitSc.m)<br>
&emsp;-> Functions to solve the 2-D diffusive temperature equation for constant thermal parameters in a dimensional or non-dimensional (Sc) form using an explicit method.

[SolveDiff2Dimplicit](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dimplicit.m)[(Sc).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2DimplicitSc.m)<br>
&emsp;-> Functions to solve the 2-D diffusive temperature equation for constant thermal parameters in a dimensional or non-dimensional (Sc) form using an implicit method.

[SolveDiff2Dimplicit_opt](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dimplicit_opt.m)[(Sc).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolveDiff2Dimplicit_optSc.m)<br>
&emsp;-> Functions to solve the 2-D diffusive temperature equation for constant thermal parameters in a dimensional or non-dimensional (Sc) form using an optimized implicit method. See functions for more details.

[SolvePoisson2D](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolvePoisson2D.m)[(varyk).m](https://github.com/LukasFuchs/FDCSGm/blob/main/DiffusionProblem/SolvePoisson2Dvaryk.m)<br>
&emsp;-> Functions to solve the 2-D steady state diffusive temperature equation for constant and variable thermal parameters.

[comment]: <> (Needs, in detail:- discretization of explicit variable thermal parameters, at some point!)

