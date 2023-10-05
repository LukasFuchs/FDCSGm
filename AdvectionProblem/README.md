# General Information 

&emsp;This directory contains several files to calculate the advection of a certain property within a two-dimensional domain by the e.g., *upwind*, *semi-lagrangian*, or *passive tracer* method. In the [thermal convection](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems) code I do prefer, so far, the semi-lagrangian method. However, I only advect the **absolut** temperature and not its increments, so far. In case of a temperature independent material (e.g., a [Rayleigh-Taylor Instability](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RTI)), advection with passive tracers, wich advect different properties (so far, density and viscosity), should be used. For more details and examples see the [benchmark directory](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark). For the sake of simplicity, I do focus on the advection of the absolute temperature in the following, however, in theory, one can chose any kind of property for advection (with some conditions to keep in mind). 

## Advection equation 
&emsp;In case the material is not moving, one can solve the energy equation only for the diffusive part (e.g., an intrusion problem or a non-deforming lithosphere). Generally, however, the material is moving and certain properties need to be advected with the flow (e.g., the temperature). Thermal mantle convection is a perfect example on how to transport heat with both diffusion (especiall in the thermal boundary layers) and advection (mainly within the interior). 
   
&emsp;The energy equation can be solved simultaneously with the diffusive and convective part using different discretization methods (interestingly, *FTCS* is stable with some numerical diffusion, which is always unstable for pure advection). However, for the sake of simplicity and a more conveniant way to teach both mechanisms (at least in my opinion), I do prefer, so far, the operator-splitting method, that is, I first solve for the convective part of the energy equation, followed by the conductive part. The conducitve part can be solved by the different discretization methods as described in the [diffusion problem directory](https://github.com/LukasFuchs/FDCSGm/tree/main/DiffusionProblem) and the convective part by the e.g., *upwind*, *semi-lagrangian*, or *passive tracer* method. 

&emsp;In general, advection describes the transport of a property, here the temperature, from one point to another, where one can assume different reference frames for the given point of interest. If we assume a not moving reference frame (that is an *Eulerian* grid), the change in temperature at a certain point can be described by (i.e. the *eulerian advective transport equation*): 

$\frac{\partial T}{\partial t} = - \overrightarrow{v} \cdot \nabla T$,&emsp;&emsp;&emsp;(1)

or in a Lagrangian reference frame (along a moving point; i.e., the *substantive* derivative) as: 

$\frac{DT}{Dt}$, &emsp;&emsp;&emsp;(2)

where both are related by: 

$\frac{DT}{Dt} = \frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T$.&emsp;&emsp;&emsp;(3)

&emsp;For a Lagrangian reference point advection is given by a simple *ordinary differential equation* particle advection scheme, where changes in its coordinates are related with the material velocities as: 

$\frac{Dx_i}{Dt} = v_i$,&emsp;&emsp;&emsp; (4)

where *i* is the coordinate index and *x<sub>i</sub>* is a spatial coordinate. 

--------------------------

## Discretization Schemes

&emsp;As simple as the advection problem sounds, it is rather difficult to properly solve advection without some kind of numerical diffusion or inaccuracies due to interpolation of properties between the tracers and the regular grid. The particle advection is used for the *tracer/marker in cell* method (either passive or active) and can be solved using different numerical methods, e.g., Euler integration or Runge-Kutta. The Eulerian form of the advection equation can also be solved in different ways. In the following, I would like to focus on *four* different methods to advect material. Different advection methods are used within the individual benchmarks and all can be tested in the [*Rigid Body Rotation*](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RigidBodyRotation) benchmark. 

### The Upwind Scheme

&emsp;The idea is that the flux into the local cell will only depend on the gradient of temperature in the direction upstream. The upwind scheme is similar to a forward in time and centered in space discretization, however, the central spacial derivatives are replaced by single-sided forward and backward finite differences and one needs to consider the advection velocity as well, to ensure that the discretization in space is always upstream. In 2-D the advection equation is then given as: 

$\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t}=-v_{x;i,j}\cases{\matrix{\frac{T_{i,j}^{n}-T_{i,j-1}^n}{\Delta x} \quad \text{if} \quad v_{x;i,j} > 0 \\\ \frac{T_{i,j+1}^{n}-T_{i,j}^n}{\Delta x} \quad \text{if} \quad v_{x;i,j}<0}} 
-v_{z;i,j}\cases{\matrix{\frac{T_{i,j}^{n}-T_{i-1,j}^n}{\Delta z} \quad \text{if} \quad v_{z;i,j} > 0 \\\ \frac{T_{i+1,j}^{n}-T_{i,j}^n}{\Delta z} \quad \text{if} \quad v_{z;i,j}<0}}$, &emsp;&emsp;&emsp;(5)

where *T* is the temperature, *v* the velocity, *n* is the current time step, $\Delta t$ the time step increment, and *i* and *j* are the indices in *z*- and *x*- direction, respectively. For more details see [*UpwindAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/UpwindAdvection2D.m).

&emsp;This is a stable and effective way, however, with a certain amount of numerical diffusion if the *courant criteria* is not fulfilled and only first order accurate in space. The courant criteria implies that the time step is smaller than the minimum grid spacing divided by the maximum velocity, that is, a property should not be advected over a distance larger than the grid spacing, or:

$\Delta t \le \frac{\Delta x}{max(|v|)}$.&emsp;&emsp;&emsp; (6)
   
### The Staggered Leap Frog (SLF) Scheme 

&emsp;This method considers a centered in time and centered in space discretization of the partial differentials, thus it has a higher order of accuracy in space (second order) and is suppose to not have any numerical diffusion. In 2-D the advection equation discretizes to:

$\frac{T_{i,j}^{n+1} - T_{i,j}^{n+1}}{2\Delta t}=-v_{x;i,j}\frac{T_{i,j+1}^{n} - T_{i,j-1}^{n}}{2\Delta x}-v_{z;i,j}\frac{T_{i+1,j}^{n} - T_{i-1,j}^{n}}{2\Delta z}$. &emsp;&emsp;&emsp; (7)

For more details see [*SLFAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/SLFAdvection2D.m).

### The semi-lagragian scheme 
&emsp;This method is related to the tracer-based advection by solving ordinary differential equations (*ODEs*), where it assumes that *imaginary tracers* are located at certain positions and land directly at the finite difference grid nodes after advection within one time step. Thus, one needs to calculate the *origin points* for each grid node back in time (e.g., one Euler time step) with a given velocity field (e.g., using an *iterative mid-point scheme*, i.e. one uses the velocity at a point half a time step backward in time) and then interpolate the property from the regular grid points to the determined *origin points*. This scheme assumes that no heat-sources were active during the advection. The method does not have any numerical diffusion but shows inaccuracies due to the interpolation method. For more details see [*SemiLagAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/SemiLagAdvection2D.m).
   
### Passive tracers

&emsp;Here, one assumes that the model domain is completely filled with so-called *tracers* or *markers*. These tracers are then advected by solving the *ODE* of a particle advection using a certain method (e.g., Euler or Runge Kutta) and they transport any property stored on them. However, care needs to be taken when interpolating those properties from the regular grid onto the tracers and back. This is even more complex if the property advected does have an effect on parameters controlling the governing equations (e.g., the viscosity in continuum euqation).<br>
&emsp;Here, I advect the tracers using Runge-Kutta fourth order; the tracers can transport the absolute temperature and the composition (so far only for two compositions with a constant viscosity and density). The property is then interpolated back to the regular grid points every time step. For more details see [*AdvectMarker2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/AdvectMarker2D.m) and [*TracerInterp.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/TracerInterp.m).

### Examples 
&emsp;For the [thermal convection](https://github.com/LukasFuchs/FDCSGm/tree/main/MixedHeatedSystems) code or in the [Blankenbach benchmark](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/Blanckenbach) I do prefer, so far, the semi-lagrangian method. However, I only advect the absolut temperature and not the increments, so far.

&emsp;In case of a temperature independent material (e.g., a [Rayleigh-Taylor Instability](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RTI) or the [falling block problem](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/FallingBlock)), advection with passive tracers should be used. The tracers advect the composition (so far, expressed by density and viscosity) and are advected using fourth order Runge-Kutta. The advection of composition partly works with a semi-lagrangian advection method, but problems do occur along the compositional boundaries.

&emsp;Besides for advection, the tracers are also very useful to define a certain model setup, like the in [viscous inclusion model](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/ViscousInclusion). The tracers are only used to define the geometry of the model, where the properties are then interpolated on the regular grid to solve for the steady state solution of the momentum equation.

&emsp;For more details and examples see the [benchmark directory](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark).

### Directory Content
For more details check description in the files.

[*AdvectMarker2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/AdvectMarker2D.m)<br>
&emsp;-> Script to advect passive tracers within a 2-D environment using Runge-Kutta 4th order.
   
[*Advection.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/Advection.m)<br>
&emsp;-> General function to chose a certain advection scheme, to be used in the main program.
   
[*SLFAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/SLFAdvection2D.m)<br>
&emsp;-> Staggered Leap Frog advection scheme for a 2-D environment. Advected is any kind of field *A*.

[*SemiLagAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/SemiLagAdvection2D.m)<br>
&emsp;-> Semi-Lagrangian advection scheme for a 2-D environment. Advected is any kind of field *A*.

[*TracerInterp.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/TracerInterp.m)<br>
&emsp;-> Script to interpolate certain properties from the grid to the tracers and back. 
   
[*UpwindAdvection2D.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/AdvectionProblem/UpwindAdvection2D.m)<br>
&emsp;-> Script to advect any kind of field *A* with a 2-D upwind scheme. 
