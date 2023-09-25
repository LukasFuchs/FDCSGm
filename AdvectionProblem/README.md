# General Information 

This directory contains several files to calculate the advection of a certain property within a two-dimensional domain by the e.g., upwind, semi-lagrangian, or passive tracer method. For the sake of simplicity, I do focus on the advection of the absolute temperature in the following, however, in theory, one can chose any kind of property for advection (with some conditions to keep in mind). 

In general, advection describes the transport of a property, here the temperature, from one point to another, where one can assume different reference frames for the given point of interest. If we assume a not moving reference frame (that is an Eulerian grid), the change in temperature at a certain point can be described by (*eulerian advective transport equation*): 

$\frac{\partial T}{\partial t} = - \overrightarrow{v} \cdot \nabla T$,&emsp;&emsp;&emsp;(1)

or in a Lagrangian reference frame (i.e., along a moving point) as: 

$\frac{DT}{Dt}$, where both are related by: $\frac{DT}{Dt} = \frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T$.&emsp;&emsp;&emsp;(2)

For a Lagrangian reference point advection is given by a simple *ODE* particle advection scheme, where changes in its coordinates are related with the material velocities as: 

$\frac{Dx}{Dt} = v_i$,&emsp;&emsp;&emsp; (3)

where i is the coordinate index and x<sub>i</sub> is a spatial coordinate. 

--------------------------

## Discretization Schemes

As simple as the advection problem sounds, it is rather difficult to properly solve advection without some kind of numerical diffusion or inaccuracies due to interpolation of properties between the tracers and the regular grid. The particle advection is used for the *tracer/marker in cell* method (either passive or active) and can be solved using different numerical methods, e.g. Euler integration or Runge-Kutta. The Eulerian form of the advection equation can also be solved in different ways. In the following, I would like to focus on *four* different methods to advect material. Different advection methods are used within the individual benchmarks and all can be tested in the *Rigid Body Rotation* benchmark. 

### The upwind scheme

The idea is that the flux into the local cell will only depend on the gradient of temperature in the direction upstream. The upwind scheme is similar to a forward in time and centered in space discretization, however, the central spacial derivative is replaced by single-sided forward and backward finite differences and one needs to consider the advection velocity as well, to ensure that the discretization in space is always upstream. In 2-D the advection equation is then given as: 

![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/09474b14-e6c8-4cc4-a5c8-d97eb0e84406),

where *i* and *j* are the indices in *z*- and *x*- direction, respectively. 

This is a stable and effective way, however, with a certain amount of numerical diffusion if the *courant criteria* is not fulfilled and only first order accurate in space. The courant criteria implies that the time step is smaller than the minimum grid spacing divided by the maximum velocity, that is, a property should not be advected over a distance larger than the grid spacing, or:<br>
$\Delta t \le \frac{\Delta x}{max(|v|)}$<br>
   
### The Staggered Leap Frog (SLF) scheme 

This method considers a centered in time and centered in space discretization of the partial differentials, thus it has a higher order of accuracy and is suppose to not have any numerical diffusion. In 2-D the advection equation discretizes to:

![image](https://github.com/LukasFuchs/FDCSGm/assets/25866942/6b13c8ad-0ec7-4248-a114-90b1b87d3eaf).

   
### The semi-lagragian scheme 
This method assumes that an *imaginary tracer* is located at a certain position and lands directly at a finite difference grid node after advection within one time step. Thus, one needs to calculate the *origin point* from each grid node back in time with a given velocity field (using a central point iteration method) and then interpolate the property from the regular grid points to the determined *origin points*. The method does not have any numerical diffusion but shows inaccuracies due to the interpolation method.<br>
   
### Passive tracers

Here, one assumes that the model domain is completely filled with so-called tracers or markers. These tracers are then advected by a certain method (e.g., Euler or Runge Kutta) and they transport any property stored on them. However, care needs to be taken when interpolating those properties from the regular grid onto the tracers and back. This is even more complex if the property advected does have an effect on parameters controlling the governing equations (e.g., the viscosity in continuum euqation). Here, I advect the tracers using Runge-Kutta fourth order; the tracers do transport the absolute temperature, which is interpolated only to the regular grid points every time step and not back to the tracers again (since it is not suppose to change here).<br>

### Examples 
   For the thermal convection code I do prefer, so far, the semi-lagrangian method. However, I do only advect the total temperature and not the increments, so far. In case of a temperature independent material (e.g., a Rayleigh-Taylor Instability), advection with passive tracers, wich advect different properties (so far, density and viscosity), should be used (this partly works with a semi-lagrangian advection method, too). The tracers are advected using fourth order Runge-Kutta. For more details and examples see the [benchmark directory](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark).

