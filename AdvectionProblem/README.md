# General Information 

This directory contains seven files: 

1. skdfj

2. sdf

3. sdf

4. sdf

5. sdf

6. sdf

7. sdf

Advection can be described by: 

$$.&emsp;&emsp;&emsp;(1)

Euler form

$$. &emsp;&emsp;&emsp;(2)

Lagrangian form: 

$$. &emsp;&emsp;&emsp; (3)

## Discretization Schemes

As simple as this problem sounds, it is rather difficult to preserve the initial shape, mainly due to numerical diffusion or due to inaccuracies of an interpolation. Here, I use *four* different advection schemes to advect the temperature: 

1. *The upwind scheme*<br>
    -> A stable and effective way, however, with a certain amount of numerical diffusion if the *courant criteria* is not fulfilled. The courant criteria implies that the time step is smaller than the minimum grid spacing divided by the maximum velocity, that is, a property should not be advected over a distance larger than the grid spacing, or:<br>
$\Delta t \le \frac{\Delta x}{max(|v|)}$<br>
Here, I use a courant criteria of one. The upwind scheme is similar to a forward in time and centered in space discretization, however, one needs to consider the advection velocity as well, to ensure that the discretization in space is always upstream.<br>
   
2. *The Staggered Leap Frog (SLF) scheme*<br>
    -> This method considers a centered in time and centered in space discretization of the partial differentials, thus it has a higher order of error and is suppose to not have any numerical diffusion. As promissing as this scheme sounds it is not properly working here yet (I believe, could also be some boundary condition effects)!<br>
   
3. *The semi-lagragian scheme*<br>
    -> This method assumes that an *imaginary tracer* is located at a certain position and lands directly at a finite difference grid node after advection within one time step. Thus, one needs to calculate the *origin point* from each grid node back in time with a given velocity field (using a central point iteration method) and then interpolate the property from the regular grid points to the determined *origin points*. The method does not have any numerical diffusion but shows inaccuracies due to the interpolation method.<br>
   
5. *Passive tracers*<br>
    -> Here, one assumes that the model domain is completely filled with so-called tracers or markers. These tracers are then advected by a certain method (e.g., Euler or Runge Kutta) and they transport any property stored on them. However, care needs to be taken when interpolating those properties from the regular grid onto the tracers and back. This is even more complex if the property advected does have an effect on parameters controlling the governing equations (e.g., the viscosity in continuum euqation). Here, I advect the tracers using Runge-Kutta fourth order; the tracers do transport the absolute temperature, which is interpolated only to the regular grid points every time step and not back to the tracers again (since it is not suppose to change here).<br>

### Upwind

### Staggered Leap Frog

### Semi-lagrange

### Tracers
