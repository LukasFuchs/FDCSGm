# General Information

&emsp; This directory contains several functions to solve the 2-D stokes, for both, constant and variable viscosities, and continuum equation for the velocities and the pressure assuming a staggered finite difference grid. There are also some functions to interpolate the velocities from the staggered grid to the regular grid, as well as, functions to calculate the deviatoric strain-rate and stress tensors. 

------------------------------------------------------------------------------------------

## Stokes Equation 
&emsp;On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one needs consider, i.e., *inertia*, *surface*, and *volumetric* forces. A common equation to describes such motion is given by: 

$\rho \frac{D \overrightarrow{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \boldsymbol{\rho}$,&emsp;&emsp;&emsp;(1)

where *ρ* is the density [kg/m<sup>3</sup>], $\overrightarrow{v}$ is the velocity vector [m/s], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [Pa], $\boldsymbol{g}$ is the gravitational acceleration [m/s<sup>2</sup>], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial}{\partial t} + \overrightarrow{v} \cdot \nabla$. The *Cauchy stress tensor* is given by: 

$\boldsymbol{\sigma} = -\nabla P + \nabla \cdot \boldsymbol{\tau}$,&emsp;&emsp;&emsp;(2)

where *P* is the total pressure (*P = P<sub>dynamic</sub> + P<sub>hydrostatic</sub>*) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. In Eulerian form, equation (1) is given by (Navier-Stokes equation):

$\rho (\frac{\partial v_{i}}{\partial t} + v_{j}\frac{v_{i}}{\partial x_{j}}) = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(3)

where summation over repeated indices is implied. 

&emsp;To solve equation (3), one needs to define a rheology which, for a purely viscous medium, can be given by a constitutive relationship between stress and strain rate in the form of, e.g.:

$\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij}$,&emsp;&emsp;&emsp;(4)

where $\eta$ is the dynamic viscosity [Pa s] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* [1/s] and given by: 

$\dot{\varepsilon}_{ij} = \frac{1}{2} (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i})$.&emsp;&emsp;&emsp;(5)

&emsp;Assuming that the inertia forces are negligible in comparison to the gravitational forces, one can further simplify equation (3) to:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\tau_{ij}}{\partial x_j} + \rho g_{i}$,&emsp;&emsp;&emsp;(6)

or in the form of the unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*:

$0 = -\frac{\partial P}{\partial x_{i}} + \frac{\partial}{\partial x_j} \eta (\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i}) + \rho g_{i}$,&emsp;&emsp;&emsp;(7)

Assuming constant viscosity equation (7) simplifies further to (Stokes equation): 

$0 = -\frac{\partial P}{\partial x_{i}} + \eta (\frac{\partial^2 v_i}{\partial x_j^2} + \frac{\partial^2 v_j}{\partial x_i^2}) + \rho g_{i}$.&emsp;&emsp;&emsp;(8)

&emsp;Equation (8) provides us two equations for our three unknowns. Thus, one needs to also consider the mass conservation equation (i.e., we do work with a continuum), where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$\frac{\partial v_i}{\partial x_i} = 0$.&emsp;&emsp;&emsp;(9)

Equations (8) and (9) enable us to solve for the three unknowns *v<sub>x</sub>*, *v<sub>z</sub>*, and *P*. 

&emsp;The conservation equations of *momentum* and *mass* are solved properly in two dimensions (*x* and *z*, *z* pointing downwards in negative direction) using a staggered finite difference grid, based on the description by *Gerya* (2009), where the viscosity and the density are defined on the regular grid points, the horizontal and vertical velocity in between the regular grid points, and the pressure within a finite difference cell. A staggered grid enables the conservation of the stress between adjacent grid points and one can solve equations (8) and (9) for the unknows. For more details on the discretization and the solving equations see *Gerya* (2009). 

### Equation of State
&emsp;The buoyance term on the right-hand side of equation (7), that is the density term which is temperature dependent (and pressure, but I do neglect this effect here so far), can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$\rho = \rho_0 (1-\alpha T)$,&emsp;&emsp;&emsp;(10)

where *ρ<sub>0</sub>* is the reference density and *α* the thermal expansion coefficient [1/K]. 

## Internal Structure

&emsp; I initially wrote this code during my [PhD](https://onedrive.live.com/?authkey=%21AHBJz3oP27xd42Y&cid=775D7EE18409BBFB&id=775D7EE18409BBFB%21674908&parId=775D7EE18409BBFB%21623998&o=OneUp) and, unfortunately, the documentation was quite brief. However, the main scripts to solve for the velocities and the pressure are the functions [solveSECE](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE.m)[Sc.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECESc.m) and [solveSECE_const_Eta](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE_const_Eta.m)[Sc.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE_const_EtaSc.m) in a dimensional and non-dimensional (files ending with *Sc.m) form. Each function follows the same procedure, which is given as 

- setting up the $rhs$ vector, using a vertical running scheme and ordering the unknowns per grid point as *P*, *v<sub>x</sub>*, and *v<sub>z</sub>*,
- setting up the coefficient matrix,
- solving the sysmtem of equations, and 
- updating the velocity and pressure field.

&emsp;Setting up the coefficient matrix, I also include the boundary nodes, and the order of the unknows for the inner grid points is *P*, *v<sub>x</sub>*, and *v<sub>z</sub>*, too. The function is structured in such a way, that I first collect all coefficients for the *x*-stokes equation, second for the *z*-stokes equation, and finally for the continuum euqation. The coefficients are stored in an array with 3 * (3 * 4 + 2 ) columns, which is the total number of used grid points for the finite difference stencils or the number of non-zero diagonals, and 3 * nx * nz rows, i.e., the total number of equations. I finally use this array to setup the coefficient matrix in MATLAB with the *spdiags* command, i.e., I do use a sparse matrix format. The system of equation is then solved with the right-array division in MATLAB. 

&emsp;Since this structure of the code is working and has been tested with respect to analytical solutions for different [channel flows](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/ChannelFlow) and additional [benchmarks](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark), I never intended to optimize *my* structure for the stokes solver and, moreover, I do not let the students write this solver. Hoever, I am open for suggestions (highly appreciated and thanks in advance!). 

### Value Interpolation

&emsp;The velocity on the staggered grid is interpolated (harmonic average) to the regular grid points for visualization purposes, as well as, to advect the material. This will be optimized in the near future. 

# References 

*Gerya, T. (2019). Introduction to numerical geodynamic modelling. Cambridge University Press.*

# Directory Content
[GetStrainRate.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/GetStrainRate.m)<br>
&emsp;-> Function to solve for the strain-rate on the **regular** grid

[GetStrainRateStress.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/GetStrainRateStress.m)<br>
&emsp;-> Function to solve for the strain-rate and stress on the **staggered** grid. (Not correctly implemented yet!)

[InterpStaggered.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/InterpStaggered.m)<br>
&emsp; -> Function to interpolate the viscosities from the staggered onto the regular grid.

[solveSECE.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE.m)<br>
&emsp; -> Function to solve the stokes and continuum equation for a variable viscosity and in dimensional form. 

[solveSECESc.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECESc.m)<br>
&emsp; -> Function to solve the stokes and continuum equation for a variable viscosity and in a non-dimensional form. 

[solveSECE_const_Eta.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE_const_Eta.m)<br>
&emsp; -> Function to solve the stokes and continuum equation for a constant viscosity and in a dimensional form. 

[solveSECE_const_EtaSc.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/solveSECE_const_EtaSc.m)<br>
&emsp; -> Function to solve the stokes and continuum equation for a constant viscosity and in a non-dimensional form. 

[velgrad.m](https://github.com/LukasFuchs/FDCSGm/blob/main/StokesProblem/velgrad.m)<br>
&emsp; -> Function to calculate a horizontal or vertical velocity gradient.

