# Routines to Solve the Diffusion Problem

   This directory contains all rountines to solve the diffusive part of the *temperature conservation equation* (1- and 2-D) using different numerical discretization methods. For more details on the discretization methods see the Numerical Methods section in *FDCSGm/README.md*. The routines are avaible for a dimensional or non-dimensional (files ending with *Sc.m) form of the equation (so far the 1-D routines are only availabe for a dimensional version!). 

## 1-D Geotherms
   The 1-D temperature profile is calculated by solving diffusive parte of the 1-D temperature conservation equation (so far only with a radiogenic heat source) for variable thermal parameters with a proper conserving finite difference scheme. That is, the heat flow is calculated on the centered and the remaining parameters on the regular grid points, respectively. The discretization scheme for variable thermal parameters is picked to solve for a temperature profile of a continental lithosphere with upper, lower crust, and mantle. 
The 1-D heat equation is given by: 

$\rho c_{p} \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}(k \frac{partial T}{\partial z}) + \rho H$

### Oceanic Geotherms

### Continental Geotherms

## Steady State Solution
