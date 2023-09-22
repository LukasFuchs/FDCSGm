# General Information

This directory contains two main scripts: 



Function to solve the two-dimensional advection problem with different
finite difference discrectization schemes. One can choose from:       
    'upwind'    - Upwind
    'slf'       - Staggered Leaped Frog
   'semi-lag'  - Semi-Lagrangian
   'tracers'   - Passive tracer method

 For the initial temperature anomaly one can choose:
   'block'     - Rectangular block
   'gaussian'  - Gaussian temperature distribution
   'circle'    - elliptical constant temperature anomaly

 For the constant background velocity field, one can choose from:
   'RigidBody' - Rotational field with constant rotation velocity
   'ShearCell' - Convection velocity field with shear deformation

