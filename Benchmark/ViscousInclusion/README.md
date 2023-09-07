This directory contains three scripts: 

1. ViscousInclusion.m<br>
   -> Script to calculate the instantaneous solution of the stokes equation for a viscous spherical or elliptical includion within a viscous medium under pure shear or simple shear deformation (define model and inclusion geometry, numerical resolution, viscosity contrast, and deformation field). In case the includion is a sphere, the solution is compared to the analytical solution (Schmidt and Podladchikov, 2003). The parameters are scaled by the constants as defined in */FDCSGm/ScaleParam/ScaleParameters.m*.

2. ViscousInclusion_Deta.m<br>
   -> Script to calculate the instantaneous stokes solution for an elliptical inclusion with a certain orientation (*angle*) under pure shear or simple shear for a range of different viscosity contrats. The final plot shows the mean of the second invariant of the deviotoric stress, the strain-rate, and the "dissipative energy" $(\tau_{II} \cdot \varepsilon_{II})$ of the inclusion over the viscosity contrast for each orientaion.
   
3. ViscousInclusion_Rot.m<br>
   -> Script to calculate the instantaneous stokes solution for an elliptical inclusion with a certain viscosity contrast under pure shear or simple shear for a range of different orientations. The final plot shows the mean of the second invariant of the deviotoric stress, the strain-rate, and the "dissipative energy" $(\tau_{II} \cdot \varepsilon_{II})$ of the inclusion over the orientation.

Background strain-rate <br>
$\dot{\varepsilon}_{bg} = 10^{-15} 1/s$

nx = nz = 201 <br>

$\eta_{inc} = 10^{21}$ [Pa s] <br>
$\eta_{mat} = 10^{23}$ [Pa s] <br>

## Pure Shear <br>
![ViscousInclusionFields_PureShear](https://github.com/LukasFuchs/FDCSGm/assets/25866942/cd39d4ee-4a5a-4ce8-93dd-0961baba62c3)

![ViscousInclusionCompAna_PureShear](https://github.com/LukasFuchs/FDCSGm/assets/25866942/e13fdff3-522c-4545-9d62-fea1bf1888db)

## Simple Shear <br>
![ViscousInclusionFields_SimpleShear](https://github.com/LukasFuchs/FDCSGm/assets/25866942/541fff64-28bf-42de-b7e9-fef9efc10d0d)

![ViscousInclusionCompAna_SimpleShear](https://github.com/LukasFuchs/FDCSGm/assets/25866942/e087915d-6771-4819-9bc7-0bc915230e90)
