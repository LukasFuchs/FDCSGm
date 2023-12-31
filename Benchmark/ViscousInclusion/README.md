# General information 

This directory contains three scripts: 

1. [*ViscousInclusion.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/ViscousInclusion/ViscousInclusion.m)<br>
   -> Script to calculate the instantaneous solution of the stokes equation for a viscous spherical or elliptical inclusion within a viscous medium under pure shear or simple shear deformation (define model and inclusion geometry, numerical resolution, viscosity contrast, and deformation field). In case the inclusion is a sphere, the solution is compared to the analytical solution (Schmidt and Podladchikov, 2003). The parameters are scaled by the constants as defined in */FDCSGm/ScaleParam/ScaleParameters.m*.

2. [*ViscousInclusion_Deta.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/ViscousInclusion/ViscousInclusion_Deta.m)<br>
   -> Script to calculate the instantaneous stokes solution for an elliptical inclusion with a certain orientation (*angle*) under pure shear or simple shear deformation for a range of different viscosity contrats. The script stores (*data* directory) the mean of the second invariant of the deviotoric stress, the strain-rate, and the "dissipative energy" $(\tau_{II} \cdot \varepsilon_{II})$ of the inclusion and the matrix for each orientaion. The data can be visualized with the script *DataVis.m* by defining the data directory and the deformation field.
   
3. [*ViscousInclusion_Rot.m*](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/ViscousInclusion/ViscousInclusion_Rot.m)<br>
   -> Script to calculate the instantaneous stokes solution for an elliptical inclusion with a certain viscosity contrast under pure shear or simple shear for a range of different orientations. The script plots the mean of the second invariant of the deviotoric stress, the strain-rate, and the "dissipative energy" $(\tau_{II} \cdot \varepsilon_{II})$ of the inclusion over the orientation.

-----------------------------------------------------------------

## Example of a spherical inclusion model

### **Constants**

Model height **H** [m]: 1000

Model Width **W** [m]: 1000

nx = 201

nz = 201

Background strain-rate $\dot{\varepsilon}_{bg}$ [1/s]: 10<sup>-15</sup>

Matrix viscosity $\eta_{mat}$ [Pa s]: 10<sup>23</sup> 

Inclusion viscosity $\eta_{inc}$ [Pa s]: 10<sup>21</sup>

Major half axis **a** [m]: 200

Minor half axis **b** [m]: 200

... Calculation of the deviation of the numerical solution from the analytical one ...

## Pure Shear

![Pure_Shear_Solution](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f441f7d4-80b3-417a-beb3-1cec1e82451a)<br>
**Figure 1.** Instantaneous solution of the stokes equation. **Top Left:** Logarithm of the viscosity. **Top Right:** Logartihm of the dissipative energy $\psi$. **Bottom Left:** Logarithm of the second invariant of the strain-rate tensor. **Bottom Right:** Logartihm of the second invariant of the deviatoric stress tensor. 

![Pure_Shear_Deviation](https://github.com/LukasFuchs/FDCSGm/assets/25866942/d1ec1d33-11e9-47d4-938b-5263a7fada02)<br>
**Figure 2.** Comparison with the analytical solution. **Top Row:** Analytical and numerical solution of the horizontal velocity and their deviation (in ?). **Middle Row:** Analytical and numerical solution of the vertical velocity and their deviation (in ?). **Bottom Row:** Analytical and numerical solution of the dynamic pressure and their deviation (in ?).

## Simple Shear

![Simple_Shear_Solution](https://github.com/LukasFuchs/FDCSGm/assets/25866942/10e04dc8-2c64-4f20-ba62-b244ff991b19)<br>
**Figure 3.** Same as in **Figure 1.** but for a simple shear background deformation field. See caption of **Figure 1.** for more detail.

![Simple_Shear_Deviation](https://github.com/LukasFuchs/FDCSGm/assets/25866942/06945e5a-2f73-4497-896c-cc24af13140d)<br>
**Figure 4.** Same as in **Figure 2.** but for a simple shear background deformation field. See caption of **Figure 2.** for more detail.

-----------------------------------------------------------------

## Example of an elliptical inclusion model

### **Constants**

Model height **H** [m]: 1000

Model Width **W** [m]: 1000

nx = 201

nz = 201

Background strain-rate $\dot{\varepsilon}_{bg}$ [1/s]: 10<sup>-15</sup>

Matrix viscosity $\eta_{mat}$ [Pa s]: 10<sup>23</sup> 

Inclusion viscosity $\eta_{inc}$ [Pa s]: 10<sup>18</sup> - 10<sup>28</sup>

Major half axis **a** [m]: 300

Minor half axis **b** [m]: 100

Orientation **α** [°]: 45

![Inclusion_Deta_PureShear_45](https://github.com/LukasFuchs/FDCSGm/assets/25866942/fc5038f8-388b-4048-b782-3e0a115930c5)<br>
**Figure 5.** Variation of each field with increasing inclusion viscosity for an orientaion of the elliptical inclusion of 45° and a pure shear backgraound deforamtion. For more details see caption of **Figure 1.**

### References
Schmid, D. W., & Podladchikov, Y. Y. (2003). Analytical solutions for deformable elliptical inclusions in general shear. Geophysical Journal International, 155(1), 269-288.

# Directory Content
