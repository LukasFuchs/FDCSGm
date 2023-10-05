# General Information
&emsp; This directory contains several function to set up a two dimensional model by defining the initial temperature or compositional field, initial velocity field, or allocating the memory of the used fields. There are also some functions to visualize the generated data and a function to control the break in the time loop. For more details on each function see the description in the files. In the following I would like to present the most common initial conditions implemented in the code. 

--------------------------------------------------------------------------------------------------------

# Initial Setup

The initial setup can be chosen for the temperature or the compositional field by simply defining the region for the different conditions. In some cases, tracers are used to define the geometry, where the property of the tracer are the interpolated back to the grid. For more details, see each individual function. 

## Temperature and Composition Condition

### Circular or Elliptical Anomaly

&emsp; A circular or elliptical anomaly can be defined by simply considering the ellipse/circle equation: 

$1 = \frac{(x-x_c)^2}{a^2} + \frac{z-z_c)^2}{b^2}$, &emsp; &emsp; &emsp; (1)

where x<sub>c</sub> and z<sub>c</sub> are the centre coordinates of the anomaly, and *a* and *b* the major and minor half axis, respectively (or in case of a circle, it's radius). If the equation is smaller than one for any grid point, or marker, then the corresponding point is located within the anomaly, and if it is larger than one, it is located outside. One can assign each point the corresponding conditions (e.g., certain temperature, viscosity, or density). 

### Rectangular Anomaly

&emsp; Similar to the elliptical anomaly, one simply needs to define the region, e.g., by defining the center, width, and height, and to check if a coordinate point is located inside or outside this region. 

### Gaussian Perturbation

&emsp; The initial Gaussian temperature anomaly is defined by:

$T=T_0 + A exp(-\frac{(x-x_c)^2 + (z-z_c)^2}{2\sigma^2 /\pi})$, &emsp; &emsp; &emsp; (2)

where x<sub>c</sub> and z<sub>c</sub> are the centre coordinates of the anomaly, *A* is the amplitude, *σ* its width, and T<sub>0</sub> is the background temperature.

### Linear Temperature Increase

&emsp; A linear temperature increase is simply defined with a certain depth gradient as: 

$T=T_0 + \frac{\partial T}{\partial z}z$, &emsp; &emsp; &emsp; (3)

where T<sub>0</sub> is the surface temperature, $\frac{\partial T}{\partial z}$ the temperature gradient, and *z* the depth. 

## Velocity Condition

### Rigid Body Rotation

The velocity for a rigid body rotation at every point in our two dimensional model domain can be defined using the definition of the tangential component of the angular velocity for any point within the rotating body: 

$v_T = \omega r$, &emsp; &emsp; &emsp; (4)

where $\omega$ is the rotational velocity and *r* the radius, i.e., the distance of the point to the rotational center. Considering equation (4) and some trigonometry, one can derive the horizontal and vertical velocity at each corresponding grid point as: 

$v_x = -\omega z_n = \omega \frac{(z-\Delta z/2)-H/2}{H}$, &emsp;&emsp;&emsp;(5) 
$v_z = \omega x_n = \omega \frac{(x-\Delta x/2)-L/2}{L}$, &emsp;&emsp;&emsp;(6)

where *x<sub>n</sub>* and *z<sub>n</sub>* are the normalized coordinates of each grid point, $\Delta z$ and $\Delta x$ are the vertical and horizontal grid space, *H* is the model height, and *L* the model width, assuming a staggered grid for the velocity field. 

### Convectional Shear Cell

### Pure Shear 

$\dot\varepsilon_{ij} = \left( \matrix{\dot\varepsilon_{xx} \quad 0\\\ 0 \quad \dot\varepsilon_{zz}} \right)$

$\dot\varepsilon_{xx}=\frac{\partial v_x}{\partial x}$

$\dot\varepsilon_{zz}=\frac{\partial v_z}{\partial z}$

$v_x=x \cdot \dot\varepsilon_{bg}$

$v_z=z \cdot \dot\varepsilon_{bg}$

$\dot\varepsilon_{bg} = \dot\varepsilon_{xx} = \dot\varepsilon_{zz}$

### Simple Shear

$\dot\varepsilon_{ij} = \left( \matrix{0 \quad \dot\varepsilon_{xz}\\\ \dot\varepsilon_{zy} \quad 0} \right)$

$\dot\varepsilon_{xz}=\dot\varepsilon_{zx}=\frac{1}{2} \frac{\partial v_x}{\partial z}$

$v_x = 2 \dot\varepsilon_{bg} \cdot z$

$\dot\varepsilon_{bg}=\dot\varepsilon_{xz}$

# Directory Content
[CheckBreakCriteria.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/CheckBreakCriteria.m)<br>
&emsp;-> Function to controll the break in the time loop, e.g., once steady state or the maximum time is reached. 

[PlotData.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/PlotData.m)<br>
&emsp; -> General function to plot temperature, velocity, and viscosity field for a thermal convection. 

[PlotProfile.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/PlotProfile.m)<br>
&emsp; -> Function to plot a profile at a certain position. 

[PlotTimeSerieses.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/PlotTimeSerieses.m)<br>
&emsp; -> Function to plot the time series of the Nusselt number and the V<sub>RMS</sub>, the mean temperature profile, and the heat fluxes at the corners. 

[SetUpFields.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/SetUpFields.m)<br>
&emsp; -> Function to initialize the memory of the used fields.

[SetUpInitialConditions.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/SetUpInitialConditions.m)<br>
&emsp; -> Function to set up initial temperature and compositianal conditions. 

[plotfield.m](https://github.com/LukasFuchs/FDCSGm/blob/main/SetUp/plotfield.m)<br>
&emsp; -> General function to plot any two dimensional field generated by the code.
