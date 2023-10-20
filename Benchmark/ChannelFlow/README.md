# General Information

&emsp; Assuming the horizontal pressure gradient is constant and flow within a channel is only driven by the pressure and/or by a constant horizonal velocity at the surface (or at the bottom, or both), the stokes equation describes the horizontal flow velocity within the channel and simplifies to: 

$\frac{\partial P}{\partial x} = \frac{\partial \tau_{x,z}}{\partial z}$, &emsp;&emsp;&emsp; (1)

where *P* is the pressure and $\tau_{x,z}$ is the deviatoric shear stress, which is defined as: 

$\tau_{x,z} = 2 \eta \dot{\varepsilon}_{x,z} = \eta \frac{\partial v_x}{\partial z}$, &emsp;&emsp;&emsp; (2)

where $\eta$ is the dynamic viscosity and $\dot{\varepsilon}_{x,z}$ is the deviatoric shear strain-rate, which is defined as: 

$\dot{\varepsilon}_{x,z} = \frac{1}{2}(\frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x})$.&emsp;&emsp;&emsp; (3)

For the given setup I can assume that the vertical velocity is zero and thus equation (3) simplifies to the last expression of equation (2).

&emsp; This directory contains a script to calculate the horizontal velocity for a two-dimensional Couette(-Poiseuille) channel flow with constant and logarithmically, with depth varying viscosity and to compare the numerical solution with its analytical solution. The depth-dependent viscosity is defined as: 

$\eta = \eta_0 exp(log(m) \frac{H-z}{H})$,&emsp;&emsp;&emsp; (4)

where *m* is the viscosity ratio of $\frac{\eta_1}{\eta_0}$, $\eta_0$ and $\eta_1$ are the bottom and surface viscosities, respectively, *H* is the model height, and *z* the depth. 

&emsp;Considering the definition of the viscosity as given in equation (4), one can derive an analytical solution of the horizontal velocity from the 1-D stokes equation in *x*-direction by twice integrating equation (1). The analytical solution with depth depends on the viscosity ratio *m*, the horizontal pressure gradient $\frac{\partial P}{\partial x}$, and the shear velocity at the surface $v_{x,0}$. For an upward pointing coordinate system (*z* positive) the analytical solution is given as: 

$v_{x,ana} =-\frac{1}{2 \eta_0} \frac{\partial P}{x} (Hz - z^2) + v_{x,0}\frac{z}{H}$,&emsp;&emsp; if $m = 1$, and &emsp;&emsp;&emsp; (5)

$v_{x,ana} = -\frac{\partial P}{\partial x} \frac{H}{\eta_0 log(m)} (\frac{m^{-\frac{z}{H}}}{m-1}(z(m-1)+H) - \frac{H}{m-1})-m^{-\frac{z}{H}} m \frac{v_{x,0}}{m-1} + \frac{v_{x,0}m}{m-1}$, &emsp;&emsp; if $m \neq 0$.&emsp;&emsp;&emsp; (6)

&emsp;The numerical solution is calculated using fixed boundary velocities, which are defined by the analytical solution of the horizontal velocity as defined in equations (5) and (6) and I simply flip the analytical solution so that it fits to the downward point coordinate system I use in the code. 

---------------------------------------------------------------

## Examples 

![Couette_ISO](https://github.com/LukasFuchs/FDCSGm/assets/25866942/b9bc9e28-c145-44b9-b0e9-bbdbd627cd4e)<br>
**Figure 1.** Solution of an isoviscous Couette flow. TOP: Viscosity field overlain by the velocity vectors. BOTTOM LEFT: Mean horizontal velocity of the numerical solution (solid black line) overlain by its analytical solution (dashed yellow line). BOTTOM RIGHT: Root mean square devation of the numerical solution from its analytical solution. 

![Couette_EXP](https://github.com/LukasFuchs/FDCSGm/assets/25866942/f228cfef-680f-4b5b-8e35-b52134924c2d)<br>
**Figure 2.** Solution of a Couette flow with logarithmically varying viscosity $(m = 10^6)$. For more details see caption of **Figure 1**.

![CouettePoiseuille_ISO](https://github.com/LukasFuchs/FDCSGm/assets/25866942/7c48f1ba-d593-4c3e-a121-a3c3a3bca115)<br>
**Figure 3.** Solution of an isoviscous Couette-Poiseuille flow $(\frac{\partial P}{\partial x} = -10)$. For more details see caption of **Figure 1**.

![CouettePoiseuille_EXP](https://github.com/LukasFuchs/FDCSGm/assets/25866942/6f10316d-22a2-4205-9ac0-307e2f407bf3)<br>
**Figure 4.** Solution of a Couette-Poiseuille flow $(\frac{\partial P}{\partial x} = -200)$ with logarithmically varying viscosity $(m = 10^6)$. For more details see caption of **Figure 1**.


# Directory Content

[NSE_test_deriv.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/ChannelFlow/NSE_test_deriv.m)<br>
&emsp; -> Script to calculate the 1-D stokes solution for a Couette-(Poiseuille) channel flow with variable viscosity using different finite discretization schemes (staggered, centered 2nd order, centered 4th order, chain rule) and compare it to its analytical solution.

[ChannelFlow.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Benchmark/ChannelFlow/ChannelFlow.m)<br>
&emsp; -> Script to calculate the 2-D solution of a horizontal Couette(-Poiseuille) channel flow for constant or depth-dependent viscosity condition and compares the numerical solution with its analytical. 
