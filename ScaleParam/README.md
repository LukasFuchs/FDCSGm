# General Information

&emsp;To better compare different kinds of thermal convection, one can scale the governing equations. The equation can be scaled by the following reference parameters: 

$h_{sc} = h$,&emsp;&emsp;&emsp;(1)

$\eta_{sc} = \eta_0$,&emsp;&emsp;&emsp;(2)

$t_{sc} = \frac{h^2}{\kappa}$,&emsp;&emsp;&emsp;(3)

$v_{sc} = \frac{\kappa}{h}$,&emsp;&emsp;&emsp;(4)

$\tau_{sc} = \frac{\eta_0 \kappa}{h^2}$,&emsp;&emsp;&emsp;(5)

$T_{sc} = \Delta T$,&emsp;&emsp;&emsp;(6)

$H_{sc} = \frac{c_p \Delta T \kappa}{h^2}$,&emsp;&emsp;&emsp;(7)

where the subscript *sc* stands for the scaling parameters, and *h*, *η<sub>0</sub>*, *t*, *κ*, *v*, *τ*, *T*, *H*, *c<sub>p</sub>*, are the height, the reference viscosity, the time, the thermal diffusivity, the velocity, the stress, the temperature, the heat generation source, and the specific heat capacit, respectively. 

&emsp;For the given scaling parameters, the non-dimensional governing equations are given as (assuming a constant viscosity; scaling for a variable viscosity is applicable in the same way):

$\frac{\partial v_x}{\partial x} + \frac{\partial v_z}{\partial z} = 0$,&emsp;&emsp;&emsp;(8)

$\frac{\partial T}{\partial t} + v_x \frac{\partial T}{\partial x} + v_z \frac{\partial T}{\partial z} = (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2} + H)$,&emsp;&emsp;&emsp;(9)

$-\frac{\partial P}{\partial x} + \eta \frac{\partial^2 v_x}{\partial x^2} + \eta \frac{\partial^2 v_x}{\partial z^2} = 0$,&emsp;&emsp;&emsp;(10)

$-\frac{\partial P}{\partial z} + \eta \frac{\partial^2 v_z}{\partial z^2} + \eta \frac{\partial^2 v_z}{\partial x^2} - RaT = 0$,&emsp;&emsp;&emsp;(11)

where *Ra* is the so-called thermal *Rayleigh number* and *P* the *dynamic pressure*. In case of a basally heated thermal convection, the convective vigor is defined by the Rayleigh number, which describes a relationship between heat transported by buoyancy and conduction, and the effect of the layers thickness and bulk viscosity.
