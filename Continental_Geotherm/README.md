# General Information

&emsp;This is an example model for a two dimensional setup of a continental lithosphere with or without a (upper and lower) crust including radioactive heat sources. The script does calculate the steady-state solution of a 2-D temperature field for assuming variable thermal parameters (varying for each layer) and stores a temperature-depth profile in an additional *.txt file. 

One needs to define:
- the geometry of the setting, that is, the thickness of each layer,
- the thermal parameters for each layer, that are,
  - $\rho$, the density [kg/m<sup>3</sup>],
  - *k*, the conductivity [W/m/K], and
  - *H*, the heat production rate per mass [W/kg], and
- the thermal boundary conditions, that are,
  - Dirichlet, or
  - Neumann.
 
The script was used, e.g., to calculate the initial temperature field of a certain geodynamic lithosphere model with an upper and lower crust and is listed here simply as an example of how to use the code. For more details, see [ContinentalGeotherm.m](https://github.com/LukasFuchs/FDCSGm/blob/main/Continental_Geotherm/ContinentalGeotherm.m).

![Example](https://github.com/LukasFuchs/FDCSGm/assets/25866942/8aa3d185-e993-4ac6-baee-0b185381ba62)<br>
**Figure 1.** Example of a 2-D continental lithosphere setup with upper and lower crust. 
