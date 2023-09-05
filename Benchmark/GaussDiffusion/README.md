

explicit, nx = nz = 41 <br>
![Evolution](https://github.com/LukasFuchs/FDCSGm/assets/25866942/3cff6778-028d-48ce-b63d-7afff14b8c2c)

Gaussian initial temperature anomaly; <br>

$T=T_0 + A exp(-\frac{(x-0.5L)^2 + (z-0.5H)^2}{2\sigma^2 /\pi})$

A = 200 [K], <br>
σ = 10 [% of L], <br>
T<sub>0</sub> = 1000 [K]. <br>

**Analytical Solution**<br>
$T_{ana} = T_0 + \frac{A}{(1 + 2\pi t \kappa / \sigma^2)} exp(-\frac{(x-0.5L)^2 + (z-0.5H)^2}{2\sigma^2 / \pi + 4t\kappa})$


## Resolution Test<br>
explicit, implicit, CNV, ADI<br>
after 10 Million years<br>

![Comparison](https://github.com/LukasFuchs/FDCSGm/assets/25866942/b4bfe7a0-96e1-43b5-8656-02269bf06e67)
