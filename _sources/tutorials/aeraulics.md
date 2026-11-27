# Aeraulics

This chapter deals with modeling airflow using flow networks, which are similar to the thermal networks presented in the section  on [Thermal networks for heat transfer](https://cghiaus.github.io/dm4bem_book/tutorials/T01_Heat_transfer.html). These networks are oriented weighted graphs that have conductances and pressure sources on branches and flow sources in nodes. This section addresses only steady-state flow networks. What is specific to aeraulic networks is that the constitutive equations are (generally) non-linear and relate volumetric flow to pressure drop, while the balance equations are based on mass conservation. Since mass flow rate is related to volumetric flow through density, and the density is that of the upstream air—which generally differs from the density of the downstream air—changing the direction of flows changes the values of the mass flow rates (it is not simply a change of sign).

## Summary of theory

### Sources
#### Pressure sources

Pressure sources [(Rock et al. 2017)](https://ia802804.us.archive.org/23/items/ashraehandbook2017hvacfundamentalssi/2017%2C%20ashrae%20fundamental%20handbook%20SI.pdf):
- stack pressure: $P_s = P_0 - \rho \, g \, z$;
- wind pressure: $P_w = C_p \frac{1}{2} \rho_o \, v^2$

where:
- $P_s$ stack pressure, Pa
- $P_0$ stack pressure at reference height, Pa
- $g$ gravitational acceleration, m/s²
- $\rho$ indoor or outdoor air density, kg/m³
- $z$ height above the reference level, m
- $P_w$ wind pressure difference relative to the outdoor pressure in undisturbed flow, Pa
- $C_p$ pressure coefficient, dimensionless
- $\rho_o$ outside air density, kg/m³
- $v$ wind speed, m/s

#### Flow sources
Flow sources inject a controlled volumetric airflow rate into the control volume, as is the case in [balanced mechanical ventilation systems](https://en.m.wikipedia.org/wiki/Ventilation_(architecture)).

### Constitutive equation

Each opening of a building can be modelled by a [constitutive equation](https://en.m.wikipedia.org/wiki/Constitutive_equation) or power law [(Rock et al. 2017)](https://ia802804.us.archive.org/23/items/ashraehandbook2017hvacfundamentalssi/2017%2C%20ashrae%20fundamental%20handbook%20SI.pdf):

$$ \dot V_{i,j} = K_{i,j} \, \Delta P^n_{i,j} = K_{i,j} \, (P_i - P_j)^n$$

where:
- $\dot V_{i,j}$ volumetric airflow through opening from control volume $i$ to control volume $j$, m³/s
- $n$ pressure exponent, dimensionless; typical value about 0.65
- $K_{i,j}$ flow coefficient, m³/(s·Paⁿ)

The values of the pressure exponent $n$ and flow coefficient $K_{i,j}$ can be obtained experimentally by using the [blower door](https://en.wikipedia.org/wiki/Blower_door) test.

### Mass conservation equation

Mass air flow rate from control volume $i$ to control volume $j$ is:

$$ \dot m_{i,j} = \rho_i \, \dot V_{i,j}$$

where:
- $\dot m_{i,j}$ mass flow rate of air from control volume $i$ to control volume $j$, kg/s 
- $\rho_i$ air density in control volume $i$, kg/m³
- $\dot V_{i,j}$ volumetric air flow rate from control volume $i$ to control volume $j$, m³/s

In steady-state, the [mass conservation](https://en.m.wikipedia.org/wiki/Conservation_of_mass) law states that the algebraic sum of mass flow rates entering in a control volume $j$ is zero:

$$\sum_i \dot m_{i,j} = 0$$

**Notes:**
1. The constitutive law relates volumetric airflow to pressure difference.
2. The balance equation is done for mass flow rate.
3. Generally, $\dot m_{i,j} \neq \dot m_{j,i}$, even if $\dot V_{i,j} = \dot V_{j,i}$. This implies that if the mass flow rates obtained by solving the aeraulic circuit are negative, than the direction of the flows need to be changed and a new solution to be searched for.
