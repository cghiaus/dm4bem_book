# Airflow networks

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


## Physical principles: balance equations and constitutive laws

An airflow network represents the movement of air through openings, cracks, ducts, or other flow paths connecting different pressure nodes. The network is governed by two fundamental types of equations:

1.  **Balance equations**, expressing conservation of mass at the
    pressure nodes.
2.  **Constitutive laws**, describing the relationship between the
    airflow rate through a branch and the pressure difference across it.

These equations form a nonlinear system because the airflow through an opening generally depends nonlinearly on the pressure difference.

### Mass conservation at pressure nodes

Consider a network composed of (N) internal pressure nodes connected by (B) branches. Each branch connects two nodes or an internal node to a reference-pressure node.

At each internal node, conservation of mass requires that the algebraic sum of the mass flow rates through the connected branches, together with any imposed mass-flow source, be zero.

Using a consistent orientation for the branches, the mass balance
equations can be written as

$$\sum \dot{m}_i + f_i = 0$$

where:

-   $\sum$ is the summation operator;
-   $\dot{m}_i$ are the branch mass flow rates (kg/s);
-   $f_i$ are the imposed mass-flow sources
    at internal nodes (kg/s).

The balance equations express the fact that air cannot accumulate at a pressure node in a steady-state airflow network.

### Nonlinear constitutive law: volumetric flow and pressure difference

The airflow through an opening or leakage path is commonly described by
a power-law relationship between the volumetric flow rate and the
pressure difference across the opening:

$$\dot{V}_{ij} = K_{ij} \operatorname{sign} (\Delta p_{ij}) |\Delta p_{ij}|^{n_{ij}}$$

where:

-   $\dot{V}_{ij}$ is the volumetric flow rate, $m^3/h$;
-   $K_{ij}$ is the flow coefficient, $m^3 /(h·Pa^{n_{ij}})$;
-   $\Delta p_{ij}$ is the pressure difference across the branch, $Pa$;
-   $n_{ij}$ is the flow exponent, with $0 < n_{ij} \leq 1$;
-   $\operatorname{sign} (\Delta p_{ij})$ indicates the direction of the flow.

The flow exponent depends on the characteristics of the opening and the
flow regime. Typical values include ([Orme & Leksmono, 2002](https://www.aivc.org/sites/default/files/members_area/medias/pdf/Guides/Ventilation%20Modelling%20Data%20Guide_web.pdf), Appendix E):

-   $n = 0.5$: turbulent flow, commonly used for openings;
-   $n$ approaching 1: predominantly laminar flow.

The $\operatorname{sign}$ function is essential because the pressure difference can
change sign, causing the flow direction to reverse.

### Relationship between mass flow and volumetric flow

Mass flow rate and volumetric flow rate are related through air density:

$$\dot{m}_{ij} = \rho_{ij} \dot{V}_{ij}$$

where:

- $\dot{m}_{ij}$ is the  mass flow, kg/s; 
- $\rho_{ij}$ is the air density, kg/m³;
- $\dot{V}_{ij}$ the volumetric flow rate, m³/s.

Substituting the power-law constitutive relationship gives

$$\dot{m}_{ij} = \rho_{ij}K_{ij} \operatorname{sign} (\Delta{p}) |\Delta p_{ij}| ^{n_{ij}}$$

Note that air density depends on the location from which the air originates.

### Upstream air density

In an airflow network, the density used to calculate mass flow must correspond to the upstream node of the branch, that is, the node from which the air enters the branch.

For a branch oriented from node $i$ to node $j$:

-   if $\Delta p_{ij}$ produces positive flow, air travels
    from $i$ to $j$, and the density is $\rho_i$;
-   if $\Delta p_{ij}$ produces negative flow, air travels
    from $j$ to $i$, and the density is $\rho_j$.

Therefore, the mass flow relationship becomes

$$\dot{m}_{ij} = K_{ij}
\rho_{up,ij} \operatorname{sign} (\Delta p_{ij})
|\Delta p_{ij}|^{n_{ij}}$$

where $\rho_{up,ij}$ is the density at the
upstream end of the branch.

This distinction is particularly important when a branch connects  environments with different air densities.

In the Python implementation, the upstream density is selected according to the sign of the branch pressure difference. The code stores two possible densities for each branch, `rho_forward` and `rho_backward`, corresponding to the two possible flow directions, and selects the appropriate one during each iteration.

## Pressure sources and flow sources

An airflow is driven by imposed pressure differences and imposed mass flows.

### Pressure sources

Pressure sources contribute to the pressure difference across a branch
independently of the unknown pressure differences between its internal
nodes.

**Wind pressure**

Wind creates a pressure distribution on the building envelope. The
pressure at a particular façade location is commonly expressed as

$$p_{wind} = C_p \frac{\rho_e v^2}{2}$$

where:

-   $C_p$ is the wind pressure coefficient;
-   $\rho_e$ is the outdoor air density, kg/m³;
-   $v$ is the wind velocity, m/s.

The wind pressure coefficient depends on the building geometry, wind direction, and position on the façade. It can be positive or negative, representing positive or suction pressure relative to the reference outdoor pressure ([Orme & Leksmono, 2002](https://www.aivc.org/sites/default/files/members_area/medias/pdf/Guides/Ventilation%20Modelling%20Data%20Guide_web.pdf), Appendix F).

**Static pressure**

The hydrostatic pressure variation with height results from the weight of the air column. For a constant density, the pressure variation
between a reference height and a height $z$ is

$$p_{static}(z) = -\rho g z$$

where:

- $\rho$ is the air density, kg/m³;
- $g$ is gravitational acceleration, m/s²;
- $z$ is the vertical coordinate measured relative to the reference
    level, m.

When indoor and outdoor air densities differ, the hydrostatic pressure distributions also differ. This produces the [stack effect](https://en.wikipedia.org/wiki/Stack_effect), which can drive infiltration and exfiltration even in the absence of wind.

**Combined pressure sources**

The pressure sources acting on a branch can combine wind and static
pressure contributions. Their algebraic sum is represented by the branch pressure-source $b$. The [sign of each source](https://en.wikipedia.org/wiki/Passive_sign_convention) depends on the orientation assigned to the mass flow of the branch:

- if the flow is from $-$ to $+$ of the source, then the source is positif;
- if the flow is from $+$ to $-$ of the source, then the source is negative.

### Flow sources

Flow sources, $f$, represent imposed mass flows at pressure nodes (such as mechanical ventilation extraction or supply ventilation).

If the mechanical ventilation flow is specified as a volumetric  rate, then it needs to be converted into a mass-flow source using
the indoor or outdoor air density for the extraction and supply, respectively.

The sign convention must be consistent with the mass balance equation: extraction is represented by a positive sink contribution in the balance equation, while supply is represented by a negative sink contribution.


## Graph representation of the airflow network

An airflow network can be represented as an **oriented, weighted
graph**:

-   pressure nodes are the vertices of the graph;
-   airflow paths are the edges (branches) connecting the vertices;
-   each branch has an assigned orientation, which defines the positive direction of airflow;
-   each branch is associated with a flow conductance that depends on its pressure difference and upstream air density.

The graph representation provides a systematic way to formulate the mass
balance equations and the pressure-flow relationships using matrices and
vectors.

### Incidence matrix (A)

The topology of the network is described by the branch-node incidence matrix $A$, with dimensions $B \times N$, where $B$ is the number of branches and $N$ is the number of internal pressure nodes.

The elements of the matrix are defined as follows:

$$A_{ij} =
\begin{cases}
-1 & \text{if branch } i \text{ starts at node } j\\
+1 & \text{if branch } i \text{ ends at node } j\\
0 & \text{if branch } i \text{ is not connected to node } j
\end{cases}
$$

The orientation of a branch is arbitrary but must remain consistent throughout the formulation.

The reference-pressure node is not included among the unknown pressure nodes. Consequently, a branch connecting an internal node to the reference node has only one nonzero entry in its row of $A$.

The incidence matrix is used to express pressure
differences along the branches, $\Delta p = -Ap$, and the mass sum in the nodes, $\sum \dot{m} = A^T \dot{m}$.

### Pressure vector and pressure-source vector

Let $p$ denote the vector of unknown pressures at
the internal nodes, relative to the reference pressure.

The pressure difference across all branches is calculated using

$$\Delta p = -A p + b$$

where:

- $p$ is the vector of internal node pressures, Pa;
- $b$ is the vector of branch pressure sources, Pa;
- $\Delta p$ is the vector of branch pressure differences, Pa.

The vector $b$ includes the wind and static pressure
contributions associated with each branch.

The minus sign in $-A p$ follows the chosen incidence-matrix convention: a branch oriented from node $i$ to node $j$
has a pressure difference $(p_i - p_j)$, before adding the
imposed pressure source.

### Conductance matrix (G)

For a linear flow network, the mass flow through a branch is
proportional to the pressure difference across it. The proportionality coefficient is its mass-flow conductance.

For a nonlinear airflow network, the equivalent conductance depends on the current pressure difference and the upstream air density.

For branch $i$, the equivalent mass-flow conductance is

$$G_i = \rho_{{up},i} K_i
|\Delta p_i|^{n_i-1}$$

The conductance matrix $G$ is diagonal:

$$G = \operatorname{diag}(G_1, G_2, \ldots, G_B)$$

Its diagonal elements have units of kg/(s·Pa).

The mass-flow vector can then be expressed as

$$\dot{m} = G \Delta p$$

Unlike a linear network, the airflow conductance changes with the pressure difference and can also change when the flow
direction reverses and a different upstream density applies.

### Complete matrix formulation

Using the incidence matrix and conductance matrix, the mass conservation equations become

$$A^T G \Delta p + f = 0$$

Substituting the expression for branch pressure differences gives

$$A^T G(-A p + b) + f=0$$

or, equivalently,

$$A^T G A p=A^T G b + f$$

This equation is the fundamental matrix equation used by the fixed-point algorithm.

It is linear in the unknown pressure vector $p$ if $G$ is known. However, in the actual airflow network, $G$ depends on the pressure differences and upstream densities, which themselves depend on the unknown pressures. The original problem is therefore nonlinear.

### Role of the vectors and matrices

The airflow network is completely described by the following quantities:

| Symbol | Description | Dimensions |
|--------|-------------|------------|
| $A$ | Branch-node incidence matrix | $B \times N$ |
| $G$ | Diagonal equivalent conductance matrix | $B \times B$ |
| $b$ | Branch pressure-source vector | $B$ |
| $f$ | Node mass-flow-source vector | $N$ |
| $\rho$ | Internal-node air-density vector | $N$ |
| $\rho_{ref}$ | Reference-node air density | Scalar |
| $p$ | Unknown node pressure vector | $N$ |

The density vector and reference-node density are used to determine the upstream density for each branch.

The matrices $A$ and $G$, together with the vectors $b$ and $f$, provide the mathematical description of the network topology, branch characteristics, and external conditions.


## References

1. [Orme, M., & Leksmono, N. (2002)](https://www.aivc.org/sites/default/files/members_area/medias/pdf/Guides/Ventilation%20Modelling%20Data%20Guide_web.pdf). AIVC Guide 5: Ventilation modelling data guide. International Energy Agency, Air Infiltration Ventilation Center. AIC-GUI, 5.
