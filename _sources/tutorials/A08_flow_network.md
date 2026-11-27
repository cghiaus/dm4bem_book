# Airflow networks

## 1. Physical principles: balance equations and constitutive laws

An airflow network represents the movement of air through openings, cracks, ducts, or other flow paths connecting different pressure nodes. The network is governed by two fundamental types of equations:

1.  **Balance equations**, expressing conservation of mass at the
    pressure nodes.
2.  **Constitutive laws**, describing the relationship between the
    airflow rate through a branch and the pressure difference across it.

These equations form a nonlinear system because the airflow through an opening generally depends nonlinearly on the pressure difference.

### 1.1. Mass conservation at pressure nodes

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

### 1.2. Nonlinear constitutive law: volumetric flow and pressure difference

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

### 1.3. Relationship between mass flow and volumetric flow

Mass flow rate and volumetric flow rate are related through air density:

$$\dot{m}_{ij} = \rho_{ij} \dot{V}_{ij}$$

where:

- $\dot{m}_{ij}$ is the  mass flow, kg/s; 
- $\rho_{ij}$ is the air density, kg/m³;
- $\dot{V}_{ij}$ the volumetric flow rate, m³/s.

If the flow coefficient is expressed in
m³/h, then

$$\dot{m}_{ij} =
\frac{\rho_{ij}}{3600} \dot{V}_{ij}$$

Substituting the power-law constitutive relationship gives

$$\dot{m}_{ij} = \frac{\rho_{ij}K_{ij}}{3600} \operatorname{sign} (\Delta{p}) |\Delta p_{ij}| ^{n_{ij}}$$

Note that air density depends on the location from which the air originates.

### 1.4. Upstream air density

In an airflow network, the density used to calculate mass flow must correspond to the upstream node of the branch, that is, the node from which the air enters the branch.

For a branch oriented from node $i$ to node $j$:

-   if $\Delta p_{ij}$ produces positive flow, air travels
    from $i$ to $j$, and the density is $\rho_i$;
-   if $\Delta p_{ij}$ produces negative flow, air travels
    from $j$ to $i$, and the density is $\rho_j$.

Therefore, the mass flow relationship becomes

$$\dot{m}_{ij} = \frac{K_{ij}}{3600}
\rho_{up,ij} \operatorname{sign} (\Delta p_{ij})
|\Delta p_{ij}|^{n_{ij}}$$

where $\rho_{up,ij}$ is the density at the
upstream end of the branch.

This distinction is particularly important when a branch connects  environments with different air densities.

In the Python implementation, the upstream density is selected according to the sign of the branch pressure difference. The code stores two possible densities for each branch, `rho_forward` and `rho_backward`, corresponding to the two possible flow directions, and selects the appropriate one during each iteration.

## 2. Pressure sources and flow sources

An airflow is driven by imposed pressure differences and imposed mass flows.

### 2.1. Pressure sources

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

### 2.2. Flow sources

Flow sources, $f$, represent imposed mass flows at pressure nodes (such as mechanical ventilation extraction or supply ventilation).

If the mechanical ventilation flow is specified as a volumetric  rate, then it needs to be converted into a mass-flow source using
the indoor or outdoor air density for the extraction and supply, respectively.

The sign convention must be consistent with the mass balance equation: extraction is represented by a positive sink contribution in the balance equation, while supply is represented by a negative sink contribution.


## 3. Graph representation of the airflow network

An airflow network can be represented as an **oriented, weighted
graph**:

-   pressure nodes are the vertices of the graph;
-   airflow paths are the edges (branches) connecting the vertices;
-   each branch has an assigned orientation, which defines the positive direction of airflow;
-   each branch is associated with a flow conductance that depends on its pressure difference and upstream air density.

The graph representation provides a systematic way to formulate the mass
balance equations and the pressure-flow relationships using matrices and
vectors.

### 3.1. Incidence matrix (A)

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

### 3.2. Pressure vector and pressure-source vector

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

### 3.3. Conductance matrix (G)

For a linear flow network, the mass flow through a branch is
proportional to the pressure difference across it. The proportionality coefficient is its mass-flow conductance.

For a nonlinear airflow network, the equivalent conductance depends on the current pressure difference and the upstream air density.

For branch $i$, the equivalent mass-flow conductance is

$$G_i = \frac{\rho_{{up},i} K_i}{3600}
|\Delta p_i|^{n_i-1}$$

The conductance matrix $G$ is diagonal:

$$G = \operatorname{diag}(G_1, G_2, \ldots, G_B)$$

Its diagonal elements have units of kg/(s·Pa).

The mass-flow vector can then be expressed as

$$\dot{m} = G \Delta p$$

Unlike a linear network, the airflow conductance changes with the pressure difference and can also change when the flow
direction reverses and a different upstream density applies.

### 3.4. Complete matrix formulation

Using the incidence matrix and conductance matrix, the mass conservation equations become

$$A^T G \Delta p + f = 0$$

Substituting the expression for branch pressure differences gives

$$A^T G(-A p + b) + f=0$$

or, equivalently,

$$A^T G A p=A^T G b + f$$

This equation is the fundamental matrix equation used by the fixed-point algorithm.

It is linear in the unknown pressure vector $p$ if $G$ is known. However, in the actual airflow network, $G$ depends on the pressure differences and upstream densities, which themselves depend on the unknown pressures. The original problem is therefore nonlinear.

### 3.5. Role of the vectors and matrices

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


## 4. Fixed-point method for solving the nonlinear airflow network

### 4.1. Principle of the method

The difficulty in solving the airflow network comes from the dependence
of the conductance matrix (G) on the unknown pressure differences.

The fixed-point method addresses this difficulty by repeatedly solving a
linear system in which the conductances are calculated using the
pressure estimates from the preceding iteration.

The procedure is:

1.  Assume initial values for the internal node pressures.
2.  Calculate the pressure differences across all branches.
3.  Determine the flow direction and upstream air density for each
    branch.
4.  Calculate the equivalent conductance of each branch.
5.  Assemble and solve the linear pressure system.
6.  Update the pressure estimates, possibly using under-relaxation.
7.  Calculate the exact nonlinear mass flows and check convergence.
8.  Repeat until the pressure and mass-balance convergence criteria are
    satisfied.

The method transforms the nonlinear problem into a sequence of linear
problems.

### 4.2. Initialisation

The algorithm starts with an initial estimate of the internal node
pressures:

\[ `\boldsymbol `{=tex}p\^{(0)}=`\boldsymbol `{=tex}p_0 \]

The initial pressures can be estimated from physical considerations,
previous calculations, or simply assigned reasonable values.

The branch orientations are defined by the incidence matrix (A).

For each branch, the code identifies the two possible upstream
densities:

-   Density at the starting node for positive flow;
-   Density at the ending node for negative flow.

If one end of the branch is the reference-pressure node, the density at
that end is taken from (`\rho`{=tex}\_{`\mathrm{ref}`{=tex}}).

These two densities are stored before starting the iterations.

### 4.3. Calculation of pressure differences

At iteration (k), the branch pressure differences are calculated from
the current pressure estimates:

\[ `\Delta`{=tex}`\boldsymbol `{=tex}p\^{(k)} =
-A`\boldsymbol `{=tex}p\^{(k)}+`\boldsymbol `{=tex}b \]

The sign of each component determines the direction of airflow.

A positive pressure difference corresponds to flow in the assigned
branch orientation, whereas a negative pressure difference corresponds
to flow in the opposite direction.

### 4.4. Selection of upstream density

For each branch, the upstream density is selected according to the sign
of its pressure difference:

\[ `\rho`{=tex}\_{`\mathrm{up}`{=tex},i}\^{(k)} =
```{=tex}
\begin{cases}
\rho_{\mathrm{forward},i} &
\Delta p_i^{(k)}\geq 0\\
\rho_{\mathrm{backward},i} &
\Delta p_i^{(k)}<0
\end{cases}
```
\]

This operation is important because the density used to calculate mass
flow changes when the airflow reverses.

The density selection is implemented in Python using:

``` python
rho_up = np.where(d >= 0, rho_forward, rho_backward)
```

### 4.5. Calculation of equivalent conductances

Using the pressure differences and upstream densities obtained at
iteration (k), the equivalent conductance of each branch is calculated
as

\[ G_i\^{(k)} = `\frac{\rho_{\mathrm{up},i}^{(k)}K_i}{3600}`{=tex}
`\left`{=tex}\|`\Delta `{=tex}p_i^{(k)}`\right`{=tex}\|^{n_i-1} \]

For numerical stability, the implementation introduces a minimum
pressure difference (d\_{`\mathrm{floor}`{=tex}}):

\[ G_i\^{(k)} = `\frac{\rho_{\mathrm{up},i}^{(k)}K_i}{3600}`{=tex}
`\left[
\max\left(
|\Delta p_i^{(k)}|,d_{\mathrm{floor}}
\right)
\right]`{=tex}\^{n_i-1} \]

This avoids excessively large conductances when the pressure difference
approaches zero, particularly for flow exponents smaller than one.

The resulting diagonal conductance matrix is

\[ G\^{(k)} = `\operatorname{diag}`{=tex} `\left`{=tex}(
G_1^{(k)},`\ldots`{=tex},G_B^{(k)} `\right`{=tex}) \]

The conductance matrix is therefore updated at every iteration.

### 4.6. Solution of the linear pressure system

Once the conductance matrix is known, the nonlinear problem is
approximated by a linear problem.

The mass balance equations are written as

\[ A^TG^{(k)}A`\boldsymbol `{=tex}p\_{`\mathrm{linear}`{=tex}}\^{(k+1)}
= A^TG^{(k)}`\boldsymbol `{=tex}b+`\boldsymbol `{=tex}f \]

The new pressure estimate is obtained by solving this system:

\[ `\boldsymbol `{=tex}p\_{`\mathrm{linear}`{=tex}}\^{(k+1)} =
`\left`{=tex}(A^TG^{(k)}A`\right`{=tex})\^{-1}
`\left`{=tex}(A^TG^{(k)}`\boldsymbol `{=tex}b+`\boldsymbol `{=tex}f`\right`{=tex})
\]

In the Python implementation, the inverse is not explicitly calculated.
Instead, the linear system is solved using:

``` python
M = A.T @ G @ A
rhs = A.T @ G @ b + f

p_linear = np.linalg.solve(M, rhs)
```

This is preferable numerically to explicitly calculating the inverse
matrix.

### 4.7. Under-relaxation

Directly replacing the previous pressure estimate with the solution of
the linear system can sometimes produce oscillations or convergence
difficulties.

To improve numerical stability, the algorithm uses under-relaxation:

\[ `\boldsymbol `{=tex}p\^{(k+1)} =
(1-`\omega`{=tex})`\boldsymbol `{=tex}p\^{(k)} +
`\omega`{=tex}`\boldsymbol `{=tex}p\_{`\mathrm{linear}`{=tex}}\^{(k+1)}
\]

where (`\omega`{=tex}) is the relaxation coefficient:

-   (0\<`\omega`{=tex}\<1): under-relaxation;
-   (`\omega=1`{=tex}): no relaxation;
-   Smaller values of (`\omega`{=tex}) produce more conservative
    pressure updates.

The code uses (`\omega=0.5`{=tex}) by default.

Under-relaxation can improve convergence when successive pressure
estimates oscillate, although excessively small values can make
convergence unnecessarily slow.

### 4.8. Calculation of exact nonlinear mass flows

After updating the pressures, the algorithm calculates the mass flow
rates using the original nonlinear constitutive relationship, rather
than the equivalent linear conductances.

For each branch:

\[ `\dot `{=tex}m_i\^{(k+1)} =
`\frac{\rho_{\mathrm{up},i}^{(k+1)}K_i}{3600}`{=tex}
`\operatorname{sign}`{=tex}(`\Delta `{=tex}p_i\^{(k+1)})
\|`\Delta `{=tex}p_i^{(k+1)}\|^{n_i} \]

The use of the exact nonlinear relationship is essential for checking
whether the solution satisfies the original airflow network equations.

The mass-balance residual is then calculated:

\[ `\boldsymbol `{=tex}r_m\^{(k+1)} =
A^T`\dot{\boldsymbol m}`{=tex}^{(k+1)} +`\boldsymbol `{=tex}f \]

The maximum absolute residual is

\[ `\varepsilon`{=tex}*m\^{(k+1)} =
\|`\boldsymbol `{=tex}r_m\^{(k+1)}\|*`\infty`{=tex} \]

This measures the largest mass imbalance among all internal pressure
nodes.

### 4.9. Convergence criteria

The algorithm uses two independent convergence criteria.

**Pressure convergence**

The relative change in pressure between two successive iterations is
calculated as

\[ `\varepsilon`{=tex}\_p\^{(k+1)} = `\frac{
\|\boldsymbol p^{(k+1)}-\boldsymbol p^{(k)}\|_\infty
}{
\max\left(
p_{\mathrm{scale}},
\|\boldsymbol p^{(k+1)}\|_\infty
\right)
}`{=tex} \]

The parameter (p\_{`\mathrm{scale}`{=tex}}) prevents the relative
pressure error from becoming excessively restrictive when all node
pressures are close to zero.

**Mass conservation convergence**

The mass-balance residual must satisfy

\[ `\varepsilon`{=tex}\_m\^{(k+1)}\<`\mathrm{tol}`{=tex}\_m \]

The pressure change must also satisfy

\[ `\varepsilon`{=tex}\_p\^{(k+1)}\<`\mathrm{tol}`{=tex}\_p \]

The algorithm terminates only when both conditions are satisfied:

\[ `\varepsilon`{=tex}\_p\<`\mathrm{tol}`{=tex}\_p
`\quad`{=tex}`\text{and}`{=tex}`\quad`{=tex}
`\varepsilon`{=tex}\_m\<`\mathrm{tol}`{=tex}\_m \]

This dual criterion is important. A small pressure change does not
necessarily mean that the nonlinear mass-balance equations are
satisfied. Checking both criteria ensures that the solution is
sufficiently stable and conserves mass.

If the convergence criteria are not satisfied within the maximum number
of iterations, the algorithm raises an error and reports the final
pressure error and mass-balance residual.

## 5. Summary of the algorithm

The complete fixed-point procedure can be summarised as follows.

``` text
INPUT:
    A, b, f, K, n, rho, rho_ref
    Initial pressures p0
    Relaxation coefficient omega
    Convergence tolerances tol_p, tol_m

INITIALISATION:
    p = p0
    Determine forward and backward densities
    for each branch

ITERATION:

    1. Calculate branch pressure differences
       d = -A @ p + b

    2. Determine upstream density
       rho_up = upstream_density(d)

    3. Calculate equivalent conductances
       g = rho_up * K / 3600
           * max(abs(d), d_floor)**(n - 1)

    4. Assemble the linear system
       M = A.T @ G @ A
       rhs = A.T @ G @ b + f

    5. Solve for new pressures
       p_linear = solve(M, rhs)

    6. Apply under-relaxation
       p_new = (1 - omega)*p + omega*p_linear

    7. Check pressure convergence
       err_p = relative_pressure_change(p_new, p)

    8. Calculate exact nonlinear mass flows
       m_dot = mass_flow(p_new)

    9. Check mass conservation
       residual = A.T @ m_dot + f
       err_m = max(abs(residual))

   10. Update pressures
       p = p_new

   11. Stop if both convergence criteria
       are satisfied.

OUTPUT:
    Node pressures p
    Branch mass flows m_dot
    Branch volumetric flows V_dot
```

Once convergence is achieved, the volumetric flow rates are calculated
using the original power-law relationship:

\[ `\dot `{=tex}V_i = K_i`\operatorname{sign}`{=tex}(`\Delta `{=tex}p_i)
\|`\Delta `{=tex}p_i\|\^{n_i} \]

The function returns three quantities:

-   (`\boldsymbol `{=tex}p): internal node pressures (Pa);
-   (`\dot{\boldsymbol m}`{=tex}): branch mass flow rates (kg/s);
-   (`\dot{\boldsymbol V}`{=tex}): branch volumetric flow rates (m³/h).

## 6. Interpretation and characteristics of the method

The fixed-point method can be interpreted as a successive linearisation
of the nonlinear airflow network.

At each iteration, the current pressure differences determine the
equivalent conductances. These conductances define a linear network,
which is solved to obtain updated node pressures. The procedure is
repeated until the linearised solution is consistent with the original
nonlinear constitutive laws and the mass conservation equations.

Its main characteristics are:

1.  **Physical formulation:** Mass conservation is enforced at every
    internal pressure node, while the constitutive laws describe the
    nonlinear behaviour of individual airflow paths.

2.  **Graph-based formulation:** The incidence matrix describes the
    network topology, while the diagonal conductance matrix describes
    the current hydraulic (or, more precisely, pneumatic) properties of
    its branches.

3.  **Variable density:** The upstream density is updated according to
    the direction of airflow, allowing the method to account for
    different indoor and outdoor air densities.

4.  **Linear systems:** Each iteration requires solving a linear system,
    rather than solving the entire nonlinear system simultaneously.

5.  **Under-relaxation:** The relaxation coefficient provides a means of
    controlling the pressure updates and improving numerical stability.

6.  **Physical convergence check:** The algorithm verifies both pressure
    convergence and mass conservation using the exact nonlinear
    mass-flow equations.

The method is particularly useful for building airflow networks
involving infiltration, exfiltration, stack effect, wind pressure, and
mechanical ventilation.

**An important numerical consideration:** The method is a fixed-point
iteration, not a Newton method. Its convergence is not guaranteed for
every network or every initial pressure estimate. The choice of
relaxation coefficient, initial pressures, and treatment of very small
pressure differences can influence convergence. Furthermore, the
pressure matrix must be nonsingular, which requires the internal nodes
to be adequately connected to a reference-pressure node through the
network.

------------------------------------------------------------------------

*Source: Description based on the uploaded Python implementation,
`fix_point_up.py`, particularly the function
`flow_network_fixed_point()` and its example airflow network.*

## References

1. [Orme, M., & Leksmono, N. (2002)](https://www.aivc.org/sites/default/files/members_area/medias/pdf/Guides/Ventilation%20Modelling%20Data%20Guide_web.pdf). AIVC Guide 5: Ventilation modelling data guide. International Energy Agency, Air Infiltration Ventilation Center. AIC-GUI, 5.
