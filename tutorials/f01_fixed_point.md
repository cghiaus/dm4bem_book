# Airflow networks analysis: fixed-point method 

## Principle of the method

The difficulty in solving the airflow network comes from the dependence of the conductance matrix $G$ on the unknown pressure differences.

The fixed-point method addresses this difficulty by repeatedly solving a linear system in which the conductances are calculated using the pressure estimates from the preceding iteration.

The procedure is:

1. Assume initial values for the internal node pressures.
2. Construct the forward and backward upstream-density vectors.
3. Define the pressure-difference function.
4. Define the upstream air density for each branch.
5. Define the mass-balance.
6. Iterations:
    - Calculate pressure differences, select upstream density, calculate conductances;
    - solve the linear system for pressure;
    - under-relaxation update of pressures;
    - calculate the mass flows;
    - check convergence criteria;
    - repeat until the pressure and mass-balance convergence criteria are satisfied.

## Initialisation

The algorithm starts with an initial estimate of the node pressures:

$$p^{(0)} = p_0$$

The initial pressures can be estimated from physical considerations, previous calculations (e.g., by using a linear lodel), or simply assigned reasonable values.

## Forward and backward upstream density

The branch orientations are defined by the incidence matrix $A$.

For each branch, the code identifies the two possible upstream densities:

- density at the starting node for positive flow, `rho_forward`;
- density at the ending node for negative flow, `rho_backward`.

If one end of the branch is the reference-pressure node, the density at that end is taken from $\rho_{ref}$.

These two densities are stored before starting the iterations.

## Pressure differences

At iteration $(k)$, the branch pressure differences are calculated from the current pressure estimates:

$$\Delta p^{(k)} = -A p^{(k)} + b$$

The sign of each component determines the direction of airflow.

A positive pressure difference corresponds to flow in the assigned branch orientation, whereas a negative pressure difference corresponds to flow in the opposite direction.

## Upstream densities

For each branch, the upstream density is selected according to the sign of its pressure difference:

$$\rho_{up,i}^{(k)} =
\begin{cases}
\rho_{forward,i} & \text{if } \Delta p_i^{(k)} \geq 0\\
\rho_{backward,i} & \text{if } \Delta p_i^{(k)} < 0
\end{cases}$$

This operation is important because the density used to calculate mass flow changes when the airflow reverses.

The density selection is implemented in Python using:

```python
rho_up = np.where(d >= 0, rho_forward, rho_backward)
```

## Equivalent conductances

Using the pressure differences and upstream densities obtained at iteration $k$, the equivalent conductance of each branch is calculated as

$$G_i^{(k)} = \rho_{up,i}^{(k)}K_i \left|\Delta p_i^{(k)}\right|^{n_i-1}$$

For numerical stability, the implementation introduces a minimum pressure difference $d_{floor}$:

$$G_i^{(k)} = \rho_{up,i}^{(k)}K_i \left[\max\left(|\Delta p_i^{(k)}|, d_{floor}\right)\right]^{n_i-1}$$

This avoids excessively large conductances when the pressure difference approaches zero, particularly for flow exponents smaller than one.

The resulting diagonal conductance matrix is

$$G^{(k)} = \operatorname{diag}\left(G_1^{(k)}, \ldots, G_B^{(k)}\right)$$

The conductance matrix is therefore updated at every iteration.

## Solution of the linear pressure system

Once the conductance matrix is known, the nonlinear problem is approximated by a linear problem.

The mass balance equations are written as

$$A^T G^{(k)} A p_{linear}^{(k+1)} = A^T G^{(k)} b + f$$

The new pressure estimate is obtained by solving this system:

$$p_{linear}^{(k+1)} = \left(A^T G^{(k)} A\right)^{-1} \left(A^T G^{(k)} b + f \right)$$

In the Python implementation, the inverse is not explicitly calculated. Instead, the linear system is solved using:

```python
M = A.T @ G @ A
rhs = A.T @ G @ b + f
p_linear = np.linalg.solve(M, rhs)
```
This is preferable numerically to explicitly calculating the inverse matrix.

## Under-relaxation update of pressures

Directly replacing the previous pressure estimate with the solution of the linear system can sometimes produce oscillations or convergence difficulties.

To improve numerical stability, the algorithm uses under-relaxation:

$$p^{(k+1)} = (1 - \omega) p^{(k)} + \omega p_{linear}^{(k+1)}$$

where $\omega$ is the relaxation coefficient:

- $0 < \omega < 1$: under-relaxation;
- $\omega = 1$: no relaxation;

Smaller values of $\omega$ produce more conservative pressure updates.

The code uses $\omega = 0.5$ by default.

Under-relaxation can improve convergence when successive pressure estimates oscillate, although excessively small values can make convergence unnecessarily slow.

## Exact nonlinear mass flows

After updating the pressures, the algorithm calculates the mass flow rates using the original nonlinear constitutive relationship, rather than the equivalent linear conductances.

For each branch:

$$\dot{m}_i^{(k+1)} = \frac{\rho_{up,i}^{(k+1)}K_i}{3600} \operatorname{sign}(\Delta p_i^{(k+1)}) |\Delta p_i^{(k+1)}|^{n_i}$$

The use of the exact nonlinear relationship is essential for checking whether the solution satisfies the original airflow network equations.

The mass-balance residual is then calculated:

$$r_m^{(k+1)} = A^T \dot{m}^{(k+1)} + f$$

The maximum absolute residual is

$$\varepsilon_m^{(k+1)} = \|r_m^{(k+1)}\|_\infty$$

This measures the largest mass imbalance among all internal pressure nodes.

## Convergence criteria

The algorithm uses two independent convergence criteria.

**Pressure convergence**

The relative change in pressure between two successive iterations is calculated as

$$\varepsilon_p^{(k+1)} = \frac{\|p^{(k+1)} - p^{(k)}\|_\infty}{\max\left(p_{scale}, \|p^{(k+1)}\|_\infty\right)}$$

The parameter $p_{scale}$ prevents the relative pressure error from becoming excessively restrictive when all node pressures are close to zero.

**Mass conservation convergence**

The mass-balance residual must satisfy

$$\varepsilon_m^{(k+1)} < \mathrm{tol}_m$$

The pressure change must also satisfy

$$\varepsilon_p^{(k+1)} < \mathrm{tol}_p$$

The algorithm terminates only when both conditions are satisfied:

$$\varepsilon_p < \mathrm{tol}_p \quad \text{and} \quad \varepsilon_m < \mathrm{tol}_m$$

This dual criterion is important. A small pressure change does not necessarily mean that the nonlinear mass-balance equations are satisfied. Checking both criteria ensures that the solution is sufficiently stable and conserves mass.

If the convergence criteria are not satisfied within the maximum number of iterations, the algorithm raises an error and reports the final pressure error and mass-balance residual.

## Summary of the algorithm

The complete fixed-point procedure can be summarised as follows.

```text
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
       g = rho_up * K * max(abs(d), d_floor)**(n - 1)

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

