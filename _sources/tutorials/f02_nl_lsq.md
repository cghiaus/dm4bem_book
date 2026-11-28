# Airflow networks analysis: nonlinear least-squares method 

The function `dm4bem.flow_network_least_squares()` [nonlinear least-squares method](https://en.wikipedia.org/wiki/Non-linear_least_squares) to analyse the network by determining the node pressures that minimise the mass-balance residuals.

The procedure is:

1. Construct the forward and backward upstream density vectors.
2. Define the pressure-difference function.
3. Define the nonlinear volumetric and mass-flow functions.
4. Define the mass-balance residual function.
5. Scale the residuals to improve numerical conditioning.
6. Call `scipy.optimize.least_squares` to minimise the sum of squared scaled mass-balance residuals.
7. Extract the calculated pressure vector.
8. Check the original, unscaled mass-balance residuals against the prescribed tolerance.
9. Calculate the final volumetric and mass flow rates and return them.

The essential difference from the fixed-point method is in step 6. The fixed-point method constructs an equivalent linear system for the pressures and iteratively updates its conductances. Here, the pressure vector is obtained by directly minimising the nonlinear mass-balance residuals.

The following description focuses on the implementation, assuming that the physical principles and the mathematical formulation of airflow networks have already been presented.

## Definition of the mass-balance residuals

The central element of the numerical method is the residual function:

$$r(p) = A^T \dot{m}(p) + f$$

The objective is to find the pressure vector for which

$$r(p) = 0$$

In other words, the numerical method searches for the node pressures $p$ that satisfy all mass-conservation equations $A^T \dot{m}(p) + f = 0$ simultaneously. Note that $A^T$ plays the role of $\sum$ operator for mass flows.

The function $r(p)$ is the function that the nonlinear optimisation algorithm will attempt to bring to zero.

## Scaling the residuals

Before calling the optimisation algorithm, the code introduces a characteristic mass-flow scale, $s_m$, used to scale the mass-balance residuals:

$$s_m =
\max\left[
\max_{j=1,\ldots, N}|f_j|,
\quad
\max_{i=1,\ldots, B}
\left(
\rho_{\max,i}K_i|b_i|^{n_i}
\right),
\quad
10^{-12}
\right]
$$

where:

- $N$ is the number of internal nodes;
- $B$ is the number of branches;
- $f_j$ is the mass-flow source at node $j$, in kg/s;
- $K_i$ is the flow coefficient of branch $i$;
- $n_i$ is the flow exponent of branch $i$;
- $b_i$ is the pressure source of branch $i$, in Pa;
- $\rho_{\max,i}$ is the maximum of the two possible upstream densities for branch $i$:

$$\rho_{\max,i} =
\max\left(
\rho_{forward,i},
\rho_{backward,i}
\right).$$

The three terms have the following meanings:

1. **Maximum imposed mass flow:**

$$s_1 = \max_j |f_j|.$$

This represents the largest imposed mass-flow source or sink in the network.

2. **Estimated maximum branch mass flow:**

$$s_2 =
\max_i
\left(
\rho_{\max,i}K_i|b_i|^{n_i}
\right).$$

This estimates the largest branch mass flow using the pressure source $b_i$ as a characteristic pressure difference and the maximum possible upstream density.

3. **Lower numerical bound:**

$$s_3 = 10^{-12}\ \mathrm{kg/s}.$$

This prevents the scale from being zero when the other two terms vanish or are extremely small.

The final scale is

$$s_m = \max(s_1, s_2, s_3).$$

The residual vector is then divided by this scale:

$$\widetilde{r}(p) = \frac{r(p)}{s_m}
$$

Since both $r$ and $s_m$ are expressed in kg/s, the scaled residual is dimensionless.

**Note:** $s_2$ is an estimate based on the pressure sources $b_i$, not necessarily the actual maximum branch mass flow. The actual pressure difference is $d_i = -A_i p + b_i$, which also depends on the unknown node pressures. The purpose of $s_m$ is numerical scaling, not establishing a strict upper bound on the flow rates.

The purpose of this scaling is to improve the numerical conditioning of the optimisation problem. It changes the magnitude of the residuals used by the algorithm but does not change the pressure vector that gives an exactly zero residual.

## Nonlinear least-squares solution

The function `scipy.optimize.least_squares` minimises the sum of the squares of the scaled mass-balance residuals:

$$\min_{p} \frac{1}{2} \left\| \widetilde{r}(p) \right\|_2^2.$$

The initial pressure vector `p0` provides the starting point for the optimisation.

The different arguments control how the optimisation is performed.

__Jacobian approximation: `jac='2-point'`__

The Jacobian matrix describes the sensitivity of the mass-balance residuals to changes in the internal node pressures:

$$J(p) = \frac{\partial r(p)}{\partial p}.$$

The option

```python
jac='2-point'
```

instructs SciPy to estimate this matrix numerically using finite differences.

For example, the derivative with respect to pressure component $p_j$ is approximated by

$$\frac{\partial r_i}{\partial p_j} \approx \frac{r_i(p + h_j e_j) - r_i(p)}{h_j},$$

where $h_j$ is a small perturbation and $e_j$ is the unit vector associated with pressure component $p_j$.

Thus, the Jacobian is not explicitly derived or supplied in the Python code. It is estimated automatically by SciPy.

__Variable scaling: `x_scale='jac'`__

The option

```python
x_scale='jac'
```

allows SciPy to scale the pressure variables according to the Jacobian information.

This is useful when the unknown variables have different sensitivities, because it helps the optimisation algorithm choose appropriate pressure updates.

### Convergence tolerances

The three tolerances are:

```python
ftol=1e-12
xtol=1e-12
gtol=1e-12
```

They control different stopping conditions:

- `ftol`: relative change in the cost function, which is the sum of squared scaled residuals;
- `xtol`: relative change in the pressure vector between successive iterations;
- `gtol`: condition based on the gradient of the cost function.

These criteria determine when the optimisation algorithm considers that it has converged.

The argument

```python
max_nfev=max_nfev
```

sets the maximum number of evaluations of the residual function. If this limit is reached before convergence, the solver may terminate without finding an acceptable solution.


## Verification of mass conservation

After the optimisation, the calculated pressure vector is extracted:

```python
p = sol.x
```

The mass-balance residuals are recalculated using the original, unscaled residual function:

```python
mass_balance = residual(p)
```

The maximum absolute residual is then calculated:

```python
err_m = np.linalg.norm(
    mass_balance,
    ord=np.inf)
```

Mathematically,

$$\varepsilon_m = \left\|r(p)\right\|_{\infty} = \max_i |r_i(p)|.$$

The use of the infinity norm means that the largest absolute mass-balance error among all internal nodes is considered.

The solution is accepted only if the optimisation reports success and the maximum mass-balance residual is smaller than the prescribed tolerance:

```python
if not sol.success or err_m > tol_m:
    raise RuntimeError(
        f"Newton solver failed: {sol.message}\n"
        f"Mass-balance residual: {mass_balance} kg/s\n"
        f"Maximum residual: {err_m:.3e} kg/s"
    )
```

The default tolerance is

$$\text{tol}_m = 10^{-9}\ \mathrm{kg/s}.$$

This additional verification is important. The optimisation algorithm can satisfy one of its stopping criteria without necessarily achieving the required accuracy in the physical mass-balance equations.

The code therefore checks the actual, unscaled mass-balance residuals before accepting the solution.

## Calculation of the final flow rates

Once the pressure vector $p$ has passed the mass-conservation check, the final pressure difference across each branch is calculated as

$$d = -A p + b$$

The volumetric flow rate in each branch is then calculated using the power-law relationship:

$$\dot{V} = K \circ
\operatorname{sign}(d) \circ
|d|^n $$

where $\circ$ denotes [element-wise multiplication](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)) and the exponentiation is also performed element-wise.

The upstream density for each branch is selected according to the direction of flow:

$$\rho_{up} =
\begin{cases}
\rho_{forward}, & \text{if } d_i \geq 0,\\
\rho_{backward}, & \text{if } d_i < 0.
\end{cases}$$

Here, the selection is performed independently for each branch $i$.

Finally, the branch mass flow rates are calculated from the volumetric flow rates:

$$\dot{m} = \rho_{up} \circ \dot{V}$$

The function returns three vectors:

- $p$: internal node pressures, in Pa;
- $\dot{m}$: branch mass flow rates, in kg/s;
- $\dot{V}$: branch volumetric flow rates, in m³/s.


## References

Gommers, R., Virtanen, P., Haberland, M., Burovski, E., Reddy, T., Weckesser, W., ... & Striega, K. (2024). scipy/scipy: SciPy 1.15. 0. Zenodo.
