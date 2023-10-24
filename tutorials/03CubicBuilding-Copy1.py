#!/usr/bin/env python
# coding: utf-8

# # Toy model building
# 
# [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem_book/main?labpath=%2Ftutorials%2F03CubicBuilding.ipynb)
# 
# This tutorial presents the theory of heat transfer in buildings and examplifies it on a [toy model](https://en.m.wikipedia.org/wiki/Toy_model). The model is a thermal circuit with thermal capacities in some of the nodes.
# 
# **Objectives:**
# - Analyse a cubic building with 5 identical walls & a transparent wall (glass window), air infiltration, and HVAC system controlling the indoor air temperature.
# - Model the heat transfer in the building by a thermal circuit.
# - Obtain the mathematical model as a system of Differential Algebraic Equations (DAE) from the thermal circuit.
# - Transfrom the system of DAE into state-space representation.
# - Find the steady-state solution.
# - Simulate by using Euler methods for numerical integration..

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem


# ## Physical analysis
# 
# ### Description of the building
# 
# ![cube](../figures/03_cube_principle.svg)
# > Figure 1. Simple ventilated room (5 two-layer walls and 1 glass window) equiped with an HVAC system which acts as a proportional controller.
# 
# Let’s consider a cubic building with an HVAC systems acting as a [proportional controller](https://en.m.wikipedia.org/wiki/Proportional_control).

# The dimensions and surface areas of the building are:
# - $l=3 \: \mathrm{m}$ - edge length of the cube;
# - $S_g=l^2$   - surface area of the glass window;
# - $S_c = S_i = 5 \times S_g$   - surface area of the 5 (concrete and insulation) walls.

# In[2]:


l = 3               # m length of the cubic room
Sg = l**2           # m² surface of the glass wall
Sc = Si = 5 * Sg    # m² surface of concrete & insulation of the 5 walls


# ### Thermo-physical properties
# The thermophysical properties of the air (in SI units) are:

# In[3]:


air = {'Density': 1.2,                      # kg/m³
       'Specific heat': 1000}               # J/(kg·K)
# pd.DataFrame.from_dict(air, orient='index', columns=['air'])
pd.DataFrame(air, index=['Air'])


# The [thermophysical properties](https://energieplus-lesite.be/donnees/enveloppe44/enveloppe2/conductivite-thermique-des-materiaux/) ([thermal conductivities](https://en.m.wikipedia.org/wiki/List_of_thermal_conductivities), [densities](https://en.wikipedia.org/wiki/Density) and [specific heat capacities](https://en.m.wikipedia.org/wiki/Table_of_specific_heat_capacities)) and the geometry (widths and surface areas) of the three materials (i.e., concrete, insulation, glass) in SI units are:

# In[4]:


concrete = {'Conductivity': 1.400,
            'Density': 2300.0,
            'Specific heat': 880,
            'Width': 0.2,
            'Surface': 5 * l**2}

insulation = {'Conductivity': 0.027,
              'Density': 55.0,
              'Specific heat': 1210,
              'Width': 0.08,
              'Surface': 5 * l**2}

glass = {'Conductivity': 1.4,
         'Density': 2500,
         'Specific heat': 1210,
         'Width': 0.04,
         'Surface': l**2}

wall = pd.DataFrame.from_dict({'Layer_out': concrete,
                               'Layer_in': insulation,
                               'Glass': glass},
                              orient='index')
wall


# where `Slices` is the number of meshes used in the [numerical discretization](https://en.m.wikipedia.org/wiki/Discretization). 

# ### Radiative properties
# 
# The [radiative properties](https://en.wikipedia.org/wiki/Emissivity#Absorptivity) of the surfaces are:
# - long wave [emmisivity](https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html) of concrete (between normal and rough) and glass pyrex;
# - short wave [absortivity of solar radiation](https://www.engineeringtoolbox.com/solar-radiation-absorbed-materials-d_1568.html) of white smooth surfaces;
# - short wave [transmittance](https://www.engineeringtoolbox.com/optical-properties-glazing-materials-d_1355.html) of window glass (thickness of 4 mm);
# - short wave [absortivity and transmittance](https://energieplus-lesite.be/techniques/enveloppe7/composants-de-l-enveloppe/vitrages/vitrage-permettant-le-controle-solaire/) of reflective blue window glass.

# In[5]:


# radiative properties
ε_wLW = 0.85    # long wave emmisivity: wall surface (concrete)
ε_gLW = 0.90    # long wave emmisivity: glass pyrex
α_wSW = 0.25    # short wave absortivity: white smooth surface
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass


# The [Stefan-Boltzmann constant](https://en.m.wikipedia.org/wiki/Stefan–Boltzmann_constant) is:

# In[6]:


σ = 5.67e-8     # W/(m²⋅K⁴) Stefan-Bolzmann constant
print(f'σ = {σ} W/(m²⋅K⁴)')


# ### Convection coefficients
# 
# Conventional values for the [convection coeficients](https://energieplus-lesite.be/theories/enveloppe9/echanges-chaleur-parois/resistance-thermique-d-echange-superficiel/) for indoor and outdoor convection in W/(m²⋅K) are:

# In[7]:


h = pd.DataFrame([{'in': 8., 'out': 25}], index=['h'])  # W/(m²⋅K)


# In[8]:


h


# In[9]:


# conduction
G_cd = wall['Conductivity'] / wall['Width'] * wall['Surface']
pd.DataFrame(G_cd, columns={'Conductance'})



# In[10]:


# convection
Gw = h * wall['Surface'][0]     # wall
Gg = h * wall['Surface'][2]     # glass


# In[11]:


# view factor wall-glass
Fwg = glass['Surface'] / concrete['Surface']



# In[12]:


# long wave radiation
Tm = 20 + 273   # K, mean temp for radiative exchange

GLW1 = 4 * σ * Tm**3 * ε_wLW / (1 - ε_wLW) * wall['Surface']['Layer_in']
GLW12 = 4 * σ * Tm**3 * Fwg * wall['Surface']['Layer_in']
GLW2 = 4 * σ * Tm**3 * ε_gLW / (1 - ε_gLW) * wall['Surface']['Glass']


# The equivalent conductance, in W/K, for the radiative long-wave heat exchange between the wall and the glass window is:
# 
# $$G = \frac{1}{1/G_1 + 1/G_{1,2} + 1/G_2}$$

# In[13]:


GLW = 1 / (1 / GLW1 + 1 / GLW12 + 1 / GLW2)


# *Note*: Resistances in [series or parallel](https://en.m.wikipedia.org/wiki/Series_and_parallel_circuits) can be replaced by their equivalent resistance. 

# #### Advection
# 
# The [volumetric flow rate](https://en.m.wikipedia.org/wiki/Volumetric_flow_rate) of the air, in m³/s, is:
# 
# $$\dot{V}_a = \frac{\mathrm{ACH}}{3600} V_a$$
# 
# where:
# - $\mathrm{ACH}$  ([air changes per hour](https://en.m.wikipedia.org/wiki/Air_changes_per_hour)) is the air infiltration rate, 1/h;
# - $3600$ - number of seconds in one hour, s/h;
# - $V_a$ - volume of the air in the thermal zone, m³.

# In[14]:


# ventilation flow rate
Va = l**3                   # m³, volume of air
ACH = 1                     # air changes per hour
Va_dot = ACH / 3600 * Va    # m³/s, air infiltration


# In[15]:


# ventilation & advection
Gv = air['Density'] * air['Specific heat'] * Va_dot


# In[16]:


# P-controler gain
Kp = 1e4            # almost perfect controller Kp -> ∞
Kp = 1e-3           # no controller Kp -> 0
Kp = 0



# In[17]:


# glass: convection outdoor & conduction
Ggs = float(1 / (1 / Gg['out'] + 1 / (2 * G_cd['Glass'])))


# ### Thermal capacities
# #### Walls
# The [thermal capacities](https://en.m.wikipedia.org/wiki/Heat_capacity) of the wall, in J/kg, are:
# 
# $$C_w= m_w c_w= \rho_w c_w w_w S_w$$
# 
# where:
# - $m_w = \rho_w w_w S_w$ is the mass of the wall, kg;
# - $c_w$ - [specific heat capacity](https://en.m.wikipedia.org/wiki/Specific_heat_capacity), J/(kg⋅K);
# - $\rho_w$ - [density](https://en.m.wikipedia.org/wiki/Density), kg/m³;
# - $w_w$ - width of the wall, m;
# - $S_w$ - surface area of the wall, m².

# In[18]:


C = wall['Density'] * wall['Specific heat'] * wall['Surface'] * wall['Width']



# In[19]:


C['Air'] = air['Density'] * air['Specific heat'] * Va
pd.DataFrame(C, columns={'Capacity'})



# In[20]:


A = np.zeros([12, 8])       # n° of branches X n° of nodes
A[0, 0] = 1                 # branch 0: -> node 0
A[1, 0], A[1, 1] = -1, 1    # branch 1: node 0 -> node 1
A[2, 1], A[2, 2] = -1, 1    # branch 2: node 1 -> node 2
A[3, 2], A[3, 3] = -1, 1    # branch 3: node 2 -> node 3
A[4, 3], A[4, 4] = -1, 1    # branch 4: node 3 -> node 4
A[5, 4], A[5, 5] = -1, 1    # branch 5: node 4 -> node 5
A[6, 4], A[6, 6] = -1, 1    # branch 6: node 4 -> node 6
A[7, 5], A[7, 6] = -1, 1    # branch 7: node 5 -> node 6
A[8, 7] = 1                 # branch 8: -> node 7
A[9, 5], A[9, 7] = 1, -1    # branch 9: node 5 -> node 7
A[10, 6] = 1                # branch 10: -> node 6
A[11, 6] = 1                # branch 11: -> node 6


# np.set_printoptions(suppress=False)
# pd.DataFrame(A)



# In[21]:


G = np.array(np.hstack(
    [Gw['out'],
     2 * G_cd['Layer_out'], 2 * G_cd['Layer_out'],
     2 * G_cd['Layer_in'], 2 * G_cd['Layer_in'],
     GLW,
     Gw['in'],
     Gg['in'],
     Ggs,
     2 * G_cd['Glass'],
     Gv,
     Kp]))

# np.set_printoptions(precision=3, threshold=16, suppress=True)
# pd.set_option("display.precision", 1)
# pd.DataFrame(G)


# In[22]:


neglect_air_glass = False

if neglect_air_glass:
    C = np.array([0, C['Layer_out'], 0, C['Layer_in'], 0, 0,
                 0, 0])
else:
    C = np.array([0, C['Layer_out'], 0, C['Layer_in'], 0, 0,
                 C['Air'], C['Glass']])

# pd.set_option("display.precision", 3)
# pd.DataFrame(C)


# ### b: temperature source vector
# 
# The vector of *temperature sources* is $b$, of size $n_q$, the number of branches (in this example 12). An element of the vector $b$ corresponding to a branch without a source is zero. If the flow in a source is from the low potential to the high potential of the source (i.e. from - to +), then the source is positive. If the flow rate in the temperature source is from high potential to low potential (i.e. from + to -), then the source is negative (see [passive sign convention](https://en.m.wikipedia.org/wiki/Passive_sign_convention)). 
# 
# For the thermal circuit shown in Figure 3,
# 
# $$b = [\begin{matrix}
# T_o &0  &0  &0  &0  &0  &0  &0  &T_o  &0  &T_o  &T_{i,sp} 
# \end{matrix}]^T$$
# 
# i.e. $b_0 = b_8 = b_{10} = T_o$ and $b_{11} = T_{i,sp}$ where:
# - $T_o$ is the outdoor temperature, °C;
# - $T_{i,sp}$ - set-point temperaure for the indoor air, °C.
# 
# Since the temperature sorces $T_o$ and $T_{i,sp}$ are [time series](https://en.m.wikipedia.org/wiki/Time_series), in vector $b$ the branches which contain temperature sources are designated by $1$ and the branches without any temeprature source by $0$.

# In[23]:


b = np.zeros(12)        # branches
b[[0, 8, 10, 11]] = 1   # branches with temperature sources
print(f'b = ', b)


θ = ['θ0', 'θ1', 'θ2', 'θ3', 'θ4', 'θ5', 'θ6', 'θ7']
q = ['q0', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q7', 'q8', 'q9', 'q10', 'q11']


b = pd.Series(['To', 0, 0, 0, 0, 0, 0, 0, 'To', 0, 'To', 'Ti_sp'],
              index=q)

# In[24]:


f = pd.Series(['Φo', 0, 0, 0, 'Φi', 0, 'Qa', 'Φa'], index=θ)


# ### y: output vector
# 
# The vector of outputs is $y$, of size $n_{\theta}$, the number of nodes (in this example 8). The non-zero values of $y$ indicate the nodes which are the outputs of the model.
# 
# For the thermal circuit shown in Figure 3, if the output is the indoor air temperature, then the output vector is:
# 
# $$y = [\begin{matrix}
# 0  &0  &0  &0  &0  &0  &\theta_6 &0 
# \end{matrix}]^T$$
# 
# In vector $y$, the nodes for which the temperatures are outputs are noted by $1$ and the other nodes by $0$.

# In[25]:


y = np.zeros(8)         # nodes
y[[6]] = 1              # nodes (temperatures) of interest
print(f'y = ', y)

"""
TC
nodes and branches
"""
A = pd.DataFrame(A, index=q, columns=θ)
G = pd.Series(G, index=q)
C = pd.Series(C, index=θ)
b = pd.Series(b, index=q)
f = pd.Series(f, index=θ)
y = pd.Series(y, index=θ)

TC = {"A": A,
      "G": G,
      "C": C,
      "b": b,
      "f": f,
      "y": y}


# In[26]:

"""
"""
[As, Bs, Cs, Ds, us] = dm4bem.tc2ss(TC)


bss = np.zeros(12)        # temperature sources
bss[[0, 8, 10]] = 10      # outdoor temperature
bss[[11]] = 20            # indoor set-point temperature

fss = np.zeros(8)         # flow-rate sources

diag_G = pd.DataFrame(np.diag(G), index=G.index, columns=G.index)

diag_G = pd.DataFrame(np.diag(G), index=G.index, columns=G.index)

θss = np.linalg.inv(A.T @ diag_G @ A) @ (A.T @ diag_G @ bss + fss)
print(f'θss = {θss} °C')


bT = np.array([10, 10, 10, 20])     # [To, To, To, Tisp]
fQ = np.array([0, 0, 0, 0])         # [Φo, Φi, Qa, Φa]
uss = np.hstack([bT, fQ])
print(f'uss = {uss}')

inv_As = pd.DataFrame(np.linalg.inv(As),
                      columns=As.index, index=As.index)
yss = (-Cs @ inv_As @ Bs + Ds) @ uss

yss = float(yss.values)
print(f'yss = {yss:.2f} °C')

print(f'Max error between DAE and state-space: {abs(θss[6] - yss):.2e} °C')

# Eigenvalue analysis
λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
λ = np.sort(λ)

print('Time constants:')
print([f'{T:.2f} s' for T in -1 / λ])

dt_max = 2 * min(-1. / λ)
print(f'\nMaximum time step: {dt_max:.2f} s = {dt_max / 60:.2f} min')

t_settle = 4 * max(-1. / λ)
print(f'Minimum settling time: \
{t_settle:.0f} s = \
{t_settle / 60:.1f} min = \
{t_settle / 3600:.2f} h = \
{t_settle / (3600 * 24):.2f} days')

# time step
dt = np.floor(dt_max / 60) * 60   # s
print(f'Time step Δt = {dt} s = {dt / 60:.0f} min')

# duration: next multiple of 3600 s that is larger than t_settle
duration = np.ceil(t_settle / 3600) * 3600
print(f'Duration: \
{duration:.0f} s = \
{duration / 60:.1f} min = \
{duration / 3600:.2f} h = \
{duration / (3600 * 24):.2f} days')


# Create input_data_set
# ---------------------
# time vector
n = int(np.floor(duration / dt))    # number of time steps

# Create a DateTimeIndex starting at "00:00:00" with a time step of dt
time = pd.date_range(start="2000-01-01 00:00:00",
                           periods=n, freq=f"{int(dt)}S")

To = 10 * np.ones(n)
Ti_sp = 20 * np.ones(n)
Φa = 0 * np.ones(n)
Qa = Φo = Φi = Φa

data = {'To': To, 'Ti_sp': Ti_sp, 'Qa': Qa, 'Φo': Φo, 'Φi': Φi, 'Φa': Φa}
input_data_set = pd.DataFrame(data, index=time)

# Get input from input_data_set
u = dm4bem.inputs_in_time(us, input_data_set)


# Initial conditions
θ0 = 0                      # initial temperatures
θ_exp = pd.DataFrame(index=u.index)
θ_exp[As.columns] = θ0      # Fill θ with initial valeus θ0
θ_imp = θ_exp


I = np.eye(As.shape[0])     # identity matrix

for k in range(n - 1):
    θ_exp.iloc[k + 1] = (I + dt * As)\
        @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k]
    θ_imp.iloc[k + 1] = np.linalg.inv(I - dt * As)\
        @ (θ_imp.iloc[k] + dt * Bs @ u.iloc[k])


# outputs
y_exp = (Cs @ θ_exp.T + Ds @  u.T).T
y_imp = (Cs @ θ_imp.T + Ds @  u.T).T

# plot results
y = pd.concat([y_exp, y_imp], axis=1, keys=['Explicit', 'Implicit'])
# Flatten the two-level column labels into a single level
y.columns = y.columns.get_level_values(0)
y.plot()

print('Steady-state indoor temperature obtained with:')
print(f'- DAE model: {float(θss[6]):.4f} °C')
print(f'- state-space model: {float(yss):.4f} °C')
print(f'- steady-state response to step input: \
      {y_exp["θ6"].tail(1).values[0]:.4f} °C')

# Wetaher data
start_date = '02-01 12:00:00'
end_date = '02-07 18:00:00'

start_date = '2000-' + start_date
end_date = '2000-' + end_date

filename = '../weather_data/FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data

weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[start_date:end_date]

weather = weather.resample(str(dt) + 'S').interpolate(method='linear')


# Temperature sources
To = weather['temp_air']
Ti_sp = 20 * np.ones(weather.shape[0])


# total solar irradiance
wall_out = pd.read_csv('pd/bldg/walls_out.csv')
w0 = wall_out[wall_out['ID'] == 'w0']

surface_orientation = {'slope': w0['β'].values[0],
                       'azimuth': w0['γ'].values[0],
                       'latitude': 45}

rad_surf = dm4bem.sol_rad_tilt_surf(
    weather, surface_orientation, w0['albedo'].values[0])

Etot = rad_surf.sum(axis=1)


# Window glass properties
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass
S_g = 9         # m2, surface area of glass

# Flow-rate sources
# solar radiation
Φo = w0['α1'].values[0] * w0['Area'].values[0] * Etot
Φi = τ_gSW * w0['α0'].values[0] * S_g * Etot
Φa = α_gSW * S_g * Etot

# auxiliary (internal) sources
Qa = 0 * np.ones(weather.shape[0])

# Input data set
input_data_set = pd.DataFrame({'To': To, 'Ti_sp': Ti_sp,
                               'Qa': Qa, 'Φo': Φo, 'Φi': Φi, 'Φa': Φa})
# Get input from input_data_set
u = dm4bem.inputs_in_time(us, input_data_set)

# Initial conditions
θ0 = 20                      # initial temperatures
θ_exp = pd.DataFrame(index=u.index)
θ_exp[As.columns] = θ0      # Fill θ with initial valeus θ0

I = np.eye(As.shape[0])     # identity matrix

for k in range(u.shape[0] - 1):
    θ_exp.iloc[k + 1] = (I + dt * As)\
        @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k]


# outputs
y_exp = (Cs @ θ_exp.T + Ds @  u.T).T

Kp = TC['G']['q11']     # controller gain
q_HVAC = Kp * (u['q11'] - y_exp['θ6'])


data = pd.DataFrame({'To': To, 'θi': y_exp['θ6'],
                     'Etot': Etot, 'q_HVAC': q_HVAC})

fig, axs = plt.subplots(2, 1)
data[['To', 'θi']].plot(ax=axs[0],
                        xticks=[],
                        ylabel='Temperature, $θ$ / °C')
axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
              loc='upper right')

data[['Etot', 'q_HVAC']].plot(ax=axs[1],
                              ylabel='Heat rate, $q$ / W')
axs[1].set(xlabel='Time')
axs[1].legend(['$Φ_{total}$', '$q_{HVAC}$'],
              loc='upper right')
plt.show()


t = dt * np.arange(data.shape[0])   # time vector

fig, axs = plt.subplots(2, 1)
# plot outdoor and indoor temperature
axs[0].plot(t / 3600 / 24, data['To'], label='$θ_{outdoor}$')
axs[0].plot(t / 3600 / 24, y_exp[0, :], label='$θ_{indoor}$')
axs[0].set(ylabel='Temperatures, $θ$ / °C',
           title='Simulation for weather')
axs[0].legend(loc='upper right')

# plot total solar radiation and HVAC heat flow
axs[1].plot(t / 3600 / 24, data['Φtot'], label='$Φ_{total}$')
axs[1].plot(t / 3600 / 24, q_HVAC, label='$q_{HVAC}$')
axs[1].set(xlabel='Time, $t$ / day',
           ylabel='Heat flows, $q$ / W')
axs[1].legend(loc='upper right')

fig.tight_layout()


"""
Commented 03CubicBuilding
"""
# θs = ['θ1', 'θ3', 'θ6', 'θ7']       # state temperature nodes
# uT = ['q0', 'q8', 'q10', 'q11']     # temperature sources
# uQ = ['θ0', 'θ4', 'θ6', 'θ7']       # flow sources
# u = uT + uQ                         # inputs
# y = ['θ6']                          # output


# # In[ ]:


# pd.DataFrame(As, index=θs, columns=θs)


# # In[ ]:


# pd.DataFrame(Bs, index=θs, columns=u)


# # In[ ]:


# pd.DataFrame(Cs, index=y, columns=θs)


# # In[ ]:


# pd.DataFrame(Ds, index=y, columns=u)


# # ## Steady-state
# # [Steady-state](https://en.m.wikipedia.org/wiki/Steady_state) means that the term $C \dot \theta = 0$ in the system of DAE.
# # 
# # In [steady-state](https://en.m.wikipedia.org/wiki/Steady_state), the model can be checked if it is incorrect. Let's consider that:
# # - the controller is not active, $K_p \rightarrow 0$,
# # - the outdoor temperature is $T_o = 10 \, \mathrm{^\circ C}$,
# # - the indoor temperature setpoint is $T_{i,sp} = 20 \, \mathrm{^\circ C}$,
# # - all flow rate sources are zero.

# # In[ ]:


# b = np.zeros(12)        # temperature sources
# b[[0, 8, 10]] = 10      # outdoor temperature
# b[[11]] = 20            # indoor set-point temperature

# f = np.zeros(8)         # flow-rate sources


# # *Note*: Steady-state analysis is a test of [falsification (refutability)](https://en.m.wikipedia.org/wiki/Falsifiability) of the model, not a [verification and validation](https://en.m.wikipedia.org/wiki/Verification_and_validation). If the model does not pass the steady-state test, it means that it is wrong. If the model passes the steady-state test, it does not mean that it is correct. For example, the values of the capacities in matrix $C$ or of the conductances in matrix $G$ can be wrong even when the steady-state test is passed. 

# # ### Steady-state from differential algebraic equations (DAE)
# # The value of temperature in [steady-state](https://en.m.wikipedia.org/wiki/Steady_state) is obtained from the system of DAE by considering that $C \dot{\theta} = 0$:
# # 
# # $$\theta_{ss} = (A^T G A)^{-1}(A^T G b + f)$$
# # 
# # For the conditions mentioned above, in steady-state, all temperatures $\theta_0 ... \theta_7$, including the indoor air temperature $\theta_6$, are equal to $T_o = 10 \, \mathrm{^\circ C}$.

# # In[ ]:


# θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
# print(f'θ = {θ} °C')


# # ### Steady-state from state-space representation
# # The input vector $u$ is obtained by stacking the vectors $b_T$ and $f_Q$:
# # 
# # $$u = \begin{bmatrix} b_T \\ f_Q\end{bmatrix}$$
# # 
# # where:
# # - $b_T$ is a vector of the nonzero elements of vector $b$ of temperature sources. For the circuit presented in Figure 3, $b_T = [T_o, T_o, T_o, T_{i,sp}]^T$ corresponding to branches 0, 8, 10 and 11, where:
# #     - $T_o$ - outdoor temperature, °C;
# #     - $T_{i,sp}$ - set-point temperaure for the indoor air, °C.
# # - $f_Q$ - vector the nonzero elements of vector $f$ of flow sources. For the circuit presented in Figure 3, $f_Q = [\Phi_o, \Phi_i, \dot{Q}_a, \Phi_a]^T$ corresponding to nodes 0, 4, 6, and 7, where:
# #     - $\Phi_o$ - solar radiation absorbed by the outdoor surface of the wall, W;
# #     - $\Phi_i$ - solar radiation absorbed by the indoor surface of the wall, W;
# #     - $\dot{Q}_a$ - auxiliary heat gains (i.e., occupants, electrical devices, etc.), W;
# #     - $\Phi_a$ - solar radiation absorbed by the glass, W.
# # 
# # *Note*: Zero in vectors $b$ and $f$ indicates that there is no source on the branch or in the node, respectively. However, a source can have the value zero.

# # In[ ]:


# bT = np.array([10, 10, 10, 20])     # [To, To, To, Tisp]
# fQ = np.array([0, 0, 0, 0])         # [Φo, Φi, Qa, Φa]
# u = np.hstack([bT, fQ])
# print(f'u = {u}')


# # The steady-state value of the output of the state-space representation is obtained when $\dot \theta_{C} = 0$:
# # 
# # $$y_{ss} = (-C_s A_s^{-1} B_s + D_s) u$$

# # In[ ]:


# yss = (-Cs @ np.linalg.inv(As) @ Bs + Ds) @ u
# print(f'yss = {yss} °C')


# # The error between the steady-state values obtained from the system of DAE, $\theta_6$, and the output of the state-space representation, $y_{ss}$, 
# # 
# # $$\varepsilon = \left | \theta_6 - y_{ss} \right |$$
# # 
# # is practically zero; the slight difference is due to [numerical errors](https://en.m.wikipedia.org/wiki/Numerical_error).

# # In[ ]:


# print(f'Max error between DAE and state-space: {max(abs(θ[6] - yss)):.2e} °C')


# # ## Dynamic simulation

# # ### Time step
# # 
# # The condition for [numerical stability](https://en.m.wikipedia.org/wiki/Euler_method#Numerical_stability) of [Euler explicit integration](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Euler_method) method is
# # 
# # $$\left |  \lambda_i \Delta t + 1 \right | < 1, \forall \lambda_i, $$
# # 
# # i.e. in the complex plane, $\lambda_i \Delta t$ is inside a circle of radius 1 centered in {-1, 0j}, where:
# # - $\lambda_i$ are the eigenvalues of matrix $A_s$,
# # - $\Delta t$ - time step.
# # 
# # For positive real eigenvalues $\left \{ \lambda \in \Re |\lambda >0  \right \}$, which is the case of thermal networks, the above condition [becomes](http://www.math.iit.edu/~fass/478578_Chapter_4.pdf)
# # 
# # $$- \lambda_i \Delta t - 1  < 1, \forall \lambda_i, $$
# # 
# # or
# # 
# # $$ 0 < \Delta t < -\frac{2}{\min \lambda_i} = 2 \min -\frac{1}{\lambda_i} = 2 \min T_i$$
# # 
# # where $T_i$ are the [time constants](https://en.m.wikipedia.org/wiki/Time_constant), $T_i = - \frac{1}{\lambda_i} $

# # In[ ]:


# λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
# λ = np.sort(λ)


# # In[ ]:


# print('Time constants:') 
# print([f'{T:.2f} s' for T in -1 / λ])

# print('\n2 x Time constants:') 
# print([f'{T:.2f} s' for T in -2 / λ])

# dtmax = 2 * min(-1. / λ)
# print(f'\nMaximum time step: {dtmax:.2f} s = {dtmax / 60:.2f} min')


# # Let's chose a time step smaller than $\Delta t_{max} = \min (-2 / \lambda_i) $.

# # In[ ]:


# # time step
# dt = np.floor(dtmax / 60) * 60   # s
# print(f'dt = {dt} s = {dt / 60:.0f} min')


# # ### Settling time
# # The [settling time](https://en.m.wikipedia.org/wiki/Step_response) is roughly 4 times the larger time constant.

# # In[ ]:


# # settling time
# time_const = np.array([int(x) for x in sorted(-1 / λ)])
# print('4 * Time constants: \n', 4 * time_const, 's \n')

# t_settle = 4 * max(-1 / λ)
# print(f'Settling time: {t_settle:.0f} s = {t_settle / 60:.1f} min = {t_settle / (3600):.2f} h = {t_settle / (3600 * 24):.2f} days')


# # ### Step response
# # Let's obtain the dynamic response of the system to a [step input](https://en.m.wikipedia.org/wiki/Step_response).
# # 
# # #### Duration
# # The duration of the simulation needs to be larger than the estimated [settling time](https://en.m.wikipedia.org/wiki/Settling_time). This requires a corresponding number of time steps in the time vector.

# # In[ ]:


# # Step response
# # -------------
# # duration: next multiple of 3600 s that is larger than t_settle
# duration = np.ceil(t_settle / 3600) * 3600
# n = int(np.floor(duration / dt))    # number of time steps
# t = np.arange(0, n * dt, dt)        # time vector for n time steps

# print(f'Duration = {duration} s')
# print(f'Number of time steps = {n}')
# # pd.DataFrame(t, columns=['time'])


# # #### Input vector
# # In dynamic simulation, the inputs are [time series](https://en.m.wikipedia.org/wiki/Time_series), e.g., the oudoor temperature will have $n$ values $T_o = [T_{o(0)}, T_{o(1)}, ..., T_{o(n-1)}]$ at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # The input vector $u$ of the state-space representation is obtained by stacking the vectors $b_T$ and $f_Q$ of the system of Differential Algebraic Equations:
# # 
# # $$u = \begin{bmatrix} b_T \\ f_Q\end{bmatrix}$$
# # 
# # where:
# # - vector $b_T$ consists of the nonzero elements of vector $b$ of temperature sources; for the circuit presented in Figure 3, 
# # 
# # $$b = [\begin{matrix}
# # T_o &0  &0  &0  &0  &0  &0  &0  &T_o  &0  &T_o  &T_{i,sp} 
# # \end{matrix}]^T$$
# # and 
# # $$b_T = [T_o, T_o, T_o, T_{i,sp}]^T$$
# # corresponding to branches 0, 8, 10 and 11; 
# # - vector $f_Q$ is the nonzero elements of vector $f$ of flow sources; for the circuit presented in Figure 3,
# # 
# # $$f = [\begin{matrix}
# # \Phi_o &0  &0  &0  &\Phi_i  &0  &\dot{Q_a} &\Phi_a 
# # \end{matrix}]^T$$
# # 
# # and
# # 
# # $$f_Q = [\Phi_o, \Phi_i, \dot{Q}_a, \Phi_a]^T$$
# # 
# # corresponding to nodes 0, 4, 6, and 7.
# # 
# # For the thermal circuit shown in Figure 3, the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the input vector, $u = [u_0, u_1, ... , u_{n-1}]^T$, is:
# # 
# # $$u = 
# # \begin{bmatrix}
# # T_o\\ 
# # T_o\\ 
# # T_o\\ 
# # T_{i,sp}\\ 
# # \Phi_o\\ 
# # \Phi_i\\ 
# # \dot{Q}_a\\ 
# # \Phi_a
# # \end{bmatrix}
# # = \begin{bmatrix}
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\\ 
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\ \\ 
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\ \\ 
# # T_{i,sp(0)} & T_{i,sp(1)}& ... & T_{i,sp(n-1)}\ \\ 
# # \Phi_{o(0)} & \Phi_{o(1)} & ... & \Phi_{o(n-1)}\\
# # \Phi_{i(0)} & \Phi_{i(1)} & ... & \Phi_{i(n-1)}\\ 
# # \dot{Q}_{a(0)} & \dot{Q}_{a(1)} & ... & \dot{Q}_{a(n-1)}\\ 
# # \Phi_{a(0)} & \Phi_{a(1)} & ... & \Phi_{a(n-1)}
# # \end{bmatrix}$$
# # 
# # where:
# # - $T_o = [T_{o(0)}, T_{o(1)}, ..., T_{o(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the oudoor temperature at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # - $T_{i, sp} = [T_{{i, sp}(0)}, T_{{i, sp}(1)}, ..., T_{{i, sp}(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the setpoint indoor temperature at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # - $\Phi_o = [\Phi_{o(0)}, \Phi_{o(1)}, ..., \Phi_{o(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the solar radiation absorbed by the outdoor surface of the wall at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # - $\Phi_i = [\Phi_{i(0)}, \Phi_{i(1)}, ..., \Phi_{i(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the solar radiation absorbed by the indoor surface of the wall at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # - $\dot{Q}_a = [\dot{Q}_{a(0)}, \dot{Q}_{a(1)}, ..., \dot{Q}_{a(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the auxiliary heat gains (i.e., occupants, electrical devices, etc.) at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # - $\Phi_a = [\Phi_{a(0)}, \Phi_{a(1)}, ..., \Phi_{a(n-1)}]$ is the [time series](https://en.m.wikipedia.org/wiki/Time_series) of the solar radiation absorbed by the glass at [discrete time](https://en.m.wikipedia.org/wiki/Discrete_time_and_continuous_time#Discrete_time) $t = [t_0, t_1, ... , t_{n-1}]$.
# # 
# # Let's consider a step response in the conditions used for steady-state analysis, i.e. $T_o = 10 \, \mathrm{^\circ C}$, $T_{i,sp} = 20 \, \mathrm{^\circ C}$, and all the flow sources zero (including the HVAC system).

# # In[ ]:


# # input vector
# #u = np.zeros([8, n])                # u = [To To To Tisp Φo Φi Qa Φa]
# #u[0:3, :] = 10 * np.ones([3, n])    # To = 10 for n time steps
# #u[3, :] = 20 * np.ones([1, n])      # Tisp = 20 for n time steps

# #pd.DataFrame(u)

# # input vector
# u = np.zeros([n, 8])                # u = [To To To Tisp Φo Φi Qa Φa]
# u[:, 0:3] = 10 * np.ones([n, 3])    # To = 10 for n time steps
# u[:, 3] = 20 * np.ones([n])         # Tisp = 20 for n time steps

# pd.DataFrame(u)


# # #### Time integration

# # The state-space model
# # 
# # $$\left\{\begin{array}{rr}
# # \dot{\theta}_C=A_s \theta_C + B_s u\\ 
# # y = C_s \theta_C + D_s u
# # \end{array}\right.$$
# # 
# # is integrated in time by using [Euler forward (or explicit) method](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Euler_method) for numerical integration:
# # 
# # $$ \theta_{s,k+1} = (I + \Delta t A) \theta_{s,k} + \Delta t B u_k $$
# # 
# # and [Euler backward (or implicit) method](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Backward_Euler_method) for numerical integration:
# # 
# # $$\theta_{s,k+1} = (I - \Delta t A)^{-1} ( \theta_{s,k} + \Delta t B u_k )$$
# # 
# # where $k = 0, ... , n - 1$.

# # In[ ]:


# # initial conditions
# n_s = As.shape[0]                      # number of state variables
# θ_exp = np.zeros([n_s, t.shape[0]])    # explicit Euler in time t
# θ_imp = np.zeros([n_s, t.shape[0]])    # implicit Euler in time t

# # time integration
# I = np.eye(n_s)                        # identity matrix

# for k in range(n - 1):
#     θ_exp[:, k + 1] = (I + dt * As) @ θ_exp[:, k]        + dt * Bs @ u[k, :]
#     θ_imp[:, k + 1] = np.linalg.inv(I - dt * As) @ (θ_imp[:, k]        + dt * Bs @ u[k, :])   


# # Then, we obtain the outputs
# # 
# # $$ y = C_s \theta_s + D_s u$$
# # 
# # for explicit and for implicit Euler methods, respectively.

# # In[ ]:


# # outputs
# y_exp = Cs @ θ_exp + Ds @  u.T
# y_imp = Cs @ θ_imp + Ds @  u.T


# # The results of explicit and implicit Euler integration are practically identical.

# # In[ ]:


# fig, ax = plt.subplots()
# ax.plot(t / 3600, y_exp.T, t / 3600, y_imp.T)
# ax.set(xlabel='Time, $t$ / h',
#        ylabel='Temperatue, $θ_i$ / °C',
#        title='Step input: outdoor temperature $T_o$')
# ax.legend(['Explicit', 'Implicit'])
# ax.grid()
# plt.show()


# # > Figure 7. Step response to outdoor temperature by using Euler
# # [implicit](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Backward_Euler_method)
# # and
# # [explicit](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Euler_method) integration.
# # 
# # The value the indoor temperature obtained after the [settling time](https://en.m.wikipedia.org/wiki/Settling_time) is almost equal to the value obtained in steady-state.

# # In[ ]:


# print('Steady-state indoor temperature obtained with:')
# print(f'- DAE model: {float(θ[6]):.4f} °C')
# print(f'- state-space model: {float(yss):.4f} °C')
# print(f'- steady-state response to step input: {float(y_exp[:, -2]):.4f} °C')


# # ### Simulation with weather data

# # #### Start and end time
# # The simulation will be done from `start_date` to `end_date` indicated in the format `MM-DD HH:MM:SS` (month, day, hour:minute:second).

# # In[ ]:


# start_date = '02-01 12:00:00'
# end_date = '02-07 18:00:00'


# # The weather data are for a year. The choice of `2000` for the year is arbitrary; it used in order to respect the format `YYYY-MM-DD HH:MM:SS`.

# # In[ ]:


# start_date = '2000-' + start_date
# end_date = '2000-' + end_date
# print(f'{start_date} \tstart date')
# print(f'{end_date} \tend date')


# # #### Inputs
# # ##### Read weather data
# # Dynamic simulation needs [time series](https://en.m.wikipedia.org/wiki/Time_series) of weather data for air temperature, direct solar radiation on a normal surface and diffuse solar radiation on an horizontal surface (see the tutorial on [Weather data and solar radiation](01WeatherData.ipynb)).
# # 
# # From the weather data, we select:
# # - hourly outdoor air temperature, °C;
# # - hourly solar [direct normal irradiance](https://en.m.wikipedia.org/wiki/Direct_insolation) (or beam radiation), W/m²;
# # - hourly solar diffuse horizontal irradiance (or [diffuse sky radiation](https://en.wikipedia.org/wiki/Diffuse_sky_radiation)), W/m².

# # In[ ]:


# filename = '../weather_data/FRA_Lyon.074810_IWEC.epw'
# [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
# weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
# del data


# # The year is set to `2000` by convention and the data is selected from start to end.

# # In[ ]:


# weather.index = weather.index.map(lambda t: t.replace(year=2000))
# weather = weather.loc[start_date:end_date]


# # ##### Solar irradiance on the walls
# # For the surface orientation given by `slope`, `azimuth`and `latitude`, and the `albedo` of the surface in front of the wall, by using the weather data, we can calculate the:
# # - direct irradiance, W/m²,
# # - diffuse irradiance, W/m²,
# # - reflected irradiance, W/m²,
# # 
# # for hourly solar [irradiance](https://en.m.wikipedia.org/wiki/Solar_irradiance) on a tilted surface.

# # In[ ]:


# surface_orientation = {'slope': 90,
#                        'azimuth': 0,
#                        'latitude': 45}
# albedo = 0.2
# rad_surf = dm4bem.sol_rad_tilt_surf(
#     weather, surface_orientation, albedo)
# # pd.DataFrame(rad_surf)


# # The total solar [irradiance](https://en.m.wikipedia.org/wiki/Solar_irradiance)  $E_{tot}$, in W/m², is the sum of direct, diffuse, and reflected components.  

# # In[ ]:


# rad_surf['Φtot'] = rad_surf.sum(axis=1)


# # ##### Resample the weather data
# # The weather data is at the time-step of 1h. It needs to be resampled at time step $\Delta t$ used for numerical integration.

# # In[ ]:


# # resample weather data
# data = pd.concat([weather['temp_air'], rad_surf['Φtot']], axis=1)
# data = data.resample(str(dt) + 'S').interpolate(method='linear')
# data = data.rename(columns={'temp_air': 'To'})
# data = data.rename_axis('Time')
# # pd.DataFrame(data)


# # ##### Other inputs
# # Let's consider the indoor temperature setpoint $T_{i,sp} = 20 \, \mathrm{^\circ C}$ and the auxiliary heat flow $\dot{Q}_a = 0 \, \mathrm{W}$ constant for the whole duration of the simulation.

# # In[ ]:


# data['Ti'] = 20 * np.ones(data.shape[0])
# data['Qa'] = 0 * np.ones(data.shape[0])
# # pd.DataFrame(data)


# # ##### Input vector in time
# # The input is formed by the vectors of time series of temperature sources $\left [ T_o, T_o ,T_o, T_{i,sp} \right ]^T$ and vectors of time series of the heat flow sources $\left [ \Phi_o, \Phi_i, \dot{Q_a}, \Phi_a \right ]^T$:
# # 
# # $$u = 
# # \begin{bmatrix}
# # T_o\\ 
# # T_o\\ 
# # T_o\\ 
# # T_{i,sp}\\ 
# # \Phi_o\\ 
# # \Phi_i\\ 
# # \dot{Q}_a\\ 
# # \Phi_a
# # \end{bmatrix}
# # = \begin{bmatrix}
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\\ 
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\ \\ 
# # T_{o(0)} & T_{o(1)}& ... & T_{o(n-1)}\ \\ 
# #  T_{i,sp(0)} & T_{i,sp(1)}& ... & T_{i,sp(n-1)}\ \\ 
# # \Phi_{o,(0)} & \Phi_{o,(1)} & ... & \Phi_{o,(n-1)}\\
# # \Phi_{i,(0)} & \Phi_{i,(1)} & ... & \Phi_{i,(n-1)}\\ 
# #  \dot{Q}_{a(0)} & \dot{Q}_{a(1)} & ... & \dot{Q}_{a(n-1)}\\ 
# # \Phi_{a,(0)} & \Phi_{a,(1)} & ... & \Phi_{a,(n-1)}
# # \end{bmatrix}$$
# # 
# # where:
# # 
# # $T_o$: the time series vector of outdoor temperatures (from weather data), °C.
# # 
# # $T_{i,sp}$: time series vector of indoor setpoint temperatures, °C.
# # 
# # $\Phi_o$: time series vector of solar (i.e. short wave) radiation, in W, absorbed by the outdoor surface of the wall:
# # 
# # $$\Phi_o = \alpha_{w,SW} S_w E_{tot}$$
# # 
# # where:
# # 
# # - $\alpha_{w,SW}$ is the absortion coefficient of the outdoor surface of the wall in short wave, $0 \leqslant \alpha_{w,SW} \leqslant 1$;
# # - $S_w$ - surface area of the wall, m²;
# # - $E_{tot}$ - total solar irradiation on the wall, W/m².
# # 
# # $\Phi_i$: time series vector of short wave (i.e. solar) radiation, in W, absorbed by the indoor surfaces of the wall:
# # 
# # $$\Phi_i = \tau_{g,SW}  \alpha_{w,SW} S_g E_{tot}$$
# # 
# # where:
# # - $\tau_{g,SW}$ is the transmission coefficient of the window glass, $0 \leqslant \tau_{g,SW} \leqslant 1$;
# # - $\alpha_{w,SW}$ - absortion coefficient of the indoor surface of the wall in short wave, $0 \leqslant \alpha_{w,SW} \leqslant 1$;
# # - $S_g$ - surface area of the window glass, m²;
# # - $E_{tot}$ - total solar radiation intensity on the wall, W/m².
# # 
# # $\dot{Q}_a$: time vector of auxiliary heat flows (from occupants, electrical devices, etc.), W.
# # 
# # $\Phi_a$: time series vector of short wave (i.e. solar) radiation, in W, absorbed by the window glass:
# # 
# # $$\Phi_a = \alpha_{g,SW} S_g E_{tot}$$
# # 
# # where:
# # - $\alpha_{g,SW}$ is the absortion coefficient of the glass window in short wave, $0 \leqslant \alpha_{w,SW} \leqslant 1$;
# # - $S_g$ - surface area of the glass window, m²;
# # - $E_{tot}$ - total solar irradiation on the wall, W/m².

# # In[ ]:


# # input vector
# To = data['To']
# Ti = data['Ti']
# Φo = α_wSW * wall['Surface']['Layer_out'] * data['Φtot']
# Φi = τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φtot']
# Qa = data['Qa']
# Φa = α_gSW * wall['Surface']['Glass'] * data['Φtot']

# u = pd.concat([To, To, To, Ti, Φo, Φi, Qa, Φa], axis=1)
# u.columns.values[[4, 5, 7]] = ['Φo', 'Φi', 'Φa']
# pd.DataFrame(u)


# # #### Initial conditions
# # The initial value of the state-vector can be zero or different from zero.

# # In[ ]:


# θ_exp = 20 * np.ones([As.shape[0], u.shape[0]])


# # #### Time integration
# # [Explicit Euler](https://en.m.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Euler_method) integration in time,
# # 
# # $$ \theta_{s,k+1} = (I + \Delta t A) \theta_{s,k} + \Delta t B u_k $$
# # 
# # where $k = 0, ... , n - 1$,

# # In[ ]:


# for k in range(u.shape[0] - 1):
#     θ_exp[:, k + 1] = (I + dt * As) @ θ_exp[:, k]        + dt * Bs @ u.iloc[k, :]


# # yields the time variation of state variable $\theta$, from which we obtain the variation of the output (i.e. indoor temperature):
# # 
# # $$y = C_s \theta_s + D_s u$$
# # 
# # and the variation of the heat flow of the HVAC system:
# # 
# # $$q_{HVAC} = K_p (T_{i,sp} - \theta_i) = K_p (T_{i,sp} - y)$$
# # 
# # where $K_p$ is the gain of the P-controller and $T_{i,sp}$ is the HVAC-setpoint for the indoor temperature.

# # In[ ]:


# y_exp = Cs @ θ_exp + Ds @ u.to_numpy().T
# q_HVAC = Kp * (data['Ti'] - y_exp[0, :])


# # In[ ]:


# data['θi_exp'] = y_exp.T
# data['q_HVAC'] = q_HVAC.T


# # In[ ]:


# fig, axs = plt.subplots(2, 1)

# data[['To', 'θi_exp']].plot(ax=axs[0],
#                             xticks=[],
#                             ylabel='Temperature, $θ$ / °C')
# axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
#               loc='upper right')

# data[['Φtot', 'q_HVAC']].plot(ax=axs[1],
#                               ylabel='Heat rate, $q$ / W')
# axs[1].set(xlabel='Time')
# axs[1].legend(['$Φ_{total}$', '$q_{HVAC}$'],
#              loc='upper right')
# plt.show()


# # > Figure 6. Simulation in free-running with weather data using Euler explicit method of integration. a) Indoor and outdoor temperatures. b) Solar and HVAC heat flow rates.

# # In[ ]:


# t = dt * np.arange(data.shape[0])   # time vector

# fig, axs = plt.subplots(2, 1)
# # plot outdoor and indoor temperature
# axs[0].plot(t / 3600 / 24, data['To'], label='$θ_{outdoor}$')
# axs[0].plot(t / 3600 / 24, y_exp[0, :], label='$θ_{indoor}$')
# axs[0].set(ylabel='Temperatures, $θ$ / °C',
#            title='Simulation for weather')
# axs[0].legend(loc='upper right')

# # plot total solar radiation and HVAC heat flow
# axs[1].plot(t / 3600 / 24, data['Φtot'], label='$Φ_{total}$')
# axs[1].plot(t / 3600 / 24, q_HVAC, label='$q_{HVAC}$')
# axs[1].set(xlabel='Time, $t$ / day',
#            ylabel='Heat flows, $q$ / W')
# axs[1].legend(loc='upper right')

# fig.tight_layout()


# # > Figure 7. Simulation in free-running with weather data using Euler explicit method of integration. a) Indoor and outdoor temperatures. b) Solar and HVAC heat flow rates.

# # ## Discussion
# # 
# # Interchange the materials  of the layers of the wall. Discuss the step responses and the simuation for weather. Give arguments for the advantages and the disadvanted of indoor and outdoor insulation.
# # 
# # The time step depends on:
# # 
# # - P-controller gain `Kp`:
# #     - if $K_p \rightarrow \infty$, then the controller is perfect and the time step needs to be small;
# #     - if $K_p \rightarrow 0$, then, the controller is ineffective and the building is in free-running.
# # - Capacities considered into the model:
# #     - if the capacities of the air $C_a =$ `C['Air']` and of the glass $C_g =$ `C['Glass']` are considered, then the time step is small;
# #     - if the capacities of the air and of the glass are zero, then the time step is large (and the order of the state-space model is reduced).
# # 
# # The controller models an HVAC system able to heat (when $q_{HVAC} > 0$) and to cool (when $q_{HVAC} < 0$).

# # ## References
# # 
# # 1. C. Ghiaus (2013) Causality issue in the heat balance method for calculating the design heating and cooling loads, *Energy* 50: 292-301, https://doi.org/10.1016/j.energy.2012.10.024, open access preprint: [HAL-03605823](https://hal.archives-ouvertes.fr/hal-03605823/document)
# # 
# # 2. C. Ghiaus (2021). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198 (ref.)
# # [DOI 10.1007/978-3-030-76477-7_5](https://doi.org/10.1007/978-3-030-76477-7_5), open access preprint: [HAL 03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)
# # 
# # 3. J.A. Duffie, W. A. Beckman, N. Blair (2020) [Solar Engineering of Thermal Processes](https://www.eng.uc.edu/~beaucag/Classes/SolarPowerForAfrica/Solar%20Engineering%20of%20Thermal%20Processes,%20Photovoltaics%20and%20Wind.pdf), 5th ed. John Wiley & Sons, Inc. ISBN 9781119540281
# # 
# # 4. [Réglementation Thermique 2005. Méthode de calcul Th-CE.](https://pdfslide.fr/documents/rt2005-methode-de-calcul-th-ce.html) Annexe à l’arrêté du 19 juillet 2006
# # 
# # 5. H. Recknagel, E. Sprenger, E.-R. Schramek (2013) Génie climatique, 5e edition, Dunod, Paris. ISBN 978-2-10-070451-4
# # 
# # 6. J.R. Howell et al. (2021) Thermal Radiation Heat Transfer 7th edition, ISBN 978-0-367-34707-0, [A Catalogue of Configuration Factors](http://www.thermalradiation.net/indexCat.html)
# # 
# # 7. J. Widén, J. Munkhammar (2019) [Solar Radiation Theory](http://www.diva-portal.org/smash/get/diva2:1305017/FULLTEXT01.pdf), Uppsala University

# # 
