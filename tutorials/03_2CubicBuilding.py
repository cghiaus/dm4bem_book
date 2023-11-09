#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem


# INPUTS
# ======
start_date = '02-01 12:00:00'
end_date = '02-07 18:00:00'


start_date = '2000-' + start_date
end_date = '2000-' + end_date
print(f'{start_date} \tstart date')
print(f'{end_date} \tend date')

# Input data set
# --------------
# temperature and flow-rate sources

# Read weather data
filename = '../weather_data/FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data

# set year 2000 for all years
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[start_date:end_date]


# Weather sources
To = weather['temp_air']

# total solar irradiance
wall_out = pd.read_csv('pd/bldg/walls_out.csv')
w0 = wall_out[wall_out['ID'] == 'w0']

surface_orientation = {'slope': w0['β'].values[0],
                       'azimuth': w0['γ'].values[0],
                       'latitude': 45}

rad_surf = dm4bem.sol_rad_tilt_surf(
    weather, surface_orientation, w0['albedo'].values[0])

Etot = rad_surf.sum(axis=1)

# Solar radiation absorbed by the outdoor surface of the wall
Φo = w0['α1'].values[0] * w0['Area'].values[0] * Etot

# window glass properties
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass
S_g = 9         # m2, surface area of glass

# solar radiation absorbed by the indoor surface of the wall
Φi = τ_gSW * w0['α0'].values[0] * S_g * Etot

# solar radiation absorbed by the glass
Φa = α_gSW * S_g * Etot

# Schedules
Ti_sp = pd.Series(20, index=To.index)

Ti_day, Ti_night = 20, 16

Ti_sp = pd.Series(
    [Ti_day if 6 <= hour <= 22 else Ti_night for hour in To.index.hour],
    index=To.index)

Qa = 0 * np.ones(weather.shape[0])      # auxiliary (internal) sources

Qa = pd.Series(0, index=To.index)

# Input data set
input_data_set = pd.DataFrame({'To': To, 'Ti_sp': Ti_sp,
                               'Φo': Φo, 'Φi': Φi, 'Qa': Qa, 'Φa': Φa,
                               'Etot': Etot})


input_data_set.to_csv('./toy_model/input_data_set.csv')

# Simulation
# ==========
globals().clear()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

controller = True
explicit_Euler = True
imposed_time_step = True
Δt = 600    # s, imposed time step

# MODEL
# =====
# Thermal circuit
TC = dm4bem.file2TC('./toy_model/TC.csv', name='', auto_number=False)

# by default TC['G']['q11'] = 0, i.e. Kp -> 0, no controller (free-floating)
if controller:
    TC['G']['q11'] = 1e3        # Kp -> ∞, almost perfect controller

# State-space
[As, Bs, Cs, Ds, us] = dm4bem.tc2ss(TC)

# Eigenvalues analysis
λ = np.linalg.eig(As)[0]        # eigenvalues of matrix As

# time step
if imposed_time_step:
    dt = Δt
else:
    dtmax = 2 * min(-1. / λ)    # max time step for Euler explicit stability
    dt = dm4bem.round_time(dtmax)

dm4bem.print_rounded_time('dt', dt)

# INPUT DATA SET
# ==============
input_data_set = pd.read_csv('./toy_model/input_data_set.csv',
                             index_col=0,
                             parse_dates=True)


# SIMULATION
# ==========
# Input vector in time from input_data_set
input_data_set = input_data_set.resample(
    str(dt) + 'S').interpolate(method='linear')
u = dm4bem.inputs_in_time(us, input_data_set)

# Initial conditions
θ0 = 20                     # °C, initial temperatures
θ = pd.DataFrame(index=u.index)
θ[As.columns] = θ0      # fill θ with initial valeus θ0


# #### Time integration

I = np.eye(As.shape[0])     # identity matrix

for k in range(u.shape[0] - 1):
    # explicit Euler
    θ.iloc[k + 1] = (I + dt * As) @ θ.iloc[k] + dt * Bs @ u.iloc[k]

    # implicit Euler
    θ.iloc[k + 1] = np.linalg.inv(
        I - dt * As) @ (θ.iloc[k] + dt * Bs @ u.iloc[k])

# outputs
y = (Cs @ θ.T + Ds @  u.T).T

Kp = TC['G']['q11']     # controller gain
S = 9                   # m², surface area of the toy house
q_HVAC = Kp * (u['q11'] - y['θ6']) / 9  # W/m²


# Plots


data = pd.DataFrame({'To': input_data_set['To'],
                     'θi': y['θ6'],
                     'Etot': input_data_set['Etot'],
                     'q_HVAC': q_HVAC})


# Plots using Pandas


fig, axs = plt.subplots(2, 1)
data[['To', 'θi']].plot(ax=axs[0],
                        xticks=[],
                        ylabel='Temperature, $θ$ / °C')

axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
              loc='upper right')

data[['Etot', 'q_HVAC']].plot(ax=axs[1],
                              ylabel='Heat rate, $q$ / (W/m²)')
axs[1].set(xlabel='Time')
axs[1].legend(['$E_{total}$', '$q_{HVAC}$'],
              loc='upper right')
plt.show()

# # ##### Plots using matplotlib


# t = dt * np.arange(data.shape[0])   # time vector

# fig, axs = plt.subplots(2, 1)
# # plot outdoor and indoor temperature
# axs[0].plot(t / 3600 / 24, data['To'], label='$θ_{outdoor}$')
# axs[0].plot(t / 3600 / 24, y.values, label='$θ_{indoor}$')
# axs[0].set(ylabel='Temperatures, $θ$ / °C',
#            title='Simulation for weather')
# axs[0].legend(loc='upper right')

# # plot total solar radiation and HVAC heat flow
# axs[1].plot(t / 3600 / 24, data['Etot'], label='$Φ_{total}$')
# axs[1].plot(t / 3600 / 24, q_HVAC, label='$q_{HVAC}$')
# axs[1].set(xlabel='Time, $t$ / day',
#            ylabel='Heat flows, $q$ / W')
# axs[1].legend(loc='upper right')

# fig.tight_layout()
