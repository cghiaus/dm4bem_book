#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem


TC = dm4bem.file2TC('./toy_model/TC.csv', name='', auto_number=False)
TC['G']['q11'] = 1e3      # Kp -> ∞, almost perfect controller
# TC['G']['q11'] = 0        # Kp -> 0, no controller (free-floating

[As, Bs, Cs, Ds, us] = dm4bem.tc2ss(TC)
λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
dtmax = 2 * min(-1. / λ)    # max time step for Euler explicit stability
dt = dm4bem.round_time(dtmax)
dm4bem.print_rounded_time('dt', dt)


start_date = '02-01 12:00:00'
end_date = '02-07 18:00:00'


start_date = '2000-' + start_date
end_date = '2000-' + end_date
print(f'{start_date} \tstart date')
print(f'{end_date} \tend date')


filename = '../weather_data/FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data

weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[start_date:end_date]

weather = weather.resample(str(dt) + 'S').interpolate(method='linear')

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


# Flow-rate sources
# solar radiation
Φo = w0['α1'].values[0] * w0['Area'].values[0] * Etot

# Window glass properties
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass
S_g = 9         # m2, surface area of glass

Φi = τ_gSW * w0['α0'].values[0] * S_g * Etot

Φa = α_gSW * S_g * Etot

# auxiliary (internal) sources
Qa = 0 * np.ones(weather.shape[0])


# Input vector in time

# Input data set
input_data_set = pd.DataFrame({'To': To, 'Ti_sp': Ti_sp,
                               'Qa': Qa, 'Φo': Φo, 'Φi': Φi, 'Φa': Φa})
# Input vector in time from input_data_set
u = dm4bem.inputs_in_time(us, input_data_set)


# SIMULATION

# Initial conditions
θ0 = 20                     # °C, initial temperatures
θ_exp = pd.DataFrame(index=u.index)
θ_exp[As.columns] = θ0      # fill θ with initial valeus θ0


# #### Time integration

I = np.eye(As.shape[0])     # identity matrix

for k in range(u.shape[0] - 1):
    θ_exp.iloc[k + 1] = (I + dt * As) @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k]

# outputs
y_exp = (Cs @ θ_exp.T + Ds @  u.T).T

Kp = TC['G']['q11']     # controller gain
q_HVAC = Kp * (u['q11'] - y_exp['θ6'])


# #### Plots


data = pd.DataFrame({'To': To, 'θi': y_exp['θ6'],
                     'Etot': Etot, 'q_HVAC': q_HVAC})


# ##### Plots using Pandas


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

# ##### Plots using matplotlib


t = dt * np.arange(data.shape[0])   # time vector

fig, axs = plt.subplots(2, 1)
# plot outdoor and indoor temperature
axs[0].plot(t / 3600 / 24, data['To'], label='$θ_{outdoor}$')
axs[0].plot(t / 3600 / 24, y_exp.values, label='$θ_{indoor}$')
axs[0].set(ylabel='Temperatures, $θ$ / °C',
           title='Simulation for weather')
axs[0].legend(loc='upper right')

# plot total solar radiation and HVAC heat flow
axs[1].plot(t / 3600 / 24, data['Etot'], label='$Φ_{total}$')
axs[1].plot(t / 3600 / 24, q_HVAC, label='$q_{HVAC}$')
axs[1].set(xlabel='Time, $t$ / day',
           ylabel='Heat flows, $q$ / W')
axs[1].legend(loc='upper right')

fig.tight_layout()
