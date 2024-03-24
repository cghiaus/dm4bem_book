#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:55:46 2024
Modified on Sun Mar  3 20:29:22 2024

@author: cghiaus

Example of modeling, input calculation and simulation

Workflow:

1. Create a working directory

2. Download dm4bem and save it the working directory

3. Create subfolder ./bldg in the working directory

4. Download TC1, TC2, TC3, wall_out, wall_types and assembly_lists files
and save them in subfolder ./bldg

4. Create subfolder weather_data

5. Download weaher file FRA_Lyon.074810_IWEC.epw
and save in subfolder weather_data

6. Run the script
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import dm4bem

# Inputs
controller = True
indoor_air_capacity = False
glass_capacity = False
insulation_width = 0.08  # m

date_start = '2000-02-01 12:00'
date_end = '2000-02-07 12:00'

# get insulation width from wall_types.csv
df = pd.read_csv("bldg/wall_types.csv")
insulation_width_file = df[df['Material'] == 'Insulation']['Width'].values[0]

# Model
# =====
# Disassembled thermal circuits
folder_bldg = 'bldg'
TCd = dm4bem.bldg2TCd(folder_bldg,
                      TC_auto_number=True)

# Assembled thermal circuit
ass_lists = pd.read_csv(folder_bldg + '/assembly_lists.csv')
ass_matrix = dm4bem.assemble_lists2matrix(ass_lists)
TC = dm4bem.assemble_TCd_matrix(TCd, ass_matrix)

"""
Modify parameters of thermal circuit
"""
if controller:
    TC['G']['c3_q0'] = 1e3  # Kp, controler gain
if not indoor_air_capacity:
    TC['C']['c2_θ0'] = 0    # indoor air heat capacity
if not glass_capacity:
    TC['C']['c1_θ0'] = 0    # glass (window) heat capacity

# insulation width
TC['G']['ow0_q3'] *= insulation_width_file / insulation_width
TC['G']['ow0_q4'] = TC['G']['ow0_q3']

# State-space
[As, Bs, Cs, Ds, us] = dm4bem.tc2ss(TC)

# Eigen-values analysis
λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
dt_max = 2 * min(-1. / λ)    # max time step for Euler explicit stability
dt = dm4bem.round_time(dt_max)

# Inputs
# ======

# Weather data
file_weather = 'weather_data/FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(file_weather, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data

# select weather data from date_start to date_end
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[date_start:date_end]

# Temperature sources
To = weather['temp_air']

Ti_day, Ti_night = 20, 16
Ti_sp = pd.Series(20, index=To.index)
Ti_sp = pd.Series(
    [Ti_day if 6 <= hour <= 22 else Ti_night for hour in To.index.hour],
    index=To.index)

# Flow-rate sources
# total solar irradiance
wall_out = pd.read_csv(folder_bldg + '/walls_out.csv')
w0 = wall_out[wall_out['ID'] == 'w0']

surface_orientation = {'slope': w0['β'].values[0],
                       'azimuth': w0['γ'].values[0],
                       'latitude': 45}

rad_surf = dm4bem.sol_rad_tilt_surf(
    weather, surface_orientation, w0['albedo'].values[0])

Etot = rad_surf.sum(axis=1)

# window glass properties
α_gSW = 0.38    # short wave absortivity: reflective blue glass
τ_gSW = 0.30    # short wave transmitance: reflective blue glass
S_g = 9         # m2, surface area of glass

# flow-rate sources:
# solar radiation
Φo = w0['α1'].values[0] * w0['Area'].values[0] * Etot
Φi = τ_gSW * w0['α0'].values[0] * S_g * Etot
Φa = α_gSW * S_g * Etot
# auxiliary (internal) sources
Qa = pd.Series(0, index=To.index)

# Input data set
input_data_set = pd.DataFrame({'To': To, 'Ti_sp': Ti_sp,
                               'Φo': Φo, 'Φi': Φi, 'Qa': Qa, 'Φa': Φa,
                               'Etot': Etot})

# Simulation
# ==========
# Resample hourly data to time step dt
input_data_set = input_data_set.resample(
    str(dt) + 'S').interpolate(method='linear')

# Get input vector in time from input_data_set
u = dm4bem.inputs_in_time(us, input_data_set)

# initial conditions
θ0 = 20                      # initial temperatures
θ_exp = pd.DataFrame(index=u.index)
θ_exp[As.columns] = θ0      # Fill θ with initial valeus θ0

# time integration
I = np.eye(As.shape[0])     # identity matrix

for k in range(u.shape[0] - 1):
    θ_exp.iloc[k + 1] = (I + dt * As)\
        @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k + 1]

# outputs
y = (Cs @ θ_exp.T + Ds @  u.T).T

Kp = TC['G']['c3_q0']     # W/K, controller gain
S = 3 * 3                 # m², surface area of the toy house
q_HVAC = Kp * (u['c3_q0'] - y['c2_θ0']) / S  # W/m²

# Plots
data = pd.DataFrame({'To': input_data_set['To'],
                     'θi': y['c2_θ0'],
                     'Etot': input_data_set['Etot'],
                     'q_HVAC': q_HVAC})

fig, axs = plt.subplots(2, 1)
data[['To', 'θi']].plot(ax=axs[0],
                        xticks=[],
                        ylabel='Temperature, $θ$ / °C')
axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
              loc='upper right')

data[['Etot', 'q_HVAC']].plot(ax=axs[1],
                              ylabel='Heat rate, $q$ / (W·m⁻²)')
axs[1].set(xlabel='Time')
axs[1].legend(['$E_{total}$', '$q_{HVAC}$'],
              loc='upper right')
axs[0].set_title(f'Time step: $dt$ = {dt:.0f} s;'
                 f' $dt_{{max}}$ = {dt_max:.0f} s')
plt.show()

# Outputs
dm4bem.print_rounded_time("Time step:", dt)
print(f"Mean outdoor temperature: {data['To'].mean():.1f} °C")
print(f"Min. indoor temperature: {data['θi'].min():.1f} °C")
print(f"Max. indoor temperature: {data['θi'].max():.1f} °C")

max_load = data['q_HVAC'].max()
max_load_index = data['q_HVAC'].idxmax()
Q_heat = q_HVAC[q_HVAC > 0].sum() * dt / 3.6e6      # kWh
Q_cool = q_HVAC[q_HVAC < 0].sum() * dt / 3.6e6      # kWh

print(f"Max. load: {max_load:.1f} W at {max_load_index}")
print(f"Energy consumption for heating: {Q_heat:.1f} kWh")
print(f"Energy consumption for cooling: {Q_cool:.1f} kWh")
