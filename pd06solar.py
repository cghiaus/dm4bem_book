#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:08:34 2024

@author: cghiaus
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem


def control_for(controller, period, dt=30, nonlinear_controller=True):
    # Obtain state-space representation
    # =================================
    # Disassembled thermal circuits
    folder_path = './pd/bldg'
    TCd = dm4bem.bldg2TCd(folder_path,
                          TC_auto_number=True)

    # Assembled thermal circuit
    ass_lists = pd.read_csv(folder_path + '/assembly_lists.csv')
    ass_matrix = dm4bem.assemble_lists2matrix(ass_lists)
    TC = dm4bem.assemble_TCd_matrix(TCd, ass_matrix)

    # TC['G']['c3_q0'] = 1e3  # Kp, controler gain
    # TC['C']['c2_θ0'] = 0    # indoor air heat capacity
    # TC['C']['c1_θ0'] = 0    # glass (window) heat capacit

    # State-space
    [As, Bs, Cs, Ds, us] = dm4bem.tc2ss(TC)

    # Eigenvaleus analysis
    # ====================
    λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As

    dt_max = 2 * min(-1. / λ)    # max time step for Euler explicit stability
    dt_max = dm4bem.round_time(dt_max)
    # dm4bem.print_rounded_time('dt_max', dt)

    # t_settle = 4 * max(-1. / λ)
    # duration: next multiple of 3600 s that is larger than t_settle
    # duration = np.ceil(t_settle / 3600) * 3600
    # dm4bem.print_rounded_time('duration', duration)

    # Simulation with weather data
    # ============================
    # Start / end time
    start_date = period[0]
    end_date = period[1]

    start_date = '2000-' + start_date
    end_date = '2000-' + end_date
    # print(f'{start_date} \tstart date')
    # print(f'{end_date} \tend date')

    # Weather
    filename = '../weather_data/FRA_Lyon.074810_IWEC.epw'
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data
    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather.loc[start_date:end_date]

    # Temperature sources
    To = weather['temp_air']

    Ti_day, Ti_night = 20, 16
    Ti_sp = pd.Series(20, index=To.index)
    Ti_sp = pd.Series(
        [Ti_day if 6 <= hour <= 22 else Ti_night for hour in To.index.hour],
        index=To.index)

    # Flow-rate sources
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
    Qa = pd.Series(0, index=To.index)

    # Input data set
    input_data_set = pd.DataFrame({'To': To, 'Ti_sp': Ti_sp,
                                   'Φo': Φo, 'Φi': Φi, 'Qa': Qa, 'Φa': Φa,
                                   'Etot': Etot})

    # Time integration
    # ----------------
    # Resample hourly data to time step dt
    input_data_set = input_data_set.resample(
        str(dt) + 'S').interpolate(method='linear')

    # Get input from input_data_set
    u = dm4bem.inputs_in_time(us, input_data_set)

    # initial conditions
    θ0 = 20                     # initial temperatures
    θ_exp = pd.DataFrame(index=u.index)
    θ_exp[As.columns] = θ0      # Fill θ with initial valeus θ0

    # time integration
    I = np.eye(As.shape[0])     # identity matrix

    for k in range(1, u.shape[0] - 1):
        if nonlinear_controller:
            exec(controller)

        θ_exp.iloc[k + 1] = (I + dt * As)\
            @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k]

    # outputs
    y = (Cs @ θ_exp.T + Ds @  u.T).T

    Kp = TC['G']['c3_q0']     # W/K, controller gain
    S = 3 * 3                 # m², surface area of the toy house

    # q_HVAC / [W/m²]
    if nonlinear_controller:
        q_HVAC = u['c2_θ0']
    else:
        q_HVAC = Kp * (u['c3_q0'] - y['c2_θ0']) / S  # W/m²

    # plot
    data = pd.DataFrame({'To': input_data_set['To'],
                         'θi': y['c2_θ0'],
                         'Etot': input_data_set['Etot'],
                         'q_HVAC': q_HVAC})

    fig, axs = plt.subplots(2, 1, sharex=True)
    data[['To', 'θi']].plot(ax=axs[0],
                            xticks=[],
                            ylabel='Temperature, $θ$ / [°C]')
    axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
                  loc='upper right')
    axs[0].grid(True)

    data[['Etot', 'q_HVAC']].plot(ax=axs[1],
                                  ylabel='Heat rate, $q$ / [W / m²]')
    axs[1].set(xlabel='Time')
    axs[1].legend(['$E_{total}$', '$q_{HVAC}$'],
                  loc='upper right')
    axs[1].grid(True)
    plt.show()


# in summer
period = ['07-01 12:00:00', '07-03 12:00:00']

solar = """
Tisp = 20   # indoor setpoint temperature, °C
Δθ = 1      # temperature deadband, °C
Kpp = 3e2   # controller gain

if θ_exp.iloc[k - 1]['c2_θ0'] < Tisp + Δθ:
    u.iloc[k]['c2_θ0'] = 0
else:
    u.iloc[k]['c2_θ0'] = Kpp * (Tisp - θ_exp.iloc[k - 1]['c2_θ0'])
    u.iloc[k]['c2_θ0'] = min(u.iloc[k]['c2_θ0'], 0)
    u.iloc[k]['ow0_θ4'] *= 0.1
"""

control_for(solar, period)

# Cooling
# =======
# in summer
cooling = """
Tisp = 20   # indoor setpoint temperature, °C
Δθ = 1      # temperature deadband, °C
Kpp = 3e2   # controller gain

if Tisp + Δθ > θ_exp.iloc[k - 1]['c2_θ0']:
    u.iloc[k]['c2_θ0'] = 0
else:
    u.iloc[k]['c2_θ0'] = Kpp * (Tisp - θ_exp.iloc[k - 1]['c2_θ0'])
    u.iloc[k]['c2_θ0'] = min(u.iloc[k]['c2_θ0'], 0)
"""
control_for(cooling, period)

# # Heating & Cooling
# # =======
# # in summer
# heat_cool = """
# Tisp = 20   # indoor setpoint temperature, °C
# Δθ = 1      # temperature deadband, °C
# Kpp = 3e2   # controller gain

# if Tisp < θ_exp.iloc[k - 1]['c2_θ0'] < Tisp + Δθ:
#     u.iloc[k]['c2_θ0'] = 0
# else:
#     u.iloc[k]['c2_θ0'] = Kpp * (Tisp - θ_exp.iloc[k - 1]['c2_θ0'])
# """
# period = ['07-01 12:00:00', '07-03 12:00:00']
# # control_for(heat_cool, period)
