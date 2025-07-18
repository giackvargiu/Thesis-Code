#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 17:02:56 2025

@author: giacomovargiu
"""
# fig.write_html("delta_v_isosurface.html")
# %matplotlib qt
# %matplotlib inline

import pickle
import numpy as np
from functions_plan import to_julian, state_vector
from Lambert import Lambert_3D
from ANGLE import turning_angle
from Gravity_Assist import compute_swingby_parameters


# Gravitational parameter μ (km^3/s^2)
mu_planets = [
    2.2032e4,    # Mercury
    3.24859e5,   # Venus
    3.986004418e5, # Earth
    4.282837e4,  # Mars
    1.26686534e8, # Jupiter
    3.7931187e7, # Saturn
    5.793939e6,  # Uranus
    6.836529e6   # Neptune
]

radius_planets = [
    2439.7,   # Mercury
    6051.8,   # Venus
    6371.0,   # Earth
    3389.5,   # Mars
    71492.0,  # Jupiter
    60268.0,  # Saturn
    25559.0,  # Uranus
    24764.0   # Neptune
]

soi_planets = [
    1.12e6,    # Mercury
    6.16e6,    # Venus
    9.25e5,    # Earth
    5.77e5,    # Mars
    4.82e7,    # Jupiter
    5.48e7,    # Saturn
    5.18e7,    # Uranus
    8.66e6     # Neptune
]

n = 1000
Planet_1 = 3
Planet_2 = 4
Planet_3 = 3
DeltaT_min1 = 100
DeltaT_Max1 = 900
DeltaT_min2 = 100
DeltaT_Max2 = 900

# Insert the initial (departure) date here
Y1 = 2025
M1 = 1
D1 = 1
utch1 = 0.0
utcm1 = 0.0
utcs1 = 0.0
JD1, JC1 = to_julian(Y1, M1, D1, utch1, utcm1, utcs1)

Y2 = 2029
M2 = 1
D2 = 1
utch2 = 0.0
utcm2 = 0.0
utcs2 = 0.0
JD2, JC2 = to_julian(Y2, M2, D2, utch2, utcm2, utcs2)


## 3 DIMENSIONS CREATED
launch_dates = np.linspace(JC1, JC2, n)

deltaT_days1 = np.linspace(DeltaT_min1, DeltaT_Max1, n)

deltaT_days2 = np.linspace(DeltaT_min2, DeltaT_Max2, n)

# Initialize  3D Delta-V matrix
Delta_V_matrix = np.zeros((n, n, n))

for i in range(n):
    
    for j in range(n):
        launch_date = launch_dates[i]
        delta_t_days = deltaT_days1[j]
        delta_t_centuries = delta_t_days / 36500
        delta_t_seconds = delta_t_days * 86400
        arrival_date = launch_date + delta_t_centuries


        # Planetary parameters
        R1, V_planet_1, theta_1 = state_vector(Planet_1, launch_date)
        R2, V_planet_2, theta_2 = state_vector(Planet_2, arrival_date)

        h_transfer, intersection_SOI_dep, intersection_SOI_arr, V1, V2, r_arr, r_dep, x_dep, y_dep, z_dep, x_arr, y_arr, z_arr, x_transfer_rot, y_transfer_rot, z_transfer_rot = Lambert_3D(launch_date, arrival_date, Planet_1, Planet_2, delta_t_seconds)
        

        for k in range(n):
            launch_date2 = arrival_date
            delta_t_days2 = deltaT_days2[k]
            delta_t_centuries2 = delta_t_days2 / 36500
            delta_t_seconds2 = delta_t_days2 * 86400
            arrival_date2 = arrival_date + delta_t_centuries2

            R3, V_planet_3, theta_3 = state_vector(Planet_3, arrival_date2)

            # Lambert for P2 → P3
            h_transfer2, intersection_SOI_dep2, intersection_SOI_arr2, V3, V4, r_arr2, r_dep2, x_dep2, y_dep2, z_dep2, x_arr2, y_arr2, z_arr2, x_transfer_rot2, y_transfer_rot2, z_transfer_rot2 = Lambert_3D(launch_date2, arrival_date2, Planet_2, Planet_3, delta_t_seconds2)
            
            # Compute Delta V
            Delta_V1_vec = V1 - V_planet_1
            Delta_V1 = np.linalg.norm(Delta_V1_vec)*1e-3
            Delta_V2_vec = V_planet_3 - V4
            Delta_V2 = np.linalg.norm(Delta_V2_vec)*1e-3
            
            # Gravity assist
            v_in = V2 - V_planet_2
            v_out = V3 - V_planet_2
            V1 = np.linalg.norm(v_in)*1e-3
            V2 = np.linalg.norm(v_out)*1e-3
            theta, deg = turning_angle(v_in, v_out)

            # Constants
            mu_planet = mu_planets[Planet_2]
            r_planet = radius_planets[Planet_2]  # km (Jupiter's radius)
            r_soi_planet = soi_planets[Planet_2] # km (Jupiter's SOI)

            r_periapsis, Delta_V3, r_SOI, a1, a2 = compute_swingby_parameters(V1, V2, mu_planet, theta, r_planet, r_soi_planet)
            
            Delta_V = Delta_V1 + Delta_V2 + Delta_V3
            #print(Delta_V)
            # Store results
            Delta_V_matrix[i, j, k] = Delta_V


# Save to pickle
with open("Ea_Ma_Ea_3D_PCP_n1000.pkl", "wb") as f:
    pickle.dump((launch_dates, deltaT_days1, deltaT_days2, Delta_V_matrix), f)
    

