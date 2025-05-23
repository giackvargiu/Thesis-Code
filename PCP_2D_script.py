#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 16:15:59 2025

@author: giacomovargiu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from functions_plan import to_julian, state_vector
from Lambert import Lambert_3D


n = 75
Planet_1 = 3
Planet_2 = 4
DeltaT_min = 50
DeltaT_Max = 700

# Insert the initial (departure) date here
Y1 = 2020
M1 = 1
D1 = 1
utch1 = 0.0
utcm1 = 0.0
utcs1 = 0.0
JD1, JC1 = to_julian(Y1, M1, D1, utch1, utcm1, utcs1)

Y2 = 2025
M2 = 1
D2 = 1
utch2 = 0.0
utcm2 = 0.0
utcs2 = 0.0
JD2, JC2 = to_julian(Y2, M2, D2, utch2, utcm2, utcs2)


launch_dates = np.linspace(JC1, JC2, n)

deltaT_days = np.linspace(DeltaT_min, DeltaT_Max, n)

# Initialize Delta-V matrix
Delta_V_matrix = np.zeros((n, n))

for i in range(n):
    
    for j in range(n):
        launch_date = launch_dates[i]
        delta_t_days = deltaT_days[j]
        delta_t_centuries = delta_t_days / 36500
        delta_t_seconds = delta_t_days * 86400
        arrival_date = launch_date + delta_t_centuries


        # Planetary parameters
        R1, V_planet_1, theta_1 = state_vector(Planet_1, launch_date)
        R2, V_planet_2, theta_2 = state_vector(Planet_2, arrival_date)

        h_transfer, intersection_SOI_dep, intersection_SOI_arr, V1, V2, r_arr, r_dep, x_dep, y_dep, z_dep, x_arr, y_arr, z_arr, x_transfer_rot, y_transfer_rot, z_transfer_rot = Lambert_3D(launch_date, arrival_date, Planet_1, Planet_2, delta_t_seconds)
        

        # Compute Delta V
        Delta_V1_vec = V1 - V_planet_1
        Delta_V1 = np.linalg.norm(Delta_V1_vec)*1e-3
        Delta_V2_vec = V_planet_2 - V2
        Delta_V2 = np.linalg.norm(Delta_V2_vec)*1e-3
        Delta_V = Delta_V1 + Delta_V2

        # Store results
        Delta_V_matrix[i, j] = Delta_V


## PLOTTING

# Create meshgrid
X, Y = np.meshgrid(deltaT_days, launch_dates)

# Clip ΔV values to max 50 km/s for visualization clarity
Delta_V_clipped = np.clip(Delta_V_matrix, 0, 50)

# Plot with contourf for smooth transitions
plt.figure(figsize=(12, 7))
cp = plt.contourf(X, Y, Delta_V_clipped, levels=100, cmap='jet')

# Color bar
cbar = plt.colorbar(cp)
cbar.set_label("Total ΔV (km/s)", fontsize=12)

# Axis labels and title
plt.xlabel("Time of Flight (days)", fontsize=12)
plt.ylabel("Launch Date", fontsize=12)
plt.title("Pork Chop Plot – ΔV vs Launch Date and Time of Flight", fontsize=14)

# Format date axis
plt.gca().yaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
plt.gcf().autofmt_xdate()

plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()