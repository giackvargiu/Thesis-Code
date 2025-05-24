#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 19:24:04 2025

@author: giacomovargiu
"""

import numpy as np
#import matplotlib.pyplot as plt

#from Plot_hyperbola import plot_hyperbolic_trajectory



def compute_swingby_parameters(vinf_minus, vinf_plus, mu_planet, delta_0, r_planet, r_soi_planet, tolerance=1e-1, max_iterations=100):
    
    
    # Compute semi-major axes for inbound and outbound hyperbolas
    a1 = mu_planet / vinf_minus**2
    a2 = mu_planet / vinf_plus**2

    if delta_0 > np.pi:
        print('AAAAAAAAAAAAh')
    ############################################################################################################################

    # Finding the first guesses for Newton Raphson

    # Define periapsis radius range (0 to SOI)
    r_pi = np.linspace(1e-3, r_soi_planet, 10000)

    # Compute half turning angles for inbound and outbound trajectories
    safe_ratio_1 = np.clip(a1 / (r_pi + a1), -1 + 1e-10, 1 - 1e-10)
    safe_ratio_2 = np.clip(a2 / (r_pi + a2), -1 + 1e-10, 1 - 1e-10)
    half_delta_1 = np.arcsin(safe_ratio_1)
    half_delta_2 = np.arcsin(safe_ratio_2)
    half_delta_0 = delta_0 / 2  # Target half turning angle

    # Convert angles from radians to degrees
    half_delta_1_deg = np.degrees(half_delta_1)
    half_delta_2_deg = np.degrees(half_delta_2)
    half_delta_0_deg = np.degrees(half_delta_0)

    # Plot with log10 scale on x-axis
    # plt.figure(figsize=(8, 6))
    # plt.plot(r_pi, half_delta_1_deg, label=r'$\delta_1/2$ (Inbound Half-Turn Angle)')
    # plt.plot(r_pi, half_delta_2_deg, label=r'$\delta_2/2$ (Outbound Half-Turn Angle)')
    # plt.axhline(y=half_delta_0_deg, color='r', linestyle='--', label=r'$\delta_0/2$ (Target Half-Turn Angle)')
    # plt.axvline(x=r_planet, color='g', linestyle='--', label='Jupiter Radius')
    # plt.axvline(x=r_soi_planet, color='b', linestyle='--', label='Jupiter SOI')

    # Set x-axis to log scale
    # plt.xscale('log')

    # Labels and Legend
    # plt.xlabel('Periapsis Radius (km) [log scale]')
    # plt.ylabel('Half Turning Angle (degrees)')
    # plt.title('Half Turning Angles vs. Periapsis Radius (Log Scale)')
    # plt.legend()
    # plt.grid(which='both', linestyle='--', linewidth=0.5)
    # plt.show()


    # Find the intersection points
    idx_1 = np.argmin(np.abs(half_delta_1_deg - half_delta_0_deg))  # Closest index where half_delta_1 intersects half_delta_0
    idx_2 = np.argmin(np.abs(half_delta_2_deg - half_delta_0_deg))  # Closest index where half_delta_2 intersects half_delta_0

    # Corresponding periapsis radii at the intersections
    r_intersect_1 = r_pi[idx_1]
    r_intersect_2 = r_pi[idx_2]

    # Corresponding turning angles at the intersections
    angle_intersect_1 = half_delta_1_deg[idx_1]
    angle_intersect_2 = half_delta_2_deg[idx_2]

    r_intersect_1, angle_intersect_1, r_intersect_2, angle_intersect_2

    r_periapsis = (r_intersect_1 + r_intersect_2) / 2

    #############################################################################################################################

    # Newton Raphson to find the common radius of periapsis

    iteration = 0

    while iteration < max_iterations:
        # Compute the function value
        safe1 = np.clip(a1 / (r_periapsis + a1), -1 + 1e-10, 1 - 1e-10)
        safe2 = np.clip(a2 / (r_periapsis + a2), -1 + 1e-10, 1 - 1e-10)
        f = np.arcsin(safe1) + np.arcsin(safe2) - delta_0
        
        # f = np.arcsin(a1 / (r_periapsis + a1)) + np.arcsin(a2 / (r_periapsis + a2)) - delta_0
        # Compute the derivative df/dr_periapsis
        df_dr = - (a1 / (r_periapsis + a1)) / np.sqrt(r_periapsis**2 + 2*a1 * r_periapsis) - (a2 / (r_periapsis + a2)) / np.sqrt(r_periapsis**2 + 2*a2 * r_periapsis)
        if np.abs(df_dr) < 1e-10:
            print("Dérivée trop faible. Abandon de l’itération.")
            break
        # Newton-Raphson step
        r_periapsis_new = r_periapsis - f / df_dr

        # Check for convergence
        if abs(r_periapsis_new - r_periapsis) < tolerance:
            r_periapsis = r_periapsis_new
            break

        r_periapsis = r_periapsis_new
        iteration += 1
           

    #############################################################################################################################

    # Computing the Delta_V
    if r_periapsis > r_planet:
        term1 = vinf_minus**2 - (2 * mu_planet / r_periapsis)
        term2 = vinf_plus**2 - (2 * mu_planet / r_periapsis)

        v_pi_minus = np.sqrt(np.abs(term1))
        v_pi_plus = np.sqrt(np.abs(term2))
        Delta_V = abs(v_pi_minus - v_pi_plus)

    else:
        #print("Warning: r_periapsis is zero or too small")
        Delta_V = 30
    # print(Delta_V)

    # Plot with log10 scale on x-axis
    # plt.figure(figsize=(8, 6))
    # plt.plot(r_pi, half_delta_1_deg, label=r'$\delta_1/2$ (Inbound Half-Turn Angle)')
    # plt.plot(r_pi, half_delta_2_deg, label=r'$\delta_2/2$ (Outbound Half-Turn Angle)')
    # plt.scatter(r_periapsis, half_delta_0_deg, color='purple', marker='o', label='Newton-Raphson final value of R_periapsis')
    # plt.axhline(y=half_delta_0_deg, color='r', linestyle='--', label=r'$\delta_0/2$ (Target Half-Turn Angle)')
    # plt.axvline(x=r_planet, color='g', linestyle='--', label='Planets Radius')
    # plt.axvline(x=r_soi_planet, color='b', linestyle='--', label='Planets SOI')

    # Set x-axis to log scale
    # plt.xscale('log')

    # Labels and Legend
    # plt.xlabel('Periapsis Radius (km) [log scale]')
    # plt.ylabel('Half Turning Angle (degrees)')
    # plt.title('Half Turning Angles vs. Periapsis Radius (Log Scale)')
    # plt.legend()
    # plt.grid(which='both', linestyle='--', linewidth=0.5)
    # plt.show()
    
    # plot_hyperbolic_trajectory(r_periapsis, a1, a2, r_soi_planet)
    
    return r_periapsis, Delta_V, r_soi_planet, a1, a2

# Constants
# mu_planet = 1.26686534e8  # km^3/s^2 (Jupiter)
# vinf_minus = 8.72772490103032 # km/s (Example value)
# vinf_plus = 15.5597674931533  # km/s (Example value)
# delta_0 = np.radians(81.64555145229677)  # Target turning angle (converted to radians)
# r_planet = 71492  # km (Jupiter's radius)
# r_soi_planet = 48.2e6  # km (Jupiter's SOI)

# r_periapsis, Delta_V = compute_swingby_parameters(vinf_minus, vinf_plus, mu_planet, delta_0, r_planet, r_soi_planet)
# print(f"Periapsis Radius: {r_periapsis:.2f} km")
# print(f"Required Delta V: {Delta_V:.2f} km/s")

# Constants
# mu_mars = 4.282837e4 
# V1 = 21.952631856564214  # km/s (Example value)
# V2 = 30.991395849151957  # km/s (Example value)
# turning_angle = np.radians(76.27789386404417)  # Target turning angle (converted to radians)
# r_mars = 3389.5  # km (Jupiter's radius)
# r_soi_mars = 0.577e6  # km (Jupiter's SOI)

# r_periapsis, Delta_V = compute_swingby_parameters(V1, V2, mu_mars, turning_angle, r_mars, r_soi_mars)
# print(f"Periapsis Radius: {r_periapsis:.2f} km")
# print(f"Required Delta V: {Delta_V:.2f} km/s")