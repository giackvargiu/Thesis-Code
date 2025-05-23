#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 18:23:01 2025

@author: giacomovargiu
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_hyperbolic_trajectory(r_pi, a1, a2, r_SOI):
    """
    Plots the hyperbolic trajectory of a spacecraft given periapsis radius (r_pi) and semi-major axis (a).
    
    Parameters:
    r_pi (float): Periapsis radius (km)
    a (float): Semi-major axis of the hyperbola (km)
    """

    n = 500
# Calculate eccentricities from r_pi = a(e - 1)
    e1 = r_pi / a1 + 1
    e2 = r_pi / a2 + 1

    # Angle range near periapsis (avoid infinite tails)
    theta = np.linspace(-np.pi, np.pi, n)
    def hyperbolic_coords(a, e, theta):
        r = (a * (e**2 - 1)) / (1 + e * np.cos(theta))
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return x, y
    # r(theta) from focus (planet at origin)

    # Inbound and outbound
    x1, y1 = hyperbolic_coords(a1, e1, theta)
    x2, y2 = hyperbolic_coords(a2, e2, theta)
    y2 = -y2  # Flip outbound for symmetry

    # Plot
    plt.figure(figsize=(8, 8))
    plt.plot(x1, y1, label="Inbound Trajectory", color='blue')
    plt.plot(x2, y2, label="Outbound Trajectory", color='red')

    # Planet (focus) at origin
    plt.scatter([0], [0], color='green', marker='x', label="Planet (Focus)")

    # Periapsis should occur at (r_pi, 0)
    plt.scatter([r_pi], [0], color='orange', label="Periapsis")

    # Decorations
    plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='k', linestyle='--', linewidth=0.5)
    plt.xlabel("Distance (km)")
    plt.ylabel("Distance (km)")
    plt.title("Accurate Hyperbolic Flyby: Trajectory through Periapsis")
    plt.legend(loc='lower right')
    plt.axis('equal')
    plt.grid()
    plt.xlim(-r_SOI/7, r_SOI/7)
    plt.ylim(-r_SOI/7, r_SOI/7)
    plt.show()

# Constants
# mu_planet = 1.26686534e8  # km^3/s^2 (Jupiter)
# vinf_minus = 5.0  # km/s
# vinf_plus = 10.0  # km/s
# delta_0 = np.radians(120)  # radians

# Semi-major axes
# a1 = mu_planet / vinf_minus**2
# a2 = mu_planet / vinf_plus**2

# Periapsis from previous result
# r_pi = 785738.5738573858
# r_SOI = 100000000

# Plot it!
# plot_hyperbolic_trajectory(r_pi, a1, a2, r_SOI)

 