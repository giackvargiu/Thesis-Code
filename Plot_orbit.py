#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 15:15:04 2025

@author: giacomovargiu
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_orbit(a, e, i, Omega, omega, r1, a1, e1, i1, Omega1, omega1, r2, a2, e2, i2, Omega2, omega2, mu=1.32712440018e20, num_points=500):
    """
    Plots a 3D orbit given its classical orbital elements.

    Parameters:
        a (float): Semi-major axis (meters or AU)
        e (float): Eccentricity (dimensionless)
        i (float): Inclination (degrees)
        Omega (float): Longitude of ascending node (degrees)
        omega (float): Argument of periapsis (degrees)
        mu (float): Standard gravitational parameter (m³/s² or AU³/day²)
        num_points (int): Number of points in the trajectory

    Returns:
        None (displays the plot)
    """
    # Convert degrees to radians
    i = np.radians(i)
    Omega = np.radians(Omega)
    omega = np.radians(omega)

    # Generate true anomaly values from 0 to 360 degrees
    theta = np.linspace(0, 2*np.pi, num_points)

    # Compute orbit in polar form
    r = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Convert to Cartesian coordinates in orbital plane
    x_orb = r * np.cos(theta)
    y_orb = r * np.sin(theta)
    z_orb = np.zeros_like(x_orb)  # The orbit starts in the orbital plane

    # Rotation matrices to transform to 3D space
    R1 = np.array([
        [np.cos(Omega), -np.sin(Omega), 0],
        [np.sin(Omega), np.cos(Omega), 0],
        [0, 0, 1]
    ])

    R2 = np.array([
        [1, 0, 0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i), np.cos(i)]
    ])

    R3 = np.array([
        [np.cos(omega), -np.sin(omega), 0],
        [np.sin(omega), np.cos(omega), 0],
        [0, 0, 1]
    ])

    # Apply rotations: 3-1-3 sequence (Ω, i, ω)
    R_total = R1 @ R2 @ R3  # Combine transformations

    # Rotate orbit points to 3D space
    orbit_3D = np.dot(R_total, np.array([x_orb, y_orb, z_orb]))

    # Extract final 3D coordinates
    x, y, z = orbit_3D[0], orbit_3D[1], orbit_3D[2]

###############################################################################

    # Convert to Cartesian coordinates in orbital plane
    r1_norm = a1 * (1 - e1**2) / (1 + e1 * np.cos(theta))
    # r1_norm = np.linalg.norm(r1)
    x1 = r1_norm * np.cos(theta)
    y1 = r1_norm * np.sin(theta)
    z1 = np.zeros_like(x1)  # The orbit starts in the orbital plane

    # Rotation matrices to transform to 3D space
    R1_1 = np.array([
        [np.cos(Omega1), -np.sin(Omega1), 0],
        [np.sin(Omega1), np.cos(Omega1), 0],
        [0, 0, 1]
    ])

    R2_1 = np.array([
        [1, 0, 0],
        [0, np.cos(i1), -np.sin(i1)],
        [0, np.sin(i1), np.cos(i1)]
    ])

    R3_1 = np.array([
        [np.cos(omega1), -np.sin(omega1), 0],
        [np.sin(omega1), np.cos(omega1), 0],
        [0, 0, 1]
    ])

    # Apply rotations: 3-1-3 sequence (Ω, i, ω)
    R_total_1 = R1_1 @ R2_1 @ R3_1  # Combine transformations

    # Rotate orbit points to 3D space
    orbit_3D_1 = np.dot(R_total_1, np.array([x1, y1, z1]))

    # Extract final 3D coordinates
    x1, y1, z1 = orbit_3D_1[0], orbit_3D_1[1], orbit_3D_1[2]    

###############################################################################

    # Convert to Cartesian coordinates in orbital plane
    r2_norm = a2 * (1 - e2**2) / (1 + e2 * np.cos(theta))
    # r2_norm = np.linalg.norm(r2)
    x2 = r2_norm * np.cos(theta)
    y2 = r2_norm * np.sin(theta)
    z2 = np.zeros_like(x2)  # The orbit starts in the orbital plane

    # Rotation matrices to transform to 3D space
    R1_2 = np.array([
        [np.cos(Omega2), -np.sin(Omega2), 0],
        [np.sin(Omega2), np.cos(Omega2), 0],
        [0, 0, 1]
    ])

    R2_2 = np.array([
        [1, 0, 0],
        [0, np.cos(i2), -np.sin(i2)],
        [0, np.sin(i2), np.cos(i2)]
    ])

    R3_2 = np.array([
        [np.cos(omega2), -np.sin(omega2), 0],
        [np.sin(omega2), np.cos(omega2), 0],
        [0, 0, 1]
    ])

    # Apply rotations: 3-1-3 sequence (Ω, i, ω)
    R_total_2 = R1_2 @ R2_2 @ R3_2  # Combine transformations

    # Rotate orbit points to 3D space
    orbit_3D_2 = np.dot(R_total_2, np.array([x2, y2, z2]))

    # Extract final 3D coordinates
    x2, y2, z2 = orbit_3D_2[0], orbit_3D_2[1], orbit_3D_2 [2]    
    
###############################################################################

    # Plot the 3D orbit
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label="Orbit", color='g')
    ax.plot(x1, y1, z1, label="Orbit", color='b')
    ax.plot(x2, y2, z2, label="Orbit", color='r')
    
    # Plot planets
    ax.scatter(r1[0], r1[1], r1[2], color='blue', marker='o', s=80, label="Planet1")
    ax.scatter(r2[0], r2[1], r2[2], color='red', marker='o', s=100, label="Planet2")
    
    # Plot central body (Sun)
    ax.scatter(0, 0, 0, color='yellow', marker='o', s=100, label="Sun")

    # Labels and legend
    ax.set_xlabel("X [AU]")
    ax.set_ylabel("Y [AU]")
    ax.set_zlabel("Z [AU]")
    ax.set_title("3D Orbital Plot")
    ax.legend()
    ax.set_box_aspect([1,1,0.02])  # Better aspect ratio

    plt.show()
