#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:04:36 2025

@author: giacomovargiu
"""

import numpy as np
import matplotlib.pyplot as plt

def rotation_matrix_from_vinfs(v_in, v_out):
    v_in = v_in / np.linalg.norm(v_in)
    cross = np.cross(v_in, v_out)
    z_hat = cross / np.linalg.norm(cross)
    y_hat = np.cross(z_hat, v_in)
    x_hat = v_in
    return np.column_stack((x_hat, y_hat, z_hat))


def generate_hyperbola_2d(a, e, theta_range):
    r = a * (e**2 - 1) / (1 + e * np.cos(theta_range))
    x = r * np.cos(theta_range)
    y = r * np.sin(theta_range)
    z = np.zeros_like(x)
    return np.vstack((x, y, z))


def plot_swingby_arc_3D(r_pi, a1, a2, v_inf_in, v_inf_out, r_planet, r_SOI):
    """
    Plot the 3D swingby arc centered on Jupiter.

    Parameters:
    - ax: matplotlib 3D axis
    - r_pi: periapsis radius (same unit as position)
    - a1: semi-major axis of inbound hyperbola
    - a2: semi-major axis of outbound hyperbola
    - v_inf_in: inbound velocity vector (3D)
    - v_inf_out: outbound velocity vector (3D)
    - r_jupiter: Jupiter's position vector (3D)
    - r_SOI: optional sphere of influence radius (for context)
    """
    e1 = 1 + r_pi / a1
    e2 = 1 + r_pi / a2

    # Define theta ranges for inbound and outbound arcs
    theta_in = np.linspace(-np.radians(89), np.radians(10), 500)
    theta_out = np.linspace(-np.radians(10), np.radians(89), 500)

    # Generate in 2D
    arc_in = generate_hyperbola_2d(a1, e1, theta_in)
    arc_out = generate_hyperbola_2d(a2, e2, theta_out)

    # Rotate into 3D
    R = rotation_matrix_from_vinfs(v_inf_in, v_inf_out)
    arc_in_3d = R @ arc_in + r_planet[:, np.newaxis]
    arc_out_3d = R @ arc_out + r_planet[:, np.newaxis]

    # Plot
    
    # Configurar la figura y el eje 3D con fondo negro
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(*arc_in_3d, color='blue', label='Inbound Trajectory')
    ax.plot3D(*arc_out_3d, color='red', label='Outbound Trajectory')
    ax.scatter(*r_planet, color='orange', s=50, label='Planet (Jupiter)')
    ax.scatter(*((R @ np.array([[r_pi], [0], [0]])).flatten() + r_planet),
               color='gold', label='Periapsis')

    if r_SOI:
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = r_SOI * np.cos(u) * np.sin(v) + r_planet[0]
        y = r_SOI * np.sin(u) * np.sin(v) + r_planet[1]
        z = r_SOI * np.cos(v) + r_planet[2]
        ax.plot_wireframe(x, y, z, color='gray', alpha=0.3, linewidth=0.5)

    # Decorations
    
    plt.axis('equal')
    plt.grid()
    plt.xlim(-r_SOI/7, r_SOI/7)
    plt.ylim(-r_SOI/7, r_SOI/7)
    plt.zlim(-r_SOI/7, r_SOI/7)
    plt.show()
    

def plot_3D_flyby_trajectory_only(r_pi, a1, a2, v_inf_in, v_inf_out, r_SOI=None):
    e1 = 1 + r_pi / a1
    e2 = 1 + r_pi / a2

    theta_in = np.linspace(-np.radians(120), np.radians(0), 500)
    theta_out = np.linspace(-np.radians(0), np.radians(120), 500)

    arc_in = generate_hyperbola_2d(a1, e1, theta_in)
    arc_out = generate_hyperbola_2d(a2, e2, theta_out)

    R = rotation_matrix_from_vinfs(v_inf_in, v_inf_out)
    arc_in_3d = R @ arc_in
    arc_out_3d = R @ arc_out
    periapsis_3d = R @ np.array([[r_pi], [0], [0]])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d', facecolor='black')

    ax.plot3D(*arc_in_3d, 'b', label='Inbound')
    ax.plot3D(*arc_out_3d, 'r', label='Outbound')
    ax.scatter(0, 0, 0, color='orange', label='Planet (centered)')
    ax.scatter(*periapsis_3d.flatten(), color='gold', label='Periapsis')
    
    
    # Personalizar el fondo y los ejes
    ax.xaxis.pane.set_facecolor((0, 0, 0, 0))
    ax.yaxis.pane.set_facecolor((0, 0, 0, 0))
    ax.zaxis.pane.set_facecolor((0, 0, 0, 0))
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.zaxis.label.set_color('white')
    ax.title.set_color('white')
    ax.tick_params(colors='white')
    ax.grid(False)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    ax.set_title("3D Gravity Assist Geometry (Planet-Centered)")

    plt.tight_layout()
    plt.show()
