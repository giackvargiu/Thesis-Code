#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:56:07 2025

@author: giacomovargiu
"""

import numpy as np

def compute_orbital_elements(r, v, mu):
    """
    Computes the classical orbital elements from position and velocity vectors.
    
    Parameters:
        r (array-like): Position vector [x, y, z] in meters or AU.
        v (array-like): Velocity vector [vx, vy, vz] in m/s or AU/day.
        mu (float): Standard gravitational parameter of the central body (m³/s² or AU³/day²).

    Returns:
        a (float): Semi-major axis (m or AU)
        e (float): Eccentricity
        i (float): Inclination (degrees)
        Omega (float): Longitude of ascending node (degrees)
        omega (float): Argument of periapsis (degrees)
        theta (float): True anomaly (degrees)
    """
    # Convert to numpy arrays
    r = np.array(r, dtype=float)
    v = np.array(v, dtype=float)

    # Compute magnitudes
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    # Compute specific angular momentum vector h = r × v
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    # Compute eccentricity vector e = (v × h)/mu - r/|r|
    e_vec = np.cross(v, h) / mu - r / r_mag
    e = np.linalg.norm(e_vec)  # Eccentricity magnitude

    # Compute semi-major axis a = 1 / (2/|r| - |v|^2/mu)
    energy = (v_mag**2) / 2 - mu / r_mag
    a = -mu / (2 * energy)  # Semi-major axis

    # Compute inclination i = arccos(h_z / |h|)
    i = np.degrees(np.arccos(h[2] / h_mag))

    # Compute node vector n = [0, 0, 1] × h
    n = np.cross([1, 1, 0], h)
    n_mag = np.linalg.norm(n)

    # Compute longitude of ascending node Omega
    Omega = np.degrees(np.arctan2(n[1], n[0]))
    if Omega < 0:
        Omega += 360  # Ensure Omega is between 0 and 360 degrees

    # Compute argument of periapsis omega
    omega = np.degrees(np.arccos(np.dot(n, e_vec) / (n_mag * e)))
    if e_vec[2] < 0:
        omega = 360 - omega  # Ensure correct quadrant

    # Compute true anomaly theta
    theta = np.degrees(np.arccos(np.dot(e_vec, r) / (e * r_mag)))
    if np.dot(r, v) < 0:
        theta = 360 - theta  # Ensure correct quadrant

    return a, e, i, Omega, omega, theta