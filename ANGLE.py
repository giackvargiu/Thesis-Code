#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 15:39:42 2025

@author: giacomovargiu
"""
import numpy as np

def angle_from_coord(r1, r2):
    """
    Returns the angle from r1 to r2 (in radians), measured counterclockwise in the orbital plane.
    Result is in [0, 2π).
    """
    r1 = np.array(r1)
    r2 = np.array(r2)

    angle = np.arctan2(r2[1], r2[0]) - np.arctan2(r1[1], r1[0])
    angle = angle % (2 * np.pi)  # Ensure result is in [0, 2π)

    return angle

import numpy as np

def turning_angle(v_in, v_out):
    """
    Returns the turning angle between incoming and outgoing velocity vectors.
    Inputs:
        v_in:  Incoming velocity vector (relative to planet)
        v_out: Outgoing velocity vector (relative to planet)
    Output:
        Turning angle in degrees and radians
    """
    v_in = np.array(v_in)
    v_out = np.array(v_out)

    cos_delta = np.dot(v_in, v_out) / (np.linalg.norm(v_in) * np.linalg.norm(v_out))
    delta_rad = np.arccos(np.clip(cos_delta, -1.0, 1.0))
    delta_deg = np.degrees(delta_rad)

    return delta_rad, delta_deg