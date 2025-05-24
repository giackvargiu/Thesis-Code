#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 24 17:26:34 2025

@author: giacomovargiu
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates



def plot_PCP_2D_array(Delta_V_slice, launch_dates, deltaT_days, ax, cbar_container):
    """
    Plot a 2D Delta-V slice on a given matplotlib Axes `ax`.
    Uses `cbar_container` to manage colorbar persistently.
    """
    ax.clear()

    X, Y = np.meshgrid(deltaT_days, launch_dates)
    cp = ax.contourf(X, Y, Delta_V_slice, levels=100, cmap='jet')

    # Update colorbar: clear old and replace
    if cbar_container[0] is not None:
        cbar_container[0].remove()
    cbar = ax.figure.colorbar(cp, ax=ax)
    cbar.set_label("Total ΔV (km/s)", fontsize=12)
    cbar_container[0] = cbar  # Save reference

    ax.set_xlabel("Time of Flight (days)", fontsize=12)
    ax.set_ylabel("Launch Date", fontsize=12)
    ax.set_title("Pork Chop Plot – ΔV vs Launch Date and Time of Flight", fontsize=14)
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    ax.figure.autofmt_xdate()
    ax.grid(True, linestyle='--', alpha=0.5)
