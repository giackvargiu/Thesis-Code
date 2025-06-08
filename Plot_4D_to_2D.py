#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 14:21:01 2025

@author: giacomovargiu
"""
import numpy as np
import matplotlib.dates as mdates

def plot_PCP_2D_from_4D(Delta_V_matrix, launch_dates, deltaT_days1, TOF2_index, TOF3_index, ax, cbar_container):
    """
    Plot a 2D Delta-V slice from the 4D Delta_V_matrix at specific TOF2 and TOF3 indices.
    
    Parameters:
        Delta_V_matrix: 4D numpy array (launch_dates × deltaT_days1 × deltaT_days2 × deltaT_days3)
        launch_dates: 1D array of launch dates
        deltaT_days1: 1D array of TOF1
        TOF2_index: int, index for deltaT_days2
        TOF3_index: int, index for deltaT_days3
        ax: matplotlib Axes object
        cbar_container: mutable container for colorbar
    """
    ax.clear()
    
    # Extract the 2D slice
    Delta_V_slice = Delta_V_matrix[:, :, TOF2_index, TOF3_index]
    
    # Create the meshgrid for plotting
    X, Y = np.meshgrid(deltaT_days1, launch_dates)
    
    # Plot the 2D contour plot
    cp = ax.contourf(X, Y, Delta_V_slice, levels=100, cmap='jet')
    
    # Remove previous colorbar if it exists
    if cbar_container[0] is not None:
        cbar_container[0].remove()
    cbar = ax.figure.colorbar(cp, ax=ax)
    cbar.set_label("Total ΔV (km/s)", fontsize=12)
    cbar_container[0] = cbar
    
    ax.set_xlabel("Time of Flight 1 (days)", fontsize=12)
    ax.set_ylabel("Launch Date", fontsize=12)
    ax.set_title(f"Pork Chop Plot – ΔV vs Launch Date and TOF1\n(TOF2 index: {TOF2_index}, TOF3 index: {TOF3_index})", fontsize=14)
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    ax.figure.autofmt_xdate()
    ax.grid(True, linestyle='--', alpha=0.5)
