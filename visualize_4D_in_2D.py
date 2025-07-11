#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive 2D Pork Chop Plot visualization from 4D Delta-V data in a pickle file.
"""

import pickle
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from Plot_4D_to_2D import plot_PCP_2D_from_4D

# Load the data
with open("pcp_4D_data.pkl", "rb") as f:
    launch_dates, deltaT_days1, deltaT_days2, deltaT_days3, Delta_V_matrix = pickle.load(f)

# Initial indices
init_k = 0
init_l = 0
n_k = Delta_V_matrix.shape[2]
n_l = Delta_V_matrix.shape[3]

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(left=0.15, right=0.9, bottom=0.35, top=0.9)  # more top margin

cbar_container = [None]

# Initial plot
plot_PCP_2D_from_4D(Delta_V_matrix, launch_dates, deltaT_days1, init_k, init_l, ax, cbar_container)

# Slider for TOF2 index (k)
ax_slider_k = plt.axes([0.25, 0.08, 0.5, 0.03], facecolor='lightgrey')
slider_k = Slider(ax_slider_k, "TOF2 index (k)", 0, n_k - 1, valinit=init_k, valstep=1)

# Slider for TOF3 index (l)
ax_slider_l = plt.axes([0.25, 0.03, 0.5, 0.03], facecolor='lightgrey')
slider_l = Slider(ax_slider_l, "TOF3 index (l)", 0, n_l - 1, valinit=init_l, valstep=1)

# Update function for both sliders
def update(val):
    k = int(slider_k.val)
    l = int(slider_l.val)
    plot_PCP_2D_from_4D(Delta_V_matrix, launch_dates, deltaT_days1, k, l, ax, cbar_container)
    fig.canvas.draw_idle()

# Connect sliders to update function
slider_k.on_changed(update)
slider_l.on_changed(update)

plt.show()
