import pickle
import numpy as np
import plotly.graph_objects as go

# fig.write_html("delta_v_isosurface.html")

# Load the data
with open("pcp_4D_data.pkl", "rb") as f:
    launch_dates, deltaT_days1, deltaT_days2, deltaT_days3, Delta_V_matrix = pickle.load(f)

# Meshgrid for 3D coordinates (same for all TOF3 slices)
L, T1, T2 = np.meshgrid(launch_dates, deltaT_days1, deltaT_days2, indexing="ij")
x = L.flatten()
y = T1.flatten()
z = T2.flatten()

# Global min/max for color scaling
global_min = Delta_V_matrix.min()
global_max = Delta_V_matrix.max()

# Initial TOF3 index
initial_idx = 0

# Create the initial isosurface plot
fig = go.Figure(data=[
    go.Isosurface(
        x=x,
        y=y,
        z=z,
        value=Delta_V_matrix[:, :, :, initial_idx].flatten(),
        isomin=global_min+15,
        isomax=global_max,
        surface_count=15,
        colorscale="Jet",
        opacity = 0.5,
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorbar=dict(title="ΔV")
    )
])

# Update layout and slider (no animation)
steps = []
for i, tof3_val in enumerate(deltaT_days3):
    step = dict(
        method="restyle",
        args=["value", [Delta_V_matrix[:, :, :, i].flatten()]],  # Update the 'value' array of the isosurface
        label=f"{tof3_val:.1f} d"
    )
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue={"prefix": "TOF3: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    title="4D ΔV Isosurface Plot (Launch Date × TOF1 × TOF2) with TOF3 Slider",
    scene=dict(
        xaxis_title="Launch Date [days]",
        yaxis_title="TOF1 [days]",
        zaxis_title="TOF2 [days]"
    ),
    sliders=sliders
)

fig.show()
