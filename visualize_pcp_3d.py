# fig.write_html("delta_v_isosurface.html")
import pickle
import plotly.graph_objects as go
import numpy as np

with open("pcp_data_1.pkl", "rb") as f:
    launch_dates, deltaT_days1, deltaT_days2, Delta_V_matrix = pickle.load(f)

X, Y, Z = np.meshgrid(launch_dates, deltaT_days1, deltaT_days2, indexing='ij')

fig = go.Figure(data=go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=Delta_V_matrix.flatten(),
    isomin=np.min(Delta_V_matrix)+20,
    isomax=np.max(Delta_V_matrix)-10,
    surface_count=40,
    colorscale='Jet',
    opacity=0.5,
    caps=dict(x_show=False, y_show=False, z_show=False),
    colorbar_title='Delta-V [km/s]'
))

fig.update_layout(
    scene=dict(
        xaxis_title='Departure Date (Julian Centuries)',
        yaxis_title='TOF1 [days]',
        zaxis_title='TOF2 [days]'
    ),
    title='3D Isosurface of Delta-V Field'
)

# Find the 4D index of the minimum Delta_V
min_index = np.unravel_index(np.argmin(Delta_V_matrix), Delta_V_matrix.shape)
i, j, k = min_index

# Get the corresponding coordinates for this minimum Delta_V
min_launch_date = launch_dates[i]
min_TOF1 = deltaT_days1[j]
min_TOF2 = deltaT_days2[k]
# Note: TOF3 is not shown in 3D plot but it's l

# Add a scatter3d trace for the minimum point
fig.add_trace(
    go.Scatter3d(
        x=[min_launch_date],
        y=[min_TOF1],
        z=[min_TOF2],
        mode="markers",
        marker=dict(size=8, color="red", symbol="circle"),
        name="Minimum ΔV"
    )
)

fig.show()