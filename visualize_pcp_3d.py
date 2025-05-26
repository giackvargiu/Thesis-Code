# fig.write_html("delta_v_isosurface.html")
import pickle
import plotly.graph_objects as go
import numpy as np

with open("pcp_data.pkl", "rb") as f:
    launch_dates, deltaT_days1, deltaT_days2, Delta_V_matrix = pickle.load(f)

#Delta_V_matrix = np.nan_to_num(Delta_V_matrix, nan=80)
Delta_V_matrix = np.clip(Delta_V_matrix, 10, 110)

X, Y, Z = np.meshgrid(launch_dates, deltaT_days1, deltaT_days2, indexing='ij')

fig = go.Figure(data=go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=Delta_V_matrix.flatten(),
    isomin=np.min(Delta_V_matrix)+20,
    isomax=np.max(Delta_V_matrix),
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

fig.show()
