# %matplotlib qt
# %matplotlib inline

import pickle
import numpy as np
import matplotlib.pyplot as plt

# === Step 1: Load from pickle ===
def load_data(filepath):
    with open("pcp_4D_data.pkl", "rb") as f:
        launch_dates, deltaT_days1, deltaT_days2, deltaT_days3, Delta_V_matrix = pickle.load(f)
    return launch_dates, deltaT_days1, deltaT_days2, deltaT_days3, Delta_V_matrix

# === Step 2: Project 4D matrix to 2D by minimizing over TOF2 and TOF3 ===
def project_min_and_argmin(delta_v_4d):
    shape = delta_v_4d.shape  # (i, j, k, l)
    delta_v_flat = delta_v_4d.reshape(shape[0], shape[1], -1)

    min_values = np.min(delta_v_flat, axis=2)
    min_indices = np.argmin(delta_v_flat, axis=2)

    k_indices = min_indices // shape[3]
    l_indices = min_indices % shape[3]

    return min_values, k_indices, l_indices

# === Step 3: Plot the projected PCP ===
def plot_projected_pcp(delta_v_2d, launch_dates, deltaT_days1):
    t0_grid, tof1_grid = np.meshgrid(launch_dates, deltaT_days1, indexing='ij')

    plt.figure(figsize=(10, 6))
    contour = plt.contourf(t0_grid, tof1_grid, delta_v_2d, levels=100, cmap='jet')
    plt.colorbar(contour, label='Min ΔV over TOF2 & TOF3')
    plt.xlabel('Launch Date [MJD or days]')
    plt.ylabel('TOF1 (P1 → P2) [days]')
    plt.title('Projected Pork Chop Plot (Min ΔV over TOF2 and TOF3)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# === Main ===
if __name__ == "__main__":
    # Adjust the path to your actual pickle file
    pickle_file = "pcp_4D_data.pkl"

    # Load data
    launch_dates, deltaT_days1, deltaT_days2, deltaT_days3, Delta_V_matrix = load_data(pickle_file)

    # Project and extract minimum Delta-V over TOF2 & TOF3
    delta_v_2d, k_opt, l_opt = project_min_and_argmin(Delta_V_matrix)

    # Plot
    plot_projected_pcp(delta_v_2d, launch_dates, deltaT_days1)
