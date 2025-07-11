import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.interpolate import griddata

# Load .npz data
data = np.load("pcp_data_1.npz")
launch_dates = data['launch_dates']
deltaT_days1 = data['deltaT_days1']
deltaT_days2 = data['deltaT_days2']
Delta_V_matrix = data['Delta_V_matrix']

# Build the grid of mission parameters
L, T1, T2 = np.meshgrid(launch_dates, deltaT_days1, deltaT_days2, indexing='ij')

# Flatten everything
L_flat = L.ravel()
T1_flat = T1.ravel()
T2_flat = T2.ravel()
DV_flat = Delta_V_matrix.ravel()

# Clean Delta-V (remove NaNs or invalids)
DV_flat = np.nan_to_num(DV_flat, nan=np.nanmax(DV_flat))

# Combine into a DataFrame
df = pd.DataFrame({
    'Launch': L_flat,
    'TOF1': T1_flat,
    'TOF2': T2_flat,
    'Delta_V': DV_flat
})

# Standardize input variables
X = df[['Launch', 'TOF1', 'TOF2']]
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Apply PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
df['PC1'] = X_pca[:, 0]
df['PC2'] = X_pca[:, 1]

# Create 2D grid for contour plot
grid_x, grid_y = np.meshgrid(
    np.linspace(df['PC1'].min(), df['PC1'].max(), 200),
    np.linspace(df['PC2'].min(), df['PC2'].max(), 200)
)

# Interpolate Delta-V values
grid_z = griddata(
    (df['PC1'], df['PC2']), df['Delta_V'],
    (grid_x, grid_y), method='linear'
)

# Plot the contourf
plt.figure(figsize=(10, 6))
contour = plt.contourf(grid_x, grid_y, grid_z, levels=100, cmap='viridis')
plt.colorbar(contour, label='Delta-V (km/s)')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Delta-V Contour in PCA-Reduced Space')
plt.grid(True)
plt.tight_layout()
plt.show()

# --- PCA Metrics ---
print("\n🔍 PCA Loading Scores (Eigenvectors):")
loadings = pd.DataFrame(
    pca.components_.T,
    columns=['PC1', 'PC2'],
    index=['Launch', 'TOF1', 'TOF2']
)
print(loadings)

print("\n📊 Explained Variance:")
print(pca.explained_variance_)

print("\n📈 Explained Variance Ratio:")
print(pca.explained_variance_ratio_)

print("\n🔢 Singular Values:")
print(pca.singular_values_)

plt.figure(figsize=(10, 6))
plt.scatter(df['PC1'], df['PC2'], c=df['Delta_V'], cmap='viridis', s=5)
plt.show()