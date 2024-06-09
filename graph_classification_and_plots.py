import umap
import pandas as pd
import torch
import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set a random seed for reproducibility
np.random.seed(40)
random_state = 40 

# Paths to the Excel files
node_identities = 'graph_features.xlsx'

# Load node features for 'system' nodes
system_features_df = pd.read_excel(node_identities)

system_features = torch.tensor(system_features_df[['# of unbound states', '# of ubd-ubd transitions', 'number of unbound states for binding', 'singly bound states not for binding ', 'singly bound states for binding ', 'singly bd-singly bd transitions']].values, dtype=torch.float32)

# Perform UMAP on filtered system features
reducer = umap.UMAP(n_components=3, random_state=random_state)
system_umap_embeddings = reducer.fit_transform(system_features)

# Define custom colors for each node identity
node_identity_colors = {
    1: 'red',
    2: 'red',
    3: 'red',
    4: 'red',
    5: 'red',
    6: 'red',
    7: 'red',
    8: 'red',
    9: 'red',
    10: 'black',
    11: 'black',
    12: 'black',
    13: 'blue',
    14: 'blue',
    15: 'blue',
    16: 'blue',
    17: 'blue',
    18: 'green',
    19: 'green',
    20: 'green',
    21: 'green',
    22: 'orange',
    23: 'orange',
    24: 'orange',
}

# Assign colors to each node identity
colors = system_features_df['Node_identity'].apply(lambda x: node_identity_colors[x]).values

# Create a DataFrame for the filtered UMAP embeddings
umap_df = pd.DataFrame(data=system_umap_embeddings, columns=['UMAP component 1', 'UMAP component 2', 'UMAP component 3'])

# Add the Node_identity column to the UMAP DataFrame
umap_df['Node_identity'] = system_features_df['Node_identity'].values

fig_system_colors = px.scatter_3d(umap_df, x='UMAP component 1', y='UMAP component 2', z='UMAP component 3', color='Node_identity', labels={'color': 'Node Identity'})
fig_system_colors.update_traces(marker=dict(size=3, color=[node_identity_colors[i] for i in system_features_df['Node_identity']]))
fig_system_colors.update_traces(customdata=system_features_df['Node_identity'], hovertemplate='<b>Node Identity:</b> %{customdata}<br><b>Cluster:</b> %{color}')

# Calculate the camera eye position
elev = 40  # Elevation angle in degrees
azim = 50  # Azimuth angle in degrees

elev_rad = np.radians(elev)
azim_rad = np.radians(azim)
distance = 1.25  # This is an arbitrary distance for the camera position

x_eye = distance * np.cos(elev_rad) * np.cos(azim_rad)
y_eye = distance * np.cos(elev_rad) * np.sin(azim_rad)
z_eye = distance * np.sin(elev_rad)

camera = dict(
    eye=dict(x=x_eye, y=y_eye, z=z_eye)
)

# Update layout for font settings, figure size, and camera angle
fig_system_colors.update_layout(
    title='UMAP Clustering for Graph Features',
    scene=dict(
        xaxis_title='UMAP component 1',
        yaxis_title='UMAP component 2',
        zaxis_title='UMAP component 3',
        aspectmode='cube'
    ),
    width=600,  # Width of the figure
    height=600  # Height of the figure
)

# Save the figure as HTML and EPS
fig_system_colors.write_html('system_features_masked.html')

# Create a 3D scatter plot using Matplotlib
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(system_umap_embeddings[:, 0], system_umap_embeddings[:, 1], system_umap_embeddings[:, 2], 
                c=[node_identity_colors[i] for i in system_features_df['Node_identity']], s=2)

ax.set_xlabel('UMAP component 1', fontsize=12, fontname='Times New Roman')
ax.set_ylabel('UMAP component 2', fontsize=12, fontname='Times New Roman')
ax.set_zlabel('UMAP component 3', fontsize=12, fontname='Times New Roman')
ax.grid(True, linewidth=0.1)

# Adjust font size for ticks
ax.tick_params(axis='both', which='major', labelsize=6)
ax.tick_params(axis='both', which='major', labelsize=6)
ax.view_init(elev=30, azim=55)

# Save the figure as EPS and PDF
plt.savefig('system_features_masked.eps', format='eps', bbox_inches='tight')


