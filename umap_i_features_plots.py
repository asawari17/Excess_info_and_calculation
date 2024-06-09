import umap
import pandas as pd
import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set a random seed for reproducibility
np.random.seed(42)

# Path to the Excel file
gnn_data = 'i_features.xlsx'

# Load edge features
edge_features_df = pd.read_excel(gnn_data)
edge_features_df['Node_identity'] = edge_features_df['Node_identity'].astype(int)  # Ensure Node_identity is int
edge_features = edge_features_df[['attr1', 'attr2', 'attr3', 'attr4', 'Node_identity']].values

# Perform UMAP on edge features
reducer = umap.UMAP(n_components=3)  # Set n_components=3 for 3D UMAP
edge_umap_embeddings = reducer.fit_transform(edge_features[:, :-1])  # Exclude Node_identity

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
    22: 'red',
    23: 'orange',
    24: 'orange',
}

# Function to plot UMAP
def plot_umap(umap_embeddings, node_identities, colors_dict, file_name):
    fig = px.scatter_3d(
        x=umap_embeddings[:, 0],
        y=umap_embeddings[:, 1],
        z=umap_embeddings[:, 2],
        color=node_identities,
        title='UMAP Clustering for performance Features',
        #labels={'color': 'Node Identity'}
    )
    fig.update_traces(marker=dict(size=1, color=[colors_dict[i] for i in node_identities]))
    fig.update_traces(customdata=node_identities, hovertemplate='<b>Node Identity:</b> %{customdata}<br><b>Cluster:</b> %{color}')
    fig.update_layout(scene=dict(
        xaxis_title='UMAP component 1',
        yaxis_title='UMAP component 2',
        zaxis_title='UMAP component 3',
        aspectmode='cube'
    ))  # Set aspect ratio to 1:1:1
    fig.write_html(file_name)


plot_umap(edge_umap_embeddings, edge_features[:, -1], node_identity_colors, 'filtered_edge_features_umap.html')

# Create a 3D scatter plot using Matplotlib
fig = plt.figure(figsize=(7.5, 7.5))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(edge_umap_embeddings[:, 0], edge_umap_embeddings[:, 1], edge_umap_embeddings[:, 2], 
                c=[node_identity_colors[i] for i in edge_features_df['Node_identity']], s=0.01)

# ax.set_title('UMAP Clustering for System Features', fontsize=6, fontname='Times New Roman')
ax.set_xlabel('UMAP component 1', fontsize=6, fontname='Times New Roman')
ax.set_ylabel('UMAP component 2', fontsize=6, fontname='Times New Roman')
ax.set_zlabel('UMAP component 3', fontsize=6, fontname='Times New Roman')
ax.grid(True, linewidth=0.1)

# Adjust font size for ticks
ax.tick_params(axis='both', which='major', labelsize=6)
ax.tick_params(axis='both', which='major', labelsize=6)
ax.view_init(elev=20, azim=120)

plt.savefig('edge_features_custom_colors_3d.eps', format='eps', bbox_inches='tight')