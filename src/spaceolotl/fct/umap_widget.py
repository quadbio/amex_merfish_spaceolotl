from shiny import ui
from spaceolotl.fct.load import get_data
from spaceolotl._constants import GENES_LABEL

import numpy as np
import plotly.graph_objects as go
import glasbey

def plot_umap(input):

    if not get_data(input):
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=False,
        autosize=True,
        scene=dict(
            xaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
            zaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False))
    )

    return fig

def add_umap_clusters(input, widget):
    
    if not (adata := get_data(input)):
        return None
    
    if not input.select_resolution():
        return None

    widget.data = [trace for trace in widget.data if not trace.name.isdigit()]

    if input.switch_clusters():
        
        cluster_ids = np.unique(adata.obs[input.select_resolution()])
        colors = glasbey.create_palette(palette_size=len(cluster_ids), lightness_bounds=(50, 100))

        for color, cluster_id in zip(colors, cluster_ids):

            x_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 0]
            y_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 1]
            z_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 2]

            widget.add_trace(
                go.Scatter3d(
                name = cluster_id,
                x = x_coords,
                y = y_coords,
                z = z_coords,
                mode='markers',
                marker=dict(
                    color = color,
                    size = input.slider_dotsize_umap()
                )
                )
            )

def add_umap_expression(input, widget):

    if not (adata := get_data(input)):
        return None

    widget.data = [trace for trace in widget.data if trace.name not in GENES_LABEL]

    if input.switch_expression() and input.select_gene():

        gene_name = input.select_gene().split('-')[1]
        gene_expression = np.array(adata[:,input.select_gene()].X.flatten())
        x_coords = adata.obsm['X_umap'][:, 0]
        y_coords = adata.obsm['X_umap'][:, 1]
        z_coords = adata.obsm['X_umap'][:, 2]
        
        widget.add_trace(
            go.Scatter3d(
                name = gene_name,
                x = x_coords,
                y = y_coords,
                z = z_coords,
                mode='markers',
                marker=dict(
                    color = gene_expression,
                    colorscale = 'Viridis',
                    showscale = True,
                    size = input.slider_dotsize_umap(),
                    colorbar=dict(
                        orientation = 'h',
                        lenmode='fraction',
                        len=0.25,
                        thickness=10,
                        y = 0.05,
                        x = 0.15)
                )
            )
        )