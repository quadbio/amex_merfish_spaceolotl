import glasbey
import numpy as np

import plotly.graph_objects as go
from spaceolotl.fct.load import get_data
from spaceolotl._constants import GENES_LABEL

def plot_space(input):

    if not get_data(input):
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=True,
        autosize=True,
        yaxis_scaleanchor="x",
        xaxis=dict(title="[mm]"),
        yaxis=dict(title="[mm]")
    )
    fig.update_xaxes(
        tickvals=[-1000, 0, 1000],
        ticktext=[-1, 0, 1]
    )
    fig.update_yaxes(
        tickvals=[-1000, 0, 1000],
        ticktext=[-1, 0, 1]
    )
    
    return fig

def add_space_clusters(input, widget):

    if not (adata := get_data(input)):
        return None
    
    if not input.select_resolution():
        return None

    widget.data = [trace for trace in widget.data if not trace.name.isdigit()]

    if input.switch_clusters():

        cluster_traces = adata.uns['traces'].get(input.select_resolution(), {})
        colors = glasbey.create_palette(palette_size=len(cluster_traces), lightness_bounds=(50, 100))

        for color, (cluster_id, cluster_trace) in zip(colors, cluster_traces.items()):
            widget.add_trace(
                go.Scatter(
                    x=cluster_trace['x'],
                    y=cluster_trace['y'],
                    name=cluster_id,
                    fill='toself',
                    mode='lines',
                    line=dict(width=1),
                    marker=dict(color=color)
                )
            )

def add_space_outlines(input, widget):

    adata = get_data(input)

    if adata is None:
        return None

    p = widget
    p.data = [trace for trace in p.data if trace.name != 'cells']
    
    if input.switch_outlines():
        
        cell_traces = adata.uns['traces']['cells']
        p.add_trace(
            go.Scatter(
                name='cells',
                x=cell_traces['x'],
                y=cell_traces['y'],
                mode='lines',
                line=dict(width=0.5, color='#cfcfcf')
            )
        )

def add_space_expression(input, widget):

    if not (adata := get_data(input)):
        return None

    widget.data = [trace for trace in widget.data if trace.name not in GENES_LABEL]

    if input.switch_expression() and input.select_gene():

        gene_name = input.select_gene().split('-')[1]
        gene_expression = np.array(adata[:,input.select_gene()].X.flatten())
        x_coords = adata.obs['x']
        y_coords = adata.obs['y']
        
        widget.add_trace(
            go.Scatter(
                name = gene_name,
                x = x_coords,
                y = y_coords,
                mode='markers',
                marker=dict(
                    color = gene_expression,
                    colorscale = 'Viridis',
                    showscale = True,
                    size = input.slider_dotsize_space(),
                    colorbar=dict(
                        orientation = 'h',
                        lenmode='fraction',
                        len=0.25,
                        thickness=10,
                        y = 0.05,
                        x = 0.15)
                ),
                

            )
        )