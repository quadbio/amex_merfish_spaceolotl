# Basic imports
import json
import numpy as np
import pandas as pd
import os

# Shiny related imports
from shiny import App, ui, reactive, render
from shinywidgets import output_widget, render_widget

# Plotting
import plotly.express as px
import plotly.graph_objects as go

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
plt.style.use('dark_background')
plt.rcParams.update({'font.size': 8})
import glasbey

# Data related imports
import scanpy as sc
from spaceolotl._constants import *

# General page setup
app_ui = ui.page_navbar(  

    # RNA
    ui.nav_panel("RNA", "Page A content"),

    # ATAC
    ui.nav_panel("ATAC", "Page B content"),

    # MERFISH
    ui.nav_panel("MERFISH",
        
        # Sidebar
        ui.page_sidebar(
            ui.sidebar(
                ui.input_selectize("select_dataset", "Select dataset", ['', *DATA], selected = None),
                ui.input_selectize("select_resolution", "Cluster resolution", LEIDEN_RESOLUTIONS),
                ui.input_switch("switch_outlines", "Show cell outlines", False),
                ui.input_switch("switch_clusters", "Show clusters", True),
                ui.input_selectize("select_gene", "Select gene", ['', *GENES], selected = None),
                ui.input_switch("switch_expression", "Plot gene expression", False),
                ui.input_slider("slider_dotsize_umap", "Slider UMAP", 1, 20, 2),
                ui.input_slider("slider_dotsize_space", "Slider Space", 1, 20, 2),
                ui.input_selectize("select_gene_expression", "Select genes", ['', *GENES], selected = None, multiple = True),
                ui.input_slider("slider_n_genes", "Slider nGenes", 1, 10, 3),
                ui.input_slider("slider_lfc", "Slider minLFC", 0.5, 2.5, 0.5, step = 0.1),

            ),
            ui.layout_columns(
                ui.card(
                    ui.card_header("UMAP Projection"),
                    output_widget('plot_umap'),
                    full_screen = True),
                ui.card(
                    ui.card_header("Spatial plot"),
                    output_widget("plot"),
                    full_screen = True),
                ui.accordion(ui.accordion_panel('Gene expression per cluster', ui.output_plot("gene_expression_plot")),
                             ui.accordion_panel('Differential gene expression', ui.output_plot("de_plot")),
                             id = 'panel',
                             open = False
                             ),
                col_widths={"sm": (5, 7, 12)}
                )
                
            ),
    ),

        # Main content

        

    # Other components of the header        
    ui.nav_spacer(),
    ui.nav_control(ui.input_dark_mode(id="mode", mode = 'dark')),

    title="A Cell Atlas of the Axolotl Brain",  
    id="page"
    
)  

def server(input, output, session):

    # Load the data
    @reactive.Calc
    def get_data():
        name = input.select_dataset()
        if not name:
            return None
        
        try:
            adata = sc.read_h5ad(os.path.join(DATA_DIR, name + '.h5ad'))
            with open(os.path.join(DATA_DIR, name + '_traces.json'), 'r') as f:
                traces = json.load(f)
            return {'a': adata, 't': traces}
        
        except (FileNotFoundError, IOError):
            print("File not found")
            return None

    # Set up the UMAP
    @render_widget
    def plot_umap():
        spatial_data = get_data()

        if spatial_data is None:
            return None

        fig = go.Figure()
        fig.update_layout(
            template="plotly_dark",
            showlegend=False,
            autosize=True,
            scene=dict(
        xaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
        yaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
        zaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False)
    )
        )

        return fig

    @reactive.effect
    def add_umap_clusters():
        spatial_data = get_data()

        if spatial_data is None:
            return

        p = plot_umap.widget
        p.data = [trace for trace in p.data if not trace.name.isdigit()]

        if input.switch_clusters():
            
            adata = spatial_data['a']
            cluster_ids = np.unique(adata.obs[input.select_resolution()])
            colors = glasbey.create_palette(palette_size=len(cluster_ids), lightness_bounds=(50, 100))

            for color, cluster_id in zip(colors, cluster_ids):

                x_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 0]
                y_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 1]
                z_coords = adata.obsm['X_umap'][adata.obs[input.select_resolution()] == cluster_id, 2]

                p.add_trace(
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

    @reactive.effect
    def add_umap_expression():
        spatial_data = get_data()

        if spatial_data is None:
            return

        adata = spatial_data['a']
        p = plot_umap.widget
        p.data = [trace for trace in p.data if trace.name not in GENES_LABEL]

        if input.switch_expression() and input.select_gene():

            gene_name = input.select_gene().split('-')[1]
            gene_expression = np.array(adata[:,input.select_gene()].X.flatten())
            x_coords = adata.obsm['X_umap'][:,0]
            y_coords = adata.obsm['X_umap'][:,1]
            z_coords = adata.obsm['X_umap'][:,2]
            
            p.add_trace(
                go.Scatter3d(
                    name = gene_name,
                    x = x_coords,
                    y = y_coords,
                    z = z_coords,
                    mode='markers',
                    marker=dict(
                        color = gene_expression,
                        colorscale = 'Viridis',
                        size = input.slider_dotsize_umap()
                    ),
                    

                )
            )

    # Set up the spatial plot widget
    @render_widget
    def plot():
        spatial_data = get_data()

        if spatial_data is None:
            return None

        fig = go.Figure()
        fig.update_layout(
            template="plotly_dark",
            showlegend=True,
            autosize=True,
            yaxis_scaleanchor="x"
        )

        return fig

    # Add the cluster traces
    @reactive.effect
    def add_clusters():
        spatial_data = get_data()

        if spatial_data is None:
            return

        p = plot.widget
        p.data = [trace for trace in p.data if not trace.name.isdigit()]

        if input.switch_clusters():

            cluster_traces = spatial_data['t'].get(input.select_resolution(), {})
            color = glasbey.create_palette(palette_size=len(cluster_traces), lightness_bounds=(50, 100))

            for color, (cluster_id, cluster_trace) in zip(color, cluster_traces.items()):
                p.add_trace(
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

    # Add the cell outlines
    @reactive.effect
    def add_outlines():
        spatial_data = get_data()

        if spatial_data is None:
            return

        p = plot.widget
        p.data = [trace for trace in p.data if trace.name != 'cells']
        
        if input.switch_outlines():
            
            cell_traces = spatial_data['t']['cells']
            p.add_trace(
                go.Scatter(
                    name='cells',
                    x=cell_traces['x'],
                    y=cell_traces['y'],
                    mode='lines',
                    line=dict(width=0.5, color='#cfcfcf')
                )
            )

    @reactive.effect
    def add_expression():
        spatial_data = get_data()

        if spatial_data is None:
            return

        p = plot.widget
        p.data = [trace for trace in p.data if trace.name not in GENES_LABEL]

        if input.switch_expression() and input.select_gene():

            gene_name = input.select_gene().split('-')[1]
            gene_expression = np.array(spatial_data['a'][:,input.select_gene()].X.flatten())
            x_coords = spatial_data['a'].obs['x']
            y_coords = spatial_data['a'].obs['y']
            
            p.add_trace(
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

    @render.plot()  
    def gene_expression_plot():  

        spatial_data = get_data()
        plot_genes = list(input.select_gene_expression())

        if spatial_data is None or not plot_genes:
            return

        plot_ex = sc.pl.dotplot(spatial_data['a'],
                                var_names = plot_genes,
                                swap_axes = True,
                                groupby = input.select_resolution(),
                                return_fig = True
                                )

        main_ax_ex = plot_ex.get_axes()['mainplot_ax']
        y_axis_labels = [tick.get_text() for tick in main_ax_ex.get_yticklabels()]
        short_labels_ex = [label.split("-")[-1] for label in y_axis_labels]
        main_ax_ex.set_yticklabels(short_labels_ex)

        plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95)
        plt.show()

    @render.plot()  
    def de_plot():  

        spatial_data = get_data()

        if spatial_data is None:
            return
        
        plot_de = sc.pl.rank_genes_groups_dotplot(
            spatial_data['a'],
            n_genes = input.slider_n_genes(),
            key = input.select_resolution(),
            min_logfoldchange=input.slider_lfc(),
            return_fig = True)

        main_ax_de = plot_de.get_axes()['mainplot_ax']
        x_axis_labels = [tick.get_text() for tick in main_ax_de.get_xticklabels()]
        short_labels_de = [label.split("-")[-1] for label in x_axis_labels]
        main_ax_de.set_xticklabels(short_labels_de)
        
        plt.subplots_adjust(top=1, bottom=0.2, left=0.05, right=0.95)
        plt.show()



app = App(app_ui, server)