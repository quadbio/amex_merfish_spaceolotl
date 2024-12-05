from shiny import ui
from shinywidgets import output_widget

from spaceolotl._constants import DATA, GENES, LEIDEN_RESOLUTIONS
from spaceolotl.js._format import DROPDOWN_CONFIG
from spaceolotl import __version__

app_ui = ui.page_navbar(  

    # MERFISH
    ui.nav_panel(ui.HTML("MERFISH<br><span style='font-size: smaller;'>data v1.0.2</span>"),
        
        # Sidebar
        ui.page_sidebar(
            ui.sidebar(
                ui.input_selectize(
                    "select_dataset",
                    "Select dataset",
                    ['', *DATA],
                    selected = None,
                    options={
                        "render": DROPDOWN_CONFIG
                        }
                    ),
                ui.input_selectize(
                    "select_resolution",
                    "Cluster resolution",
                    LEIDEN_RESOLUTIONS,
                    ),
                ui.input_switch("switch_outlines", "Show cell outlines", False),
                ui.input_switch("switch_clusters", "Show clusters", True),
                ui.input_selectize(
                    "select_gene",
                    "Select gene",
                    ['', *GENES],
                    selected=None,
                    options={
                        "render": DROPDOWN_CONFIG
                        }
                    ),
                ui.input_switch("switch_expression", "Plot gene expression", False),
                ui.input_slider("slider_dotsize_umap", "Slider UMAP", 1, 20, 2),
                ui.input_slider("slider_dotsize_space", "Slider Space", 1, 20, 2),
                ui.input_selectize(
                    "select_gene_expression",
                    "Select genes",
                    ['', *GENES],
                    selected=None,
                    multiple=True,
                    options={
                        "render": DROPDOWN_CONFIG
                        }
                    ),
                ui.input_slider("slider_n_genes", "Slider nGenes", 1, 10, 3),
                ui.input_slider("slider_lfc", "Slider minLFC", 0.5, 2.5, 0.5, step = 0.1)

            ),
            ui.layout_columns(
                 ui.card(
                     ui.card_header("UMAP Projection"),
                     output_widget('plot_umap'),
                     full_screen = True),
                ui.card(
                    ui.card_header("Spatial plot"),
                    output_widget("plot_space"),
                    full_screen = True),
                ui.accordion(ui.accordion_panel('Gene expression per cluster', ui.output_plot("plot_gene_expression")),
                             ui.accordion_panel('Differential gene expression', ui.output_plot("plot_de")),
                             id = 'panel',
                             open = False
                             ),
                col_widths={"sm": (5, 7, 12)}
                )
                
            )
    ),

    # RNA
    ui.nav_panel(ui.HTML("RNA<br><span style='font-size: smaller;'>coming soon</span>")),

    # ATAC
    ui.nav_panel(ui.HTML("ATAC<br><span style='font-size: smaller;'>coming soon</span>")),

        # Main content

        

    # Other components of the header        
    ui.nav_spacer(),
    ui.nav_control(ui.input_dark_mode(id="mode", mode = 'dark')),

    title=ui.HTML(f"A Cell Atlas of the Axolotl Brain<br><span style='font-size: 12px; display: block; '>app v{__version__}</span>"),
    id="page"
    
)  