from spaceolotl.fct.load import get_data
import matplotlib.pyplot as plt
import scanpy as sc

def plot_gene_expression(input):  

    spatial_data = get_data(input)
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

def plot_de(input):  

    spatial_data = get_data(input)

    if spatial_data is None:
        return
    
    plot_de = sc.pl.rank_genes_groups_dotplot(
        spatial_data['a'],
        n_genes = input.slider_n_genes(),
        key = 'rank_' + input.select_resolution(),
        min_logfoldchange=input.slider_lfc(),
        return_fig = True)

    main_ax_de = plot_de.get_axes()['mainplot_ax']
    x_axis_labels = [tick.get_text() for tick in main_ax_de.get_xticklabels()]
    short_labels_de = [label.split("-")[-1] for label in x_axis_labels]
    main_ax_de.set_xticklabels(short_labels_de)
    
    plt.subplots_adjust(top=1, bottom=0.2, left=0.05, right=0.95)
    plt.show()