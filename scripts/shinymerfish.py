# General imports
import os
import json
import numpy as np
import geopandas as gpd
import scanpy as sc
import inspect
import spatialdata as sd

# Project specific imports
from utility_functions.dissociated.pp.preprocessing import scwf
from utility_functions.spatial.ut.geopandas_utils import extract_polygon_coords
from amex_merfish_development._constants import *

# This is to catch sd table overwrite in memory warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# Project specific classes
class ShinyMerfish:
    """
    A class to prepare data for import into shiny app.

    Attributes:
        sample (str): The sample identifier.
        sdata (SpatialData): The spatial data read from the input path.
        angle (float): The rotation angle for the sample.
        label (str): The label for the sample.
        resolution_keys (list): The keys for different clustering resolutions.
        gdf (GeoDataFrame): The processed GeoDataFrame with clustering information.

    Methods:
        _set_resolution_keys(kwargs): Sets the resolution keys for clustering.
        _process_gdf(): Processes and rotates the GeoDataFrame.
        _generate_traces(): Generates traces for visualization.
    """
    
    def __init__(self, sample, input_path, output_path, metadata, **kwargs):

        self.sample = sample
        self.sdata = sd.read_zarr(input_path / self.sample / 'sdata.zarr')
        self.angle = metadata.loc[self.sample, 'angle'] + 180
        self.label = metadata.loc[self.sample, 'label']

        # Standard sc processing:
        self.sdata['table'] = scwf(self.sdata['table'], **kwargs).copy()

        # Define resolution keys and process data
        self._set_resolution_keys(kwargs)
        self._process_gdf()
        self._generate_traces()

        # Differential expression analysis
        self._differential_expression(self.resolution_keys)

        # Create output
        self._export(output_path)

    def _set_resolution_keys(self, kwargs: dict) -> None:
        """Set resolution keys for clustering based on provided resolutions or defaults."""
    
        resolutions = kwargs.get('resolutions') or inspect.signature(scwf).parameters['resolutions'].default
        self.resolution_keys = [f'leiden_{res}' for res in resolutions]

    def _get_centroid(self, gdf: gpd.GeoDataFrame, column: str, axis: str) -> pd.Series:
        """
        Get the centroid of a GeoDataFrame column.

        Parameters:
            gdf (GeoDataFrame): The GeoDataFrame containing the geometries.
            column (str): The column name containing the geometries.
            axis (str): The axis for which to get the centroid ('x' or 'y').

        Returns:
            Series: The centroid coordinates for the specified axis.
        """
        if axis == 'x':
            return gdf[column].apply(lambda row : row.centroid.x)
        elif axis == 'y':
            return gdf[column].apply(lambda row : row.centroid.y)

    def _process_gdf(self) -> None:
        """Process and rotate the gdf, adding clustering information from adata."""

        self.gdf = sd.match_element_to_table(
            sdata = self.sdata,
            element_name = "cell_boundaries",
            table_name = "table"
        )[0]["cell_boundaries"]

        # Rotate and translate the gdf
        self.gdf['geometry'] = self.gdf['geometry'].translate(
            xoff = -round(self._get_centroid(self.gdf, 'geometry', 'x').mean()),
            yoff = -round(self._get_centroid(self.gdf, 'geometry', 'y').mean())
            )
        self.gdf['geometry'] = self.gdf['geometry'].rotate(self.angle, origin=(0, 0))

        # Match the centroid of the gdf to the table
        self.sdata['table'].obs['x'] = self._get_centroid(self.gdf, 'geometry', 'x')
        self.sdata['table'].obs['y'] = self._get_centroid(self.gdf, 'geometry', 'y')
        
        # Add clustering columns to the gdf
        clustering_columns = self.sdata['table'].obs[self.resolution_keys]
        self.gdf = self.gdf.join(clustering_columns)

    def _generate_traces(self) -> None:
        """Generate trace data for each resolution key and store cell coordinates."""

        traces = {}
        for resolution in self.resolution_keys:
            cluster_dict = {}
            for cluster_id in np.unique(self.gdf[resolution]):
                cluster = self.gdf.loc[self.gdf[resolution] == cluster_id, 'geometry'].copy()
                cluster = cluster.geometry.buffer(3).union_all().buffer(-3)

                x_coords, y_coords = extract_polygon_coords(cluster)
                cluster_dict[cluster_id] = {'x': x_coords, 'y': y_coords}
                
            traces[resolution] = cluster_dict

        cells_x_coords, cells_y_coords = extract_polygon_coords(self.gdf['geometry'], simplify = 1)
        traces['cells'] = {'x': cells_x_coords, 'y': cells_y_coords}

        self.traces = traces
        self.sdata['table'].uns['traces'] = self.traces

    def _differential_expression(self, resolutions: list) -> None:
        """Perform differential expression analysis for each resolution key."""

        for resolution in resolutions:
            sc.tl.rank_genes_groups(self.sdata['table'], groupby=resolution, key_added= 'rank_' + resolution)

    def _export(self, output_path: str) -> None:
        """Export the processed data to the output path."""
            
        self.sdata['table'].write_h5ad(os.path.join(output_path, self.sample + '.h5ad'))