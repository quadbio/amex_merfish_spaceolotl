# General imports
import os
import json
import numpy as np
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

        # Create output
        self._export(output_path)

    def _set_resolution_keys(self, kwargs: dict) -> None:
        """Set resolution keys for clustering based on provided resolutions or defaults."""
    
        resolutions = kwargs.get('resolutions') or inspect.signature(scwf).parameters['resolutions'].default
        self.resolution_keys = [f'leiden_{res}' for res in resolutions]

    def _process_gdf(self) -> None:
        """Process and rotate the gdf, adding clustering information from adata."""

        self.gdf = sd.match_element_to_table(
            sdata = self.sdata,
            element_name = "cell_boundaries",
            table_name = "table"
        )[0]["cell_boundaries"]

        self.gdf['x_center'] = self.gdf['geometry'].apply(lambda row : row.centroid.x)
        self.gdf['y_center'] = self.gdf['geometry'].apply(lambda row : row.centroid.y)

        shift_x = -round(np.mean(self.gdf['x_center']))
        shift_y = -round(np.mean(self.gdf['y_center']))

        self.gdf['geometry'] = self.gdf['geometry'].translate(xoff = shift_x, yoff = shift_y)
        self.gdf['geometry'] = self.gdf['geometry'].rotate(self.angle, origin=(0, 0))

        self.sdata['table'].obs['x'] = self.gdf['x_center']
        self.sdata['table'].obs['y'] = self.gdf['y_center']
        
        clustering_columns = self.sdata['table'].obs[self.resolution_keys]
        self.gdf = self.gdf.join(clustering_columns)

    def _generate_traces(self) -> None:
        """Generate trace data for each resolution key and store cell coordinates."""

        traces = {}
        for resolution in self.resolution_keys:
            cluster_dict = {}
            for cluster_id in np.unique(self.gdf[resolution]):
                cluster = self.gdf.loc[self.gdf[resolution] == cluster_id, 'geometry']
                cluster = cluster.geometry.buffer(3).union_all().buffer(-3)

                x_coords, y_coords = extract_polygon_coords(cluster)
                cluster_dict[cluster_id] = {'x': x_coords, 'y': y_coords}
                
            traces[resolution] = cluster_dict

        cells_x_coords, cells_y_coords = extract_polygon_coords(self.gdf['geometry'], simplify = 1)
        traces['cells'] = {'x': cells_x_coords, 'y': cells_y_coords}

        self.traces = traces

    def _export(self, output_path) -> None:

        with open(os.path.join(output_path, self.sample + '_traces.json'), 'w') as f:
            json.dump(self.traces, f)
            
        self.sdata['table'].write_h5ad(os.path.join(output_path, self.sample + '.h5ad'))