import os
import json
from typing import Union

import numpy as np
import pandas as pd
import geopandas as gpd

import shapely
from shapely.geometry import MultiPolygon

import spatialdata as sd

from _logging import _setup_logger
logger = _setup_logger()

import scanpy as sc
import scanpy_utils

class ShinyMerfish:
    """
    A class to process and visualize spatial data from MERFISH experiments with clustering and rotation transformations.
    
    Parameters:
        path (str): Path to the main directory containing spatial data.
        metadata_path (str): Path to the metadata CSV file.
        name (str): Identifier for the dataset.
    """

    def __init__(self,
                 path: str,
                 metadata_path: str,
                 output_path: str,
                 name: str,
                 **kwargs):

        logger.info(f"Processing: {os.path.join(path, name)}")
        
        self.name = name
        self.sdata = sd.read_zarr(os.path.join(path, self.name, 'sdata.zarr'))

        metadata = pd.read_csv(metadata_path, index_col=0)
        self.angle = metadata.loc[self.name, 'angle']
        self.label = metadata.loc[self.name, 'label']

        # Standard single-cell processing
        # The assignment of the spatial coordinates is necessary to regress out
        # but is redundant with assignment after rotation and shifting. Needs to
        # be fixed eventually
        self.sdata['table'].obs['x'] = self.sdata['table'].obsm['spatial'][:,0]
        self.sdata['table'].obs['y'] = self.sdata['table'].obsm['spatial'][:,1]
        self.sdata['table'] = scanpy_utils.std_sc_wf(self.sdata['table'], **kwargs)

        # Define resolution keys and process data
        self._set_resolution_keys(kwargs)
        self._process_gdf()
        self._generate_traces()

        # Create output
        self._export(output_path)

    def _set_resolution_keys(self, kwargs: dict) -> None:
        """Set resolution keys for clustering based on provided resolutions or defaults."""
        
        resolutions = kwargs.get('resolutions') or inspect.signature(std_sc_wf).parameters['resolutions'].default
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

        # Re-evaluate the centers, this is horrible
        self.gdf['x_center'] = self.gdf['geometry'].apply(lambda row : row.centroid.x)
        self.gdf['y_center'] = self.gdf['geometry'].apply(lambda row : row.centroid.y)

        self.sdata['table'].obs['x'] = self.gdf['x_center']
        self.sdata['table'].obs['y'] = self.gdf['y_center']
        
        clustering_columns = self.sdata['table'].obs[self.resolution_keys]
        self.gdf = self.gdf.join(clustering_columns)

    def _extract_polygon_coords(self,
                                polygonal_shapes: Union[gpd.GeoSeries, MultiPolygon],
                                simplify = False
                               ) -> tuple:
        """Extract x, y coordinates from polygon shapes, separating individual polygons with NaNs."""
        
        if isinstance(polygonal_shapes, MultiPolygon):
            # This is the case where more than one disconnected patches form the cluster
            polygonal_shapes = gpd.GeoSeries(polygonal_shapes.geoms)
        else:
            # This is the case where there is only one polygon forming the cluster
            polygonal_shapes = gpd.GeoSeries(polygonal_shapes)

        if simplify:
            polygonal_shapes = polygonal_shapes.apply(lambda poly: shapely.simplify(poly, simplify, preserve_topology = True))
        
        coords = polygonal_shapes.apply(lambda poly: poly.exterior.xy)
        x_coords = np.concatenate([np.append(xs, None) for xs, ys in coords])[:-1]
        y_coords = np.concatenate([np.append(ys, None) for xs, ys in coords])[:-1]

        return list(x_coords), list(y_coords)
               
    def _generate_traces(self) -> None:
        """Generate trace data for each resolution key and store cell coordinates."""

        traces = {}
        for resolution in self.resolution_keys:
            cluster_dict = {}
            for cluster_id in np.unique(self.gdf[resolution]):
                cluster = self.gdf.loc[self.gdf[resolution] == cluster_id, 'geometry']
                cluster = cluster.geometry.buffer(3).union_all().buffer(-3)

                x_coords, y_coords = self._extract_polygon_coords(cluster, simplify = 1)
                cluster_dict[cluster_id] = {'x': x_coords, 'y': y_coords}
                
            traces[resolution] = cluster_dict

        cells_x_coords, cells_y_coords = self._extract_polygon_coords(
            self.gdf['geometry'], simplify = 1)
        traces['cells'] = {'x': cells_x_coords, 'y': cells_y_coords}

        self.traces = traces

    def _export(self, output_path) -> None:

        with open(os.path.join(output_path, self.name + '_traces.json'), 'w') as f:
            json.dump(self.traces, f)
            
        self.sdata['table'].write_h5ad(os.path.join(output_path, self.name + '.h5ad'))
        logger.info("Processing finished\n")

