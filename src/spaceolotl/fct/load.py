import os
import json
import scanpy as sc

from spaceolotl._constants import *

def get_data(input):
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