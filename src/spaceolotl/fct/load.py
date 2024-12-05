import os
import json
import scanpy as sc

from spaceolotl._constants import *

def get_data(input):
        
    if not (name := input.select_dataset()):
        return None
    
    try:
        adata = sc.read_h5ad(os.path.join(DATA_DIR, name + '.h5ad'))
        return adata
    
    except (FileNotFoundError, IOError):
        print("File not found")
        return None