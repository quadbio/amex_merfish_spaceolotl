# General imports
import numpy as np

# Project specific imports
from shinymerfish import ShinyMerfish
from amex_merfish_development._constants import *
from utility_functions._logging import _setup_logger
logger = _setup_logger()

# Procext specific variables
IN_DIR = GET_DATA_DIR('v1.0.2')
OUT_DIR = '../data/'
samples = ANNO_BASE.index

# Processing
np.save(OUT_DIR + 'genes.npy', GENE_PANEL)
np.save(OUT_DIR + 'samples.npy', samples)

for sample in samples:
    logger.info(f"Processing {sample}")
    ShinyMerfish(sample,
                 IN_DIR,
                 OUT_DIR,
                 ANNO_BASE,
                 normalize_to=800,
                 regress_out=['nCount_RNA'],
                 resolutions=[0.1, 0.2, 0.5, 1.0, 1.5],
                 n_umap=3,
                 basis='X_pca',
                 obs_keys=['nCount_RNA'])