# General imports
import numpy as np
import time

# Project specific imports
from shinymerfish import ShinyMerfish
from amex_merfish_development._constants import *
from concurrent.futures import ThreadPoolExecutor

from utility_functions._logging import _setup_logger
logger = _setup_logger()

# Procext specific variables
IN_DIR = GET_DATA_DIR('v1.0.2')
OUT_DIR = '../data/'
samples = ANNO_BASE.index

# Project specific functions
def process_sample(sample):

    ShinyMerfish(
        sample,
        IN_DIR,
        OUT_DIR,
        ANNO_BASE,
        normalize_to=800,
        regress_out=['nCount_RNA'],
        resolutions=[0.1, 0.2, 0.5, 1.0, 1.5],
        n_umap=3,
        basis='X_pca',
        obs_keys=['nCount_RNA'],
        verbose = False
    )

if __name__ == '__main__':

    np.save(OUT_DIR + 'genes.npy', GENE_PANEL)
    np.save(OUT_DIR + 'samples.npy', samples)

    logger.info("Started processing")
    start_time = time.time()
    
    # TODO: There are some warnings accessing the zarr store with this multithreading approach
    # but so far no issues have been observed. Check later
    num_threads = min(len(samples), 30)
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        executor.map(process_sample, samples)
    
    end_time = time.time()
    logger.info(f"Finished processing. Execution: {round((end_time - start_time) / 60,1)} mins.")
