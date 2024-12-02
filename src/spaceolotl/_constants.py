from pathlib import Path
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / 'data'
LEIDEN_RESOLUTIONS = {'leiden_0.1': 0.1, 'leiden_0.2': 0.2, 'leiden_0.5': 0.5, 'leiden_1.0': 1.0, 'leiden_1.5': 1.5}
GENES = np.load(DATA_DIR / 'genes.npy', allow_pickle=True).tolist()
GENES_LABEL = [x.split('-')[1] for x in GENES]
DATA = np.load(DATA_DIR / 'samples.npy', allow_pickle=True).tolist()