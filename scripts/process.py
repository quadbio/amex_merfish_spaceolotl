import tqdm
import pandas as pd

import sys
import os
path = os.path.abspath("../modules")
sys.path.insert(0, path)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import argparse
from ShinyMerfish import ShinyMerfish

def main():

    parser = argparse.ArgumentParser(
        description = "Processing quantified MERFISH data for display in PyShiny"
    )

    parser.add_argument(
        '-i',
        type=str,
        required=True,
        help='Path to the input directory',
        metavar='INPUT_DIR'
    )

    parser.add_argument(
       '-o',
       type=str,
       required=True,
       help='Path to the output directory',
       metavar='OUTPUT_DIR'
    )

    parser.add_argument(
       '-m',
       type=str,
       required=True,
       help='Path to the metadata table',
       metavar='METADATA_TABLE'
    )

    args = parser.parse_args()

    # Get the names of the samples
    samples = sorted([sample for sample in os.listdir(args.i) if 'region' in sample])

    with open(os.path.join(args.o, "data.txt"), "w") as f:
        f.write("\n".join(samples))

    # Process the data
    for sample in tqdm.tqdm(samples):
        ShinyMerfish(
            path = args.i,
            metadata_path = args.m,
            output_path = args.o,
            name = sample,
            normalize_to = 300,
            regress_out = ['nCount_RNA', 'x', 'y'],
            resolutions = [0.1, 0.2, 0.5, 1.0, 1.5],
            de = True,
            n_umap = 3,
            basis = 'X_pca',
            obs_keys = ['nCount_RNA', 'x', 'y']
        )

if __name__ == "__main__":
    main()




    