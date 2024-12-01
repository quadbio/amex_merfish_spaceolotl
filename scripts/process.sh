#!/bin/bash

METADATA="/nas/groups/treutlein/PROJECTS/axolotl/MERFISH/statistics/anno_extended.csv"
INPUT_PATH="/nas/groups/treutlein/PROJECTS/axolotl/MERFISH/quant/v1.0.2/"
OUTPUT_PATH="/home/aabouelela/projects/amex_spatial/spaceolotl/data/"

python process.py -i $INPUT_PATH -o $OUTPUT_PATH -m $METADATA