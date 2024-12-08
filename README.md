[![pytest](https://github.com/quadbio/amex_merfish/actions/workflows/pytest.yml/badge.svg)](https://github.com/quadbio/amex_merfish/actions/workflows/pytest.yml)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue)

## Requirements
- A python package manager such as conda or mamba. Installing mamba is recommended as it is lightweight and fast. More information can be found [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) and [here](https://github.com/conda-forge/miniforge). Verify by typing `conda` or `mamba` in the terminal. Help messages should appear.
- For Windows users: A program to decompress `.tar.gz` object such as 7-zip or WinRAR

Note: The following commands are using `mamba`. If you are using `conda`, simply replace `mamba` with `conda`.

## Install instructions
**a. On a remote server:** Connect to the server (e.g. through SSH)\
**b. On a local machine (MacOS):** Open a terminal (e.g by typing "terminal" in SpotlightSearch)\
**c. On a local machine (Windows):** Search for "Miniforge Prompt" in the Windows search bar

Navigate to a suitable location (a directory called amex_merfish_spaceolotl will be placed in there) and run the following commands:
```
# Create a new conda environment (if git is already available on your system, you don't need to install it again):
mamba create -n spaceolotl python=3.12 git

# Activate the environment:
mamba activate spaceolotl

# Fetch the code from github. This will create a directory called 'amex_merfish_spaceolotl'
git clone --branch main https://github.com/quadbio/amex_merfish_spaceolotl.git

# Switch to that directory
cd amex_merfish_spaceolotl

# Install (in editable mode, package will incorporate changes if code is modified)
pip install -e .
```

The following step is only necessary if you are working on a remote server:
```
# Make a directory for the data
mkdir data
```

## Downloading the data
In order to fetch the data (provided as a tarball ending with tar.gz), download it to you local machine from Polybox via the link shared on Slack. Then proceed with either of the following options:

**a. On a remote server:** The data downloaded to your local machine needs to be transferred to the remote. Use your favourite File Transfer Protocol client program (e.g. FileZilla) or simply run the following command from your terminal on your local machine:
```
# This command transfers the tarball to the data directory. You might be asked to authenicate yourself
scp /path/to/tarball/on/local/machine username@domain:/path/to/amex_merfish_spaceolotl/data
```
Note: the path to the data directory on the server should be provided as an absolute link (starting with a /). If you don't know the path, enter the data directory you created during the installation process and and type in `echo $PWD`.

Upon login to the server, you should find the tarball in amex_merfish_spaceolotl/data. We can now remove unpack the data into the data directory and remove the tarball as it is no longer needed
```
# Unpack the data and clean up
tar -xzvf data.tar.gz && rm data.tar.gz
```
**b. On a local machine (MacOS):**
Drag and drop the tarball you downloaded from Polybox to the amex_merfish_spaceolotl directory. Double click on the tarball to decompress the object, this will create an ordinary folder called "data". After this process has finished, you can delete the `data.tar.gz`.

**c. On a local machine (Windows):** 
Move the tarball you downloaded from Polybox to the amex_merfish_spaceolotl directory. Right-click on the tarball and select "extract to data". This will create an ordinary folder called "data". After this process has finished, you can delete the `data.tar.gz` (might be displayed as `data.tar`).

## Running the app
**a. On a remote server:**\
Connect to the remote server while forwarding a port:
```
# Make sure the ports are not already in use
ssh -L 12345:localhost:7600 username@domain
```
Launch the app like so:
```
# Change to the amex_merfish_spaceolotl directory
cd path/to/amex_merfish_spaceolotl

# Activate the conda environment
conda activate spaceolotl

# Run the app
shiny run --port 7600 src/spaceolotl/app.py
```

You can now access the app in your local browser by typing in `localhost:12345` in the search bar.

**b. On a local machine (MacOS or Windows):**
Open a `mamba`-enabled terminal like before and execute the following commands
```
# Move to the amex_merfish_spaceolotl directory
cd path/to/amex_merfish_spaceolotl

# Activate the Python environment
mamba activate spaceolotl

# Run the app
shiny run src/spaceolotl/app.py
```
You can now access the app through the search bar in your local browser by typing `localhost:8000`.
Note: There might be cases where the standard port is blocked and Shiny selects a different one. In this case you won't be able to access the app as described above. Check the commandline upon startup and search for a message like
```
INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
```
Copy the part `http://127.0.0.1:8000` into the search bar or your browser. This should lead you to the app.

## Using the app
`The UMAP and Spatial widgets`: Rendered with plotly, different tools are available from the top right corner of the plots. For the spatial plot, clusters or traces can be isolated by double clicking on the corresponding entry in the legend or they can be hidden by clicking once. 

`Select dataset`: Select a dataset. The data is tagged to a specific version of the processing workflow (visible in the header), that can be used to trace back observations once new workflows are developed. However, it is not expected, that the data will look much different in coming interations of the analysis\
`Cluster resolutions`: Select the leiden resolution for clustering\
`Show cell outlines`: Plot the cell outlines in the spatial plot. Toggled off by default as it increases the rendering time\
`Show clusters`: Whether or not the clusters should be visible\
`Select gene` and `Plot gene expression`: Select a gene of interest to be plotted on top of the UMAP and Spatial plot\
`Slider UMAP and Slider Space`: Point size of the dots in the UMAP or spatial plot respectively\
`Select genes`: Gene expression to plotted as a dotplot in the "Gene expression per cluster" dropdown panel. Multiple can be selected\
`Slider nGenes`: Number of DE genes to display per cluster in the "Differential gene expression" dropdown panel\
`Slider minLFC`: Filter DE genes in the "Differential gene expression" dropdown panel by this log-foldchange.
``
