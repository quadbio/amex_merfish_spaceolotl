[![pytest](https://github.com/quadbio/amex_merfish_spaceolotl/actions/workflows/install.yml/badge.svg)](https://github.com/quadbio/amex_merfish_spaceolotl/actions/workflows/install.yml)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue)

## Requirements
- A python package manager such as conda or mamba. Installing mamba is recommended as it is lightweight and fast. More information can be found [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) and [here](https://github.com/conda-forge/miniforge). Verify by typing `conda` or `mamba` in the terminal. Help messages should appear.

## Install instructions
**a. On a remote server:** Connect to the server (e.g. through SSH)\
**b. On a local machine (MacOS):** Open a terminal (e.g by typing "terminal" in SpotlightSearch)\
**c. On a local machine (Windows):** Press `Windows key` + `X`, select Windows Terminal

Navigate to a suitable location (a directory called amex_merfish_spaceolotl will be placed in there) and run the following commands:
```
# Create a new conda environment:
mamba create -n spaceolotl python=3.12

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
Note: the path to the data directory on the server should be provided as an absolute link (starting with a /). If you don't know the path, enter the data directory you created during the installation process and and type in `echo $PWD`.\

Upon login to the server, you should find the tarball in amex_merfish_spaceolotl/data. We can now remove unpack the data into the data directory and remove the tarball as it is no longer needed
```
# Unpack the data and clean up
tar -xzvf data.tar.gz && rm data.tar.gz
```
**b. On a local machine (MacOS):**
Drag and drop the tarball you downloaded from Polybox to the amex_merfish_spaceolotl directory. Double click on the tarball to decompress the object, this will create an ordinary folder called "data". After this process has finished, you can delete the `data.tar.gz`.

**c. On a local machine (Windows):** Will be added soon
 
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

**b. On a local machine (MacOS):**
Open a terminal like before and execute the following commands
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

**c. On a local machine (Windows):**
Will be added soon
