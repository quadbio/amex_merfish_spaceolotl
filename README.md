## Requirements
- A python installation and a package manager such as conda or mamba.

## Install instructions
**a. On a remote server:**\
After connecting to the server and navigating to a suitable location, run the following commands:
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

# Make a directory for the data
mkdir data
```
The package is now installed and the app can be launched, but the data directory is currently empty.
In order to fetch the data (provided as a tarball ending with tar.gz), download to you local machine it from Polybox via the link shared on Slack.
Transfer the data using your favourite File Transfer Protocol client program (e.g. FileZilla) or simply through the command line:

```
# This command transfers the tarball to the data directory. You might be asked to authenicate yourself
scp /path/to/tarball/on/local/machine username@domain:/path/to/amex_merfish_spaceolotl/data
```

Upon login to the server, you should find the tarball in amex_merfish_spaceolotl/data. We can now remove unpack the data into the data directory and remove the tarball as it is no longer needed
```
# Unpack the data and clean up
tar -xzvf archive.tar.gz && rm archive.tar.gz
```

Congratulations! Everything is set up to run the app.

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
