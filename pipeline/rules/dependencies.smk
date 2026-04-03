########
# Import required packages
########

import pandas as pd
from pathlib import Path
import multiprocessing
import hashlib
import sys
import tempfile
import os
import subprocess
import re
import glob

######
# Check genodisc conda env is activated
######

conda_env_name = os.getenv('CONDA_DEFAULT_ENV')
if not conda_env_name == 'genodisc':
  print("Error: The genodisc conda environment must be active when running the pipeline.\nFor more information: https://opain.github.io/GenoDisc/pipeline_readme.html")
  sys.exit(1)

######
# Check config file
######

# Check for missing required configuration parameters
required_config_params = ['outdir', 'config_file', 'gwas_list']
missing_or_invalid_config_params = [param for param in required_config_params if param not in config or config[param] == 'NA']

if missing_or_invalid_config_params:
    # Print an informative message
    print(f"Missing or invalid (set to 'NA') required configuration parameters: {', '.join(missing_or_invalid_config_params)}. Please specify these in the configuration file.")

    # Exit Snakemake gracefully
    sys.exit(1)

# Set outdir parameter
outdir=config['outdir']

# Set resource directory
resdir = config.get('resdir', None)
if resdir is None or resdir == 'NA':
    resdir = 'resources'

# Set chromosomes to analyse
chromosomes = config.get("chromosomes", list(range(1, 23)))

########
# Create required functions
########

def get_current_version():
    cmd = "git describe --tags"
    tag = subprocess.check_output(cmd, shell=True).decode().strip()
    match = re.match(r"v?(\d+)\.(\d+)", tag)
    if match:
        return int(match.group(1)), int(match.group(2))  # Major, Minor
    else:
        raise ValueError("Git tag does not contain a valid version format.")

def read_last_version():
    if os.path.exists(last_version_file):
        with open(last_version_file, "r") as file:
            major, minor = file.read().strip().split('.')
            return int(major), int(minor)
    return 0, 0  # Default to 0.0 if file does not exist

def write_last_version(major, minor):
    with open(last_version_file, "w") as file:
        file.write(f"{major}.{minor}")

########
# Check for repo version updates
########

# If there has been a change to the major or minor version numbers, we will rerun the entire pipeline

# Define the path for storing the last known version
os.makedirs(resdir, exist_ok=True)
last_version_file = f"{resdir}/last_version.txt"

# Access overwrite flag from config
overwrite = config.get("overwrite", "false").lower() == "true"

# Main logic to check version and decide on execution
current_major, current_minor = get_current_version()
last_major, last_minor = read_last_version()

# Check if the last version is 0.0, proceed without requiring overwrite
if last_major == 0 and last_minor == 0:
    print(f"Initial version setup detected. Updating to v{current_major}.{current_minor}.")
    write_last_version(current_major, current_minor)
else:
    # Check for both major and minor version changes
    if current_major != last_major or current_minor != last_minor:
        if not overwrite:
            print(f"Change in version of GenoDisc detected from v{last_major}.{last_minor} to v{current_major}.{current_minor}. Use --config overwrite=true to proceed.")
            sys.exit(1)
        else:
            print("Proceeding with version update due to overwrite=true config.")
            write_last_version(current_major, current_minor)  # Update the stored version

####
# Download BioMart gene annotations
####

rule download_biomart:
  output:
    f"{resdir}/data/biomart/biomart_genes_grch37.tsv"
  params:
    resdir=resdir
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/download_biomart.R --resdir {params.resdir}"

