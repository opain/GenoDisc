#!/usr/bin/Rscript

# cmapR cannot be installed from Bioconductor directly
# There are version/dependency issues that require manual preparations

# Download proto RProtoBufLib
# We have to download source code as BioConductor is installing old version for some reason
dir.create('resources/software/RProtoBufLib')
download.file(url='https://www.bioconductor.org/packages/release/bioc/src/contrib/RProtoBufLib_2.8.0.tar.gz',
              destfile='resources/software/RProtoBufLib/RProtoBufLib_2.8.0.tar.gz')

# Install RProtoBufLib from source
install.packages('resources/software/RProtoBufLib/RProtoBufLib_2.8.0.tar.gz', repos = NULL, type="source")

# Install cytolib from github
devtools::install_github("RGLab/cytolib")

# Install flowCore from github
devtools::install_github("RGLab/flowCore")

# Install cmapR
BiocManager::install("cmapR")

file.create('resources/software/install_cmapR.done')
