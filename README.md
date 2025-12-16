

# Boeken Tau aggregate interactome

This repository contains the analysis code associated with the tau aggregate interactome proteomics and single-molecule characterisation project, led by Dorothea Böken. This manuscript has been submitted for publication under the title "Remodelling of the Tau Interactome Across Aggregation Stages Reveals Stage-Specific Vulnerabilities in Alzheimer’s Disease".

This manuscript has been submitted as a preprint via BioRxiv [here](biorxiv/link). A link to the final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python 3 (v 3.10.5). For specific package requirements, see the environment_proteomics.yml and the environment_smma.yml file, or  create a new conda environment containing all packages by running ```conda create -f >>env_name<<```. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v 10.6](https://www.graphpad.com/scientific-software/prism/).

## Raw data

For convenience, example image files are provided here under the ```data``` folder. These data may be used to explore the workflows presented here as described below.

All proteomics data are availavle via [PRIDE](Pride/link).


## Workflow

Example raw images are provided here within the ```data``` folder to test the included analysis scripts.

Individual analysis and plotting scripts for proteomics and single-molecule imaging (co-localisation and super-resolution microscopy) are presented within the ```src``` folder. Where processing order is important for individual analyses, scripts have been numbered and should be run in order before unnumbered counterparts.


