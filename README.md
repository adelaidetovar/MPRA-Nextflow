# MPRA-Nextflow

## Installation

### Clone repository

    git clone https://github.com/adelaidetovar/MPRA-Nextflow.git

### Install dependencies
All required software is provided by a Singularity container within the pipeline. It is recommended to install Nextflow (version 23.1 or greater) and, if applicable, Singularity, within a new conda environment.

### Get a cloud.sylabs.io access token to access Singularity image
If you don't already have a working sylabs.io access token, create or login to your [https://cloud.sylabs.io/dashboard cloud.sylabs.io] account. In your dashboard under 'Access Tokens' create a new token. In your shell, type `singularity remote login` and paste the access token.

## Running the pipeline
This pipeline takes raw FASTQ files containing reads from sequencing runs that use custom or conventional Truseq sequencing primers. It can concatenate files from separate runs or across separate lanes.

At present, it has only been tested on a SLURM-configured server (specifically Greatlakes at the University of Michigan), which is reflected in the `nextflow.config` file. If running on any other server, make sure to update this file, all job submission scripts, and the`.nf` files to match your server. Specifically, you'll need to account names, paths to directories, and in the `.nf` files, adjust `queue` parameter to reflect the names of your server's partitions.

Make sure your sample names are in the format `SAMPLENAME_[RNA,gDNA,DNA,cDNA]`, e.g. `ATS0999_RNA` or `Sample8_gDNA`.

1. Make your project work directory.
2. Make symbolic links to your raw FASTQ files in a subdirectory. You can do this manually if you have just a single set of data; otherwise, use the example `make_symlink.sh` script to automatically make symlinks to your files if you have data from multiple sequencing runs.
3. Generate your input configuration file, by running `bin/make_config.py`.
4. Run the pipeline.

Example job submission scripts for steps 2-4 are available in `examples`.
