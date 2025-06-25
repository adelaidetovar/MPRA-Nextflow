# MPRA-Nextflow

## Installation

### Clone repository

    git clone https://github.com/adelaidetovar/MPRA-Nextflow.git

### Install dependencies
All required software is provided by a Singularity container within the pipeline. It is recommended to install Nextflow (version 22.1 or greater) and, if applicable, Singularity, within a new conda environment.

## Running the pipeline
This pipeline takes raw FASTQ files containing reads from sequencing runs that use custom or conventional Truseq sequencing primers. It can concatenate files from separate runs or across separate lanes. At present, it has only been tested on a SLURM-configured server, which is reflected in the `nextflow.config` file. In the `.nf` files, adjust `queue` parameter to reflect the names of your server's partitions.

Make sure your sample names are in the format `SAMPLENAME_[RNA,gDNA,DNA,cDNA]`, e.g. `ATS0999_RNA` or `Sample8_gDNA`.

1. Make your project work directory.
2. Make symbolic links to your raw FASTQ files in a subdirectory, either manually or by using the example `make_symlink.sh` script to automatically make symlinks to your files (especially usefully when you have data across multiple sequencing runs). 
3. Generate your input configuration file, by running `bin/make_config.py`.
4. Run the pipeline.

Example job submission scripts for steps 2-3 are available in `examples`.
