# MPRA-Nextflow

## Installation

### Clone repository

    git clone https://github.com/adelaidetovar/MPRA-Nextflow.git

### Install Nextflow and dependencies
All required software is provided by a Singularity container within the pipeline. It is recommended to install Nextflow within a new conda environment.

## Running the pipeline

This pipeline takes raw FASTQ files containing reads from sequencing runs that use custom or conventional sequencing primers and can concatenate files from separate runs or across separate lanes. At present, it has only been tested on a SLURM-configured server, which is reflected in the `nextflow.config` file.

To run the pipeline, submit as a job. The first step runs a script (`make-config.py`) to generate an `input-config.json` file containing sample-level metadata when given a list of sample names and directories where files are stored. The second step runs the NextFlow pipeline. Example formatted scripts and outputs are contained in the `examples` directory.
