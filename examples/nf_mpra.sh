#!/bin/bash
#SBATCH --time="12:00:00"
#SBATCH --mem=4000M
#SBATCH --output=/path/to/proj/logs/slurm-%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type END,FAIL
#SBATCH --signal=B:TERM@60
#SBATCH --job-name=nf_mpra
#SBATCH --mail-user user@email.com

module load openjdk
module load singularity/3.11.1
export WORKDIR=/path/to/proj/
exec nextflow run -resume -params-file $WORKDIR/input-config.json --outdir $WORKDIR/bc_count --check_samp_by_samp /home/user/github/MPRA-Nextflow/bc_count.nf
