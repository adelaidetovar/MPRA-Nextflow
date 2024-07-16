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

python make_config.py --libnames sample1_DNA,sample2_DNA,sample1_RNA,sample2_RNA --dirs /path/to/data-a/,/path/to/data-b/
