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

FQ_DIR="/path/to/in_fq/"

LIBNAMES=$(ls "$FQ_DIR" | grep -E "\.r[12]\.(fq|fastq)\.gz$" | \
    sed -E 's/(_[A-Z])?\.r[12]\.(fq|fastq)\.gz$//' | sort -u | tr '\n' ',' | sed 's/,$//')

python /path/to/bin/make_config.py --libnames "$LIBNAMES" --dir "$FQ_DIR"
