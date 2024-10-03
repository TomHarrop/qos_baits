#!/bin/bash

#SBATCH --job-name=qosbaits
#SBATCH --time=7-00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --output=sm.slurm.out
#SBATCH --error=sm.slurm.err
#SBATCH --partition=io

# Dependencies
module load apptainer/1.1.5-suid

# clobber broken apptainer config from module
export APPTAINER_BINDPATH=""

# Application specific commands:
export TMPDIR="${JOBDIR}"

printf "JOBDIR: %s\n" "${JOBDIR}"
printf "LOCALDIR: %s\n" "${LOCALDIR}"
printf "MEMDIR: %s\n" "${MEMDIR}"
printf "TMPDIR: %s\n" "${TMPDIR}"

snakemake \
	--profile petrichor_tmp \
	--keep-going \
	--retries 0 \
	--cores 64 \
	--local-cores 2 \
	output/020_overlaps/pzijinensis_min1000000.mega353.peakall/overlapping_loci.csv \
	output/020_overlaps/pzijinensis_min1000000.mega353.mega353/overlapping_loci.csv \
	output/020_overlaps/qos_min1000000.mega353.peakall/overlapping_loci.csv \
	output/020_overlaps/qos_min1000000.mega353.mega353/overlapping_loci.csv
