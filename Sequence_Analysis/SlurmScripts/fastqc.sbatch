#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=36:00:00
#SBATCH --output=fastqc.out

cd /labs/emordeca/ThermalSelectionExpSeqFiles/

module load fastqc
fastqc trimmed_fastq/*.fastq.gz

