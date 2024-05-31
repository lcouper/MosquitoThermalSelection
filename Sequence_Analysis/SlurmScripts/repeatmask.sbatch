#!/bin/bash
#SBATCH --job-name=repeatmask
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=102:00:00
#SBATCH --output=repeatmask.out

module load repeatmasker/4.1.0

cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome

RepeatMasker -pa 16 -gff -lib sierrensis-families.fa asierrensis.scaffolded.fasta
