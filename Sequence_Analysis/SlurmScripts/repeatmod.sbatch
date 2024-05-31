#!/bin/bash
#SBATCH --job-name=repeatmod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=96:00:00
#SBATCH --output=repeatmod.out

module load repeatmodeler/2.0.1

cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome

RepeatModeler -database sierrensis -pa 16 &>run2.out
