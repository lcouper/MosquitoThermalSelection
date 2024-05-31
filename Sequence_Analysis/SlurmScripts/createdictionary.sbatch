#!/bin/bash
#SBATCH --job-name=createdict
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=6:00:00
#SBATCH --output=createdict.out

module load bwa
module load java

cd /labs/emordeca/ThermalSelectionExpSeqFiles/

java -jar picard.jar CreateSequenceDictionary \
-R ref_genome/asierrensis.scaffolded.fasta \
-O ref_genome/asierrensis.dict
