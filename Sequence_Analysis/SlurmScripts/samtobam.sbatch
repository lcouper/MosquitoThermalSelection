#!/bin/bash
#SBATCH --job-name=samtobam2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=36:00:00
#SBATCH --output=samtobam2.out

module load bwa
module load samtools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/

for infile in sam/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b sam/${base}.aligned.sam > bam/"${base}.aligned.bam"
done
