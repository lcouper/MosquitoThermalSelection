#!/bin/bash
#SBATCH --job-name=computedepth2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=36:00:00
#SBATCH --output=computedepth2.out

module load bwa
module load samtools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/

for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools depth -a ${base}.bam > "${base}.depth.txt"
done
