#!/bin/bash
#SBATCH --job-name=callvars3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=150:00:00
#SBATCH --output=callvars3.out

module load bwa
module load samtools
module load bcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/subset1/

bcftools mpileup --threads 12 -f /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/asierrensis.scaffolded.fasta -q 20 -Q 20 *.bam \
| bcftools call --threads 12 -mv -Oz -o VCFFILE
