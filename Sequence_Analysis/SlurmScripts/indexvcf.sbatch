#!/bin/bash
#SBATCH --job-name=indexvcf2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=8:00:00
#SBATCH --output=indexvcf2.out

module load bwa
module load samtools
module load bcftools
module load tabix

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/subset1

tabix subset1.vcf.recode.vcf.gz

