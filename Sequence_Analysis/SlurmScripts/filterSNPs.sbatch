#!/bin/bash
#SBATCH --job-name=filterSNPs2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=8:00:00
#SBATCH --output=filterSNPs2.out

module load bwa
module load samtools
module load bcftools
module load vcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_vcffiles/

vcftools --vcf AllSamples_Sorted_WithAF.vcf --maf 0.05 --minDP 10 --recode --recode-INFO-all --out AllSamples_FinalSorted.vcf
