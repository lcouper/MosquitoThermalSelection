#!/bin/bash
#SBATCH --job-name=biallelic_filter_0.995
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=6:00:00
#SBATCH --output=biallelic_filter_0.995.out

module load bwa
module load bcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF

bcftools view -m2 -M2 -v snps Filtered_VCF_All_sorted_0.995.vcf > Filtered_VCF_All_sorted_0.995_bialleliconly.vcf

