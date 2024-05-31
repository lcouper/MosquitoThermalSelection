#!/bin/bash
#SBATCH --job-name=genotype_matrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=14:00:00
#SBATCH --output=genotype_matrix.out

module load bwa
module load bcftools
module load vcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF

vcftools --012 --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --out output_geno_maf01.vcf 

