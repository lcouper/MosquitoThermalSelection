#!/bin/bash
#SBATCH --job-name=filterSNPs_round2_take7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=16:00:00
#SBATCH --output=filterSNPs_round2_0.995.out

module load bwa
module load samtools
module load bcftools
module load vcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/

vcftools --vcf Unfiltered_VCF_All_sorted.vcf --maf 0.1 --minQ 40 --max-missing 0.995 --minDP 10 --remove Removes.txt --recode --recode-INFO-all --out Filtered_VCF_All_sorted_0.995_maf0.1.vcf
