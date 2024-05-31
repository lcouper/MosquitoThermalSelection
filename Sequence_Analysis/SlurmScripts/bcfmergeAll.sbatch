#!/bin/bash
#SBATCH --job-name=bcfmergeAll
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=14:00:00
#SBATCH --output=bcfmergeAll.out

module load bwa
module load bcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/

bcftools merge *.vcf.gz > Unfiltered_VCF_All.vcf
