#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=4:00:00
#SBATCH --output=bgzip.out

module load bwa
module load samtools
module load bcftools
module load tabix

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/

bgzip -d Unfiltered_VCF_All.vcf.gz

