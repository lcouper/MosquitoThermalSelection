#!/bin/bash
#SBATCH --job-name=sortvcf2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=16:00:00
#SBATCH --output=sortvcf2.out

module load bwa
module load java

cd /labs/emordeca/ThermalSelectionExpSeqFiles/

java -jar picard.jar SortVcf \
-I /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/Unfiltered_VCF_All.vcf \
-O /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/Unfiltered_VCF_All_sorted_16h.vcf \
-SD /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/asierrensis.dict
