#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=72:00:00
#SBATCH --output=trim.out

cd /labs/emordeca/ThermalSelectionExpSeqFiles/

java -jar trimmomatic.jar PE E-09_S5_L001_R2_001.fastq.gz E-09_S5_L001_R2_001.fastq.gz \
E-09_S5_L001_R1_001.trim.fastq.gz E-09_S5_L001_R1_001un.trim.fastq.gz \
E-09_S5_L001_R2_001.trim.fastq.gz E-09_S5_L001_R2_001un.trim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15

