#!/bin/bash
#SBATCH --job-name=alignreads2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=72:00:00
#SBATCH --output=alignreads2.out

cd /labs/emordeca/ThermalSelectionExpSeqFiles/

module load bwa/0.7.10

#first unzip trimmed fastq files if not done already
gunzip trimmed_fastq/temphide2/*.gz

#Then align to scaffolded aesierrensis ref genome

for infile in trimmed_fastq/temphide2/*_R1_001.trim.fastq
do
base=$(basename ${infile} _R1_001.trim.fastq)
bwa mem -t 12 ref_genome/asierrensis.scaffolded.fasta \
trimmed_fastq/temphide2/${base}_R1_001.trim.fastq trimmed_fastq/temphide2/${base}_R2_001.trim.fastq > results/sam/"${base}.aligned.sam"
done
