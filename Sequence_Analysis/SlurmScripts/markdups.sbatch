#!/bin/bash
#SBATCH --job-name=markdups2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=36:00:00
#SBATCH --output=markdups2.out

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/sortedbams/

module load java
module load bwa

for infile in *.aligned.aligned.sorted.bam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.aligned.sorted.bam)
java -jar picard.jar MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I ${base}.aligned.aligned.sorted.bam \
-O "${base}.deduped.bam" \
-M "${base}.dup_metrics.txt"
done
