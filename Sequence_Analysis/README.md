# Sequence Analysis Workflow 

#### 1. Obtained raw reads from Stanford Genomic Sequencing Center   
Included 470 .fastq.gz files (1 forward, 1 reverse for each of 235 samples)

#### 2. Filter adapter and trim low quality reads using trimmomatic
Trimmomatic V 0.39 (Bolger et al. 2014)  
Script: trim.sbatch

#### 3. Perform quality check on samples using fastqc
Script: [ NA ]

#### 4. Index reference genome with bwa
Note: only need to do this once   
Uses scaffolded *Ae. sierrensis* reference genome generated using Aedes aegypti AaegL5 genome and program RagTag (Alonge et al. 2022)   
Code below run on Stanford SCG Genomics Cluster  
```
gunzip ref_genome/asierrensis.V1.asm.bp.p_ctg.fasta_ragtag.scaffold.fasta.gz
module load bwa
bwa index ref_genome/asierrensis.V1.asm.bp.p_ctg.fasta_ragtag.scaffold.fasta
```
#### 5. Align sample reads to reference genome using bwa
Note: With computational resources noted in script, 1 samples took ~5 hours to run   
This produes a .sam file -- a tab-delimited text file containing information for each individual read and its alignment to the genome
Script: alignscaffast.sbatch

#### 6. Compress .sam to .bam using samtools
.bam is the compressed binary version of .sam, and enables for more efficient processing
Script: samtobam.sbatch

#### 7. Sort bam file by coordinates using samtools
Script: sortbam.sbatch

#### 8. Obtain summary stats about bam file 
Includes % mapped
```
samtools flagstat results/bam/P-050_scaffast.aligned.sorted.bam > results/bam/P-050_scaffast.aligned.sorted.bam.stats.txt
```

#### 9. Mark and remove duplicates using picard
Note picard.jar was downloaded from [the Broad Institute](https://broadinstitute.github.io/picard/)
Script: markdups.sbatch

#### 10. Index de-duplicated bam files using samtools
Script: indexbam.sbatch

#### 11. Compute depth at each position of sample
Script: computedepth.sbatch
This creates a txt file where the second and third columns are the position and coverage, respectively.  
To calculate the mean depth from this file 
```
awk 'BEGIN { total = 0; count = 0 } { total += $3; count += 1; } END { avg = total / count; print avg} ' results/bam/P-050_aligned_sorted_depth.txt
```

OR USE coverage_bcf.sbatch (determine which is better once this script finishes running)

#### 12. Compute alignment statistics using bamtools
Script: alignstats.sbatch

#### 13. Detect single nucleotide variants (SNVs) using bcftools
