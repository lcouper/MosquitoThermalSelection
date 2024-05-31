# Sequence Analysis Workflow 

The steps below outline the analysis performed on the sequence data, from raw reads to variant calling. It also includes the sample and population diversity metric calculation steps, as well as steps to identify and mask repeats in the reference genome. 


#### 1. Obtained raw reads from Stanford Genomic Sequencing Center   
Included 470 .fastq.gz files (1 forward, 1 reverse for each of 235 samples)

#### 2. Filter adapter and trim low quality reads using trimmomatic
*Software used:* Trimmomatic V 0.39 (Bolger et al. 2014)  
*Script:* trim.sbatch  
Code snippet for single sample: 
```
java -jar trimmomatic.jar PE E-014_S10_L001_R1_001.fastq.gz E-014_S10_L001_R2_001.fastq.gz \
E-014_S10_L001_R1_001.trim.fastq.gz E-014_S10_L001_R1_001un.trim.fastq.gz \
E-014_S10_L001_R2_001.trim.fastq.gz E-014_S10_L001_R2_001un.trim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```

#### 3. Perform quality check on samples using fastqc
*Script:* fastqc.sbatch   
Code snippet for single sample:
```
fastqc trimmed_fastqc/*.fastq*
```

#### 4. Index reference genome with bwa
*Software used*: RagTag (Alonge et al. 2022)   
Note: only need to do this step once. Uses scaffolded *Ae. sierrensis* reference genome generated using *Aedes aegypti* AaegL5 genome
Code snippet:
```
gunzip ref_genome/asierrensis.scaffolded.fasta.gz
module load bwa
bwa index ref_genome/asierrensis.scaffolded.fasta
```

#### 5. Align sample reads to reference genome using bwa
Note: With computational resources noted in script, 1 sample took ~5 hours to run   
*Script:* alignreads.sbatch
Code snippet for single sample:
```
bwa mem -t 12 ref_genome/asierrensis.scaffolded.fasta \
trimmed_fastq/E-014_S10_L001_R1_001.trim.fastq trimmed_fastq/E-014_S10_L001_R2_001.trim.fastq > results/sam/E-014.aligned.sam
```

#### 6. Compress .sam to .bam using samtools
*Script:* samtobam2.sbatch 
Code snippet for single sample:
```
samtools view -S -b results/sam/E-014.aligned.sam -o results/bam/E-014.aligned.bam
```

#### 7. Sort bam file by coordinates using samtools
*Script:* sortbam.sbatch   
Code snippet for single sample:
```
samtools sort -o results/bam/E-014.aligned.sorted.bam results/bam/E-014.aligned.bam
```

#### 8. Obtain summary stats about bam file 
Note: Includes % mapped    
Code snippet for single sample:
```
samtools flagstat results/bam/E-014.aligned.sorted.bam > results/bam/E-014.aligned.sorted.bam.stats.txt
```

#### 9. Mark and remove duplicates 
*Software used:* picard (Broad Institute)
Note: picard.jar was downloaded from [the Broad Institute](https://broadinstitute.github.io/picard/)
*Script:* markdups.sbatch
Code snippet for single sample:
```
java -jar picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I=results/bam/E-014.aligned.sorted.bam \
O=results/bam/E-014.aligned.sorted.deduped.bam \
M=results/bam/E-014_marked_dup_metrics.txt
```

#### 10. Index de-duplicated bam files 
*Script:* indexbam.sbatch  
Code snippet for single sample:
```
samtools index -b results/bam/E-014.aligned.sorted.deduped.bam
```

#### 11. Compute depth at each position of sample
*Script:* computedepth.sbatch  
Note: This creates a txt file where the second and third columns are the position and coverage, respectively.   
To calculate the mean depth from this file:
```
awk 'BEGIN { total = 0; count = 0 } { total += $3; count += 1; } END { avg = total / count; print avg} ' results/bam/P-050_aligned_sorted_depth.txt
```

#### 12. Compute alignment statistics
Note: This calculates statistics including total reads, mapped reads, % failed QC, % duplicates, % paired-end reads, % singletons 
*Script:* alignstats.sbatch
Code snippet for single sample:
```
bamtools stats -in results/bam/E-014.aligned.sorted.deduped.bam > results/bam/E-014_aligned_sorted_AlignStats.txt
```

#### 13. Detect single nucleotide variants (SNVs) 
Note: This process takes too long to perform on all .bam files at once. Therefore, this script was run on subsets of ~20 samples at a time. Then the resulting VCFFILEs were concatenated (see step 17)
*Script:* callvariantsSx.sbatch
```
bcftools mpileup --threads 12 -f /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/asierrensis.scaffolded.fasta -q 20 -Q 20 *.bam \
| bcftools call --threads 12 -mv -Oz -o VCFFILE
```

#### 14. First pass filter SNVs
Note: For each subset, discard all SNVs with QUAL < 30
*Script:* filterSNPs.sbatch 
```
vcftools --gzvcf VCFFILE --minQ 30 --recode --recode-INFO-all --out subset1.vcf
```

#### 15. bgzip files
*Script:* bgzip.sbatch
```
module load tabix
cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/subset1
bgzip *.vcf
```

#### 16. Generate index for all vcf files
*Script:* indexvcf.sbatch
```
module load tabix
cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/subset1
tabix *.vcf
```

#### 17. Merge vcf files generated from sample subsets
Note: merging, rather than concatenating is appropriate here since the vcf subsets were from different samples, not different portions of the genome.
*Script:* bcfmergeAll.sbatch
```
cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/
bcftools merge *.vcf.gz > Unfiltered_VCF_All.vcf
```

#### 18. Create a dictionary for reference genome
*Software used:* picard
*Script:* createdictionary.sbatch
```
module load java
cd /labs/emordeca/ThermalSelectionExpSeqFiles/
java -jar picard.jar CreateSequenceDictionary \
-R ref_genome/asierrensis.scaffolded.fasta \
-O ref_genome/asierrensis.dict
```

#### 19. Sort VCF according to reference genome dictionary
*Script:* sortvcf.sbatch
```
module load bwa
module load java
cd /labs/emordeca/ThermalSelectionExpSeqFiles/
java -jar picard.jar SortVcf \
-I /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/Unfiltered_VCF_All.vcf \
-O /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/initialfilter_vcffiles/Unfiltered_VCF_All_sorted.vcf \
-SD /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/asierrensis.dict
```

#### 20. Filter SNVs 
Note: Discard all SNVs with QUAL < 30, Minor Allele Frequency of 0.05, Minimum Depth of 10x, and a Maximum Variant Missing of 0.995. Here, also remove 4 individuals with sequencing errors (E-05, E-06, E-07, E-09)
*Script:* filterSNPs_round2_0.995.sbatch 
```
vcftools --vcf Unfiltered_VCF_All_sorted.vcf --maf 0.05 --minQ 40 --max-missing 0.995 --minDP 10 --remove Removes.txt --recode --recode-INFO-all --out Filtered_VCF_All_sorted.vcf
```
*After filtering, kept 4,290,185 out of a possible 140426729 Sites*

#### 21. Remove mutli-allelic sites 
Note: Keep only bi-allelic sites for downstream analysis
*Script:* biallelic_only_0.995.sbatch
```
bcftools view -m2 -M2 -v snps Filtered_VCF_All_sorted_indexed.vcf > Filtered_VCF_All_sorted_0.995_bialleliconly.vcf
```
*Number of SNPs retained: 3,691,363*

#### 22. Generate genotype matrix 
Note: this outputs 3 files: ‘.012’ contains the genotypes of each individual on a separate line (with 0, 1, 2 denoting the number of non-reference alleles at the site), ‘.ind’ lists the individuals included in the main file, ‘.pos’ details the site location included in the main file. 
*Script:* genotype_matrix.sbatch
```
vcftools --012 --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --out output_geno.vcf
```

## Sample and population diversity metrics 

#### Estimating heterozygosity 
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --keep controls.txt --het --out output_het
```
Note: in R, calculated observed and expected heterozygosity as:
```
hetero$O.HET <- (hetero$N_SITES - hetero$O.HOM)/hetero$N_SITES
hetero$E.HET <- (hetero$N_SITES - hetero$E.HOM)/hetero$N_SITES
```

#### Estimating nucleotide diversity 
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --keep controls.txt --window-pi  10000 --out all_samples_pi
```

## Identify and mask repeats in reference genome 

#### **Step 1:** Identify repeats 
*Software used*: RepeatModeler v 2.0.1 (Flynn et al. 2020)
Note: Using NCBI and Dfam as database for repeats. Takes several days to run (includes 6 rounds of searching for repeats)  
*Script*: repeatmod.sbatch
```
module load repeatmodeler/2.0.1
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
BuildDatabase -name sierrensis -engine ncbi asierrensis.scaffolded.fasta
RepeatModeler -database sierrensis -pa 16 -LTRStruct &>run2.out  # -pa option similar to 'threads' in later versions
```
Note: the above creates 3 files:
sierrensis-families.fa (consensus sequences), sierensis-families.stk (seed alignments), and a log file
Sample stats from run shown here:

<img width="553" alt="image" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a8c266d8-79ca-4be1-8902-03bee38c7799">

#### **Step 2:** Mask repeats 
*Software used* RepeatMasker (Smit et al. 2021)
*Script:* repeatmask.sbatch
Note: Must use the older version of RepeatMasker for this to run properly (issues accessing databases in newer versions)**
```
module load repeatmasker/4.1.0
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
RepeatMasker -pa 16 -gff -lib sierrensis-families.fa  sierrensis_norepeats.fasta
```
