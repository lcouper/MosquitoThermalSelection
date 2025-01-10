# Sequence Analysis Workflow Tutorial 

The steps below outline how to repeat the analysis performed on the sequence data including: 
- processing raw reads through to variant calling
- calculating sample and population diversity metrics
- identifying and masking repeats in the reference genome
- estimating variance explained by SNPs


## Processing raw reads through to variant calling

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


## Calculating size of each chromosome

To obtain length in bp of each scaffold:
```
cat asierrensis.scaffolded.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```
<img width="190" alt="image" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/46f27809-91e6-46f9-a0ea-1cb8221cc616">


## Calculating sample and population diversity metrics

#### 1. Estimating heterozygosity 
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --keep controls.txt --het --out output_het
```
Note: in R, calculated observed and expected heterozygosity as:
```
hetero$O.HET <- (hetero$N_SITES - hetero$O.HOM)/hetero$N_SITES
hetero$E.HET <- (hetero$N_SITES - hetero$E.HOM)/hetero$N_SITES
```

#### 2. Estimating nucleotide diversity 
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --keep controls.txt --window-pi  10000 --out all_samples_pi
```

## Identifying and masking repeats in the reference genome

#### 1. Identify repeats 
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

#### 2. Mask repeats 
*Software used* RepeatMasker (Smit et al. 2021)    
*Script:* repeatmask.sbatch  
Note: Must use the older version of RepeatMasker for this to run properly (issues accessing databases in newer versions)**
```
module load repeatmasker/4.1.0
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
RepeatMasker -pa 16 -gff -lib sierrensis-families.fa  sierrensis_norepeats.fasta
```

## Estimating variance explained by SNPs 

#### 1. Calculate the genetic relationship matrix (GRM) from all the autosomal SNPs
*Software used:* GCTA (Yang et al. 2011)
Note: all SNPs here are likely autosomal, so use all identified SNPs.  
Note: Given errors with reading "1_RagTag" as chromosome names, prior to running the commande below, I had to alter pruneddata.bim file to replace "1_RagTag" to "1" (same for 2_RagTag and 3_RagTag).  
Help resources accessed: https://yanglab.westlake.edu.cn/software/gcta/#Tutorial

```
gcta64 --bfile pruneddata --autosome --make-grm --out pruneddata --autosome-num 3 --thread-num 12
gcta64 --bfile myplink --autosome --make-grm --out unpruneddata --autosome-num 3 --thread-num 12
```

#### 2. Estimate variance explained by the SNPs
*Software used:* GCTA (Yang et al. 2011)
Note: all SNPs are used, to avoid winners curse issue in GWA approaches.  
Note: steps below were conducted on both the LD pruned and unpruned dataset.   
Note: since the phenotype here is treatment (i.e., being in the control vs heat-selected group), GCTA considers this a case-control analysis 

```
gcta64 --grm pruneddata --pheno treat.txt --reml --out pruneddata --thread-num 12 --prevalence 0.46 --covar covars_gcta.txt
gcta64 --grm unpruneddata --pheno treat.txt --reml --out unpruneddata --thread-num 12 --prevalence 0.46 --covar covars_gcta.txt
# here prevalence of 0.46 refers to proportion of "cases" (i.e., heat-selected individuals) in sample, and covariates = sex
```

Results for LD-pruned dataset: 
<img width="240" alt="image" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/fdf6ac92-35d5-4ac2-8b0d-9c80d5aca4d9">

Results for unpruned dataset:
<img width="243" alt="image" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/0b72fb93-8404-480f-95e7-6c0692fce5d6">


#### 3. Repeat for adult thermal knockdown resistance, using only the significant SNPs

```
# keep only significant SNPs
gcta64 --bfile pruneddata --autosome --make-grm --extract KD_snplist.txt --out pruneddata_subset2 --autosome-num 3 --thread-num 12

gcta64 --grm pruneddata_subset2 --pheno KD.phe.txt --reml --out pruneddata_subset2 --thread-num 12 --covar covars.txt
# covariates here are sex and treatment
```
