# Steps to identify genes from reference genome sequence and assign SNPs to genes

## Identifying genes from reference genome sequence 

### Step 1. Identify repeats
- *Software used*: RepeatModeler2 2.0.1
- Note: Use RepeatModeler2 to construct a species-specific repeat library 
- RepeatModeler2 Can identify repeats de novo, But also requires an external database. Dfam is one that is free
- Help pages accessed: https://darencard.net/blog/2022-10-13-install-repeat-modeler-masker/
- https://github.com/Dfam-consortium/RepeatModeler
- https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Basic
- Interpret results from repeat modeler: https://github.com/Dfam-consortium/RepeatModeler
- *Script*: repeatmod.sbatch

```
module load repeatmodeler/2.0.1
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
BuildDatabase -name sierrensis asierrensis.scaffolded.fasta
RepeatModeler -database sierrensis -pa 16 &>run2.out
```


### Step 2. Mask repeats
- *Software used: RepeatMasker 4.1.6
- Note: Masks the repeats identified through RepeatModeler2
- Help pages accessed: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.html#gsc.tab=0
- https://www.biostars.org/p/9495313/ 
- *Script*: repeatmask.sbatch
  
```
module load repeatmasker/4.1.0
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
RepeatMasker -pa 16 -gff -lib sierrensis-families.fa asierrensis.scaffolded.fasta
```

- Note the above generates files including:
- *“asierrensis.scaffolded.fasta.masked”*, a fasta file like the original one but in which all bases of the detected repeats are replaced (masked) with the letter “N”  
- “asierrensis.scaffolded.fasta.tbl” that contains a brief of the detected repeats, providing details of the classified families and some statistics
- “asierrensis.scaffolded.fasta.out.gff” that provides the start and end position of each detected repeat

### Step 3. Predict genes
- *Software used:* braker2/2.1.6
- Note: uses repeat masked fasta created in step 2
- Can be run with RNA-Seq and protein data if available (N/A here)
- Requirements:
- Augustus (used version: augustus/3.4.0). However, BRAKER requires the config file path to be writeable. To do this, had to install augustus into my SCG account (by installing on my local computer, then uploading)
- geneMark key (i.e., .gm_key) in your home directory. Downloaded GeneMark-ES/ET/EP+ ver 4.72_lic  LINUX 64 kernel 2.6 - 4, 64 bit from here: http://exon.gatech.edu/GeneMark/license_download.cgi. Then uploaded to ref_genome folder as "gm_key"
- following guidance here: https://github.com/Gaius-Augustus/BRAKER?tab=readme-ov-file#running-braker and here: https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_to_Braker2.html#gsc.tab=0
  
```
module load braker2/2.1.6
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
braker.pl --genome=asierrensis.scaffolded.fasta.masked --esmode --min_contig=10000 --AUGUSTUS_CONFIG_PATH=/labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/Augustus/config --GENEMARK_PATH=/labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/gmes_linux_64
```

#### Resulting output 
The final output of the BRAKER gene annotation run is a file 'braker.gtf' which contains the final predicted gene set. (i.e., for each gene: start position, stop position, whether a start_codon, CDS, intron, gene, or mRNA)
For our data set: 765,124 genes identified
braker.gtf file saved here: https://drive.google.com/drive/u/1/folders/1WdF4YUSKHwXKFdPrXm_8ap6yZSsxM7hD

## Assigning candidate SNPs to genes
Note: this was conducted in R using the script ThermalSelectionExp_IdentifyingGenesFromSNPs.R uploaded here 

In brief, this script: 
- pulls in the candidate SNP lists
- defines a function that identifies the gene to which each SNP belongs based on genomic positions
- outputs the resulting gene ID
- calculates the overlap in genes identified through the 3 approaches (Fst on Treatment, GWA on Treatment, and GWA on Knockdown time)
- calculates the expected overlap through permutation 
