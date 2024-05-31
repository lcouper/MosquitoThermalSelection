# Steps for gene function prediction 
Gene function predictions were obtained through NCBI nucleotide BLAST search.
Gene sequences for input into BLAST were obtained based on gene start and stop positions identified through Braker (see *DownstreamAnalysis_GeneAnnotation*) and are listed in the file: *GeneList_AllApproaches.csv*.  
That is, using the start and stop positions for each gene: 
```
samtools faidx asierrensis.scaffolded.fasta 1_RagTag:372060-374060
```
where 1_RagTag denotes the chromosome, and the following numbers denote the genomic position
