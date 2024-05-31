# Steps for gene function prediction 
Gene function predictions were obtained through NCBI nucleotide BLAST search.
Gene sequences for input into BLAST were obtained based on gene start and stop positions identified through Braker (see *DownstreamAnalysis_GeneAnnotation*) and are listed in the file: *GeneList_AllApproaches.csv*.  
That is, using the start and stop positions for each gene: 
```
samtools faidx asierrensis.scaffolded.fasta 1_RagTag:372060-374060
```

<img width="1182" alt="Screenshot 2023-12-22 at 10 46 28 AM" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/31b649fb-dd54-41b3-b44b-c89f7511a225">
