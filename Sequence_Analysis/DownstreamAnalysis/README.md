# Downstream sequence analysis

**The following downstream analyses were conducted in R on the sequence data, after LD-pruning**

1. Visualize overall genomic differences between treatments (heat-selected vs control, male vs female, experimental round) using PCA
2. Detect candidate loci using Fst (calculated in OutFlank) between control and heat-selected groups
3. Detect candidate loci using GWA (conducted in plink) with treatment as the phenotype
4. Detect candidate loci using GWA (conducted in plink) with knckdown time as the phenotype
5. Examine allele frequency differences in candidate loci identified above
6. Identify candidate genes located within 100K bp of candidate loci


**These steps can be replicated on a subset of the data 'GenotypeMatrix_Pruned_Subset' (as the full pruned genotype matrix is too large to upload) using the script: 'ThermalSelection_SeqAnalysis_LD_Pruned.R'**

## how to pull out sequence surrounding SNP
samtools faidx asierrensis.scaffolded.fasta 1_RagTag:372060-374060
<img width="1182" alt="Screenshot 2023-12-22 at 10 46 28 AM" src="https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/31b649fb-dd54-41b3-b44b-c89f7511a225">
