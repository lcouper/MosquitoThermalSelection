# Downstream sequence analysis

**The following downstream analyses were conducted on the sequence data, after LD-pruning**

1. Visualize overall genomic differences between treatments (heat-selected vs control, male vs female, experimental round) using PCA
2. Detect candidate loci using Fst (calculated in OutFlank) between control and heat-selected groups
3. Detect candidate loci using GWA (conducted in plink) with treatment as the phenotype
4. Detect candidate loci using GWA (conducted in plink) with knckdown time as the phenotype
5. Examine allele frequency differences in candidate loci identified above
6. Estimate Tajima's D across all loci 

**These steps can be replicated on a subset of the data 'GenotypeMatrix_Pruned_Subset' (as the full pruned genotype matrix is too large to upload) using the script: 'ThermalSelection_SeqAnalysis_LD_Pruned.R'**
