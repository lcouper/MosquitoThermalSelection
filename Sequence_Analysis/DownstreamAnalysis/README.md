# Downstream sequence analysis

**The following downstream analyses were conducted in R on the sequence data, after LD-pruning. These steps can be replicated on a subset of the data** ('GenotypeMatrix_Pruned_Subset' **using the script** 'ThermalSelection_SeqAnalysis_LD_Pruned.R'

The downstream analysis steps include:
1. Visualizing overall genomic differences between treatments (heat-selected vs control, male vs female, experimental round) using PCA
2. Detecting candidate SNPs using Fst (calculated in OutFlank) between control and heat-selected groups
3. Detecting candidate SNPs using GWA (conducted in plink) with treatment as the phenotype
4. Detecting candidate SNPs using GWA (conducted in plink) with knckdown time as the phenotype
5. Examining allele frequency differences in candidate SNPs identified above
6. Generating matched controls for each SNPs and examine allele frequency differences in these
7. Calculating and plotting selection coefficients
8. Plotting population diversity metrics

