# Downstream sequence analysis

Below, we outline the downstream analysis conducted on the sequence data, after LD-pruning

1. Visualize overall genomic differences between treatments  using PCA
2.      treatments: heat-selected vs control, male vs female, experimental round
3.  Detect candidate loci using Fst (calculated in OutFlank) between control and heat-selected groups
4.  Detect candidate loci using GWA (conducted in plink) with treatment as the phenotype
5.  Detect candidate loci using GWA (conducted in plink) with knckdown time as the phenotype
6.  Examine allele frequency differences in candidate loci identified above

