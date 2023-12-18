# Downstream sequence analysis conducted in Plink
*Pertains to analysis performed in plink  


#### Step 1. Create files for use in plink 
Followed guidance from: https://www.biostars.org/p/109690/ 
Another useful tutorial here: http://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/

As our genome involves many scaffolds, we first need to create our own chromosome mapping file. Do so, by running:
```
bcftools view -H Filtered_VCF_All_sorted_0.995_bialleliconly.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > filename.chrom-map.txt
```
Then, make mapping file by running:
```
module load vcftools/0.1.16-13-gd0c95c5
# Note that earlier versions of vcftools on SCG do not have the --chrom-map option, so must load this version specifically
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --plink --chrom-map filename.chrom-map_maf01.txt --maf 0.1 --out myplink
```

#### Step 2. Create LD-pruned data set in plink 
Followed tutorial here: https://zzz.bwh.harvard.edu/plink/tutorial.shtml
Note: running plink v 1.9 on SCG 

```
# Navigate to working directory with .map and .ped files (here, /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF)
# to check that things are working:
plink --allow-extra-chr --file myplink

# Create LD-pruned dataset
plink --allow-extra-chr --file myplink --indep 50 5 1.5
# Note this retains 609,490 SNPs

# to get list of linked and focal SNPs:
plink --allow-extra-chr --file myplink_maf1 --indep 50 5 1.5 --show-tags plink.prune.in
```

# Make new, pruned file and create binary bed file
# note the --recode argument is necessary to generate the .map and .ped files needed for the association analysis
plink --allow-extra-chr --file myplink --extract plink.prune.in --recode --make-bed --out pruneddata
```

#### Step 3. Conduct GWA with treatment as phenotype 
Here, treat.phe.txt (uploaded here) is a 3 column file specificying the family ID and individual ID (here the same thing) and the phenotype-- here '1' for control and '2' for heat-selected. Note that the --allow-no-sex flag is mandatory for this line to run.
The max(T) permutation approach is discussed further here: https://zzz.bwh.harvard.edu/plink/perm.shtml

```
# Conduct with max(T) permutation approach (as a means of obtaining corrected p-values)
plink --allow-extra-chr --file pruneddata --pheno treat.txt --allow-no-sex --assoc --mperm 5000 --out treat_assoc
```

Conduct association, but clump results to account for LD
```
plink --allow-extra-chr --file myplink --pheno treat.txt --allow-no-sex --clump pruned_treat.assoc
# Note that the 'myplink' data file specified above is used to calculate LD between SNPs in the .assoc file
# In our case, it is equivalent to running: 
plink --allow-extra-chr --file pruneddata --pheno treat.txt --allow-no-sex --clump pruned_treat.assoc
```

#### Step 4. Conduct GWA with knockdown time as phenotype 
Here, KD.phe.txt (uploaded here) contains the individual knockdown times. As above, we use a permutation approach to obtain corrected significance values for each SNP. Here, we specify that permutatons should occur within-sex clusters. This is to account for known differences in body size (and potentially heat tolerance) between adult female and male mosquitoes. KD.sex.cluster.txt note the sex of all individuals, with 1 = F, 2 = M. 
```
plink --allow-extra-chr --file myplink3 --pheno KD.phe.txt --allow-no-sex --assoc --mperm 5000 --within KD.sex.cluster.txt --out KD_assoc
```
