# Downstream sequence analysis conducted in Plink
*Pertains to analysis performed in plink.
Plink files stored in: ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF/


### Step 1. Create files for use in plink 
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

### Step 2. Create LD-pruned data set in plink 
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

# Make new, pruned file and create binary bed file
# note the --recode argument is necessary to generate the .map and .ped files needed for the association analysis
plink --allow-extra-chr --file myplink --extract plink.prune.in --recode --make-bed --out pruneddata
```

### Step 3. Conduct GWA with treatment as phenotype (sex as covariate)
Here, treat.phe.txt (uploaded here) is a 3 column file specificying the family ID and individual ID (here the same thing) and the phenotype-- here '1' for control and '2' for heat-selected. Similarly, covars.txt is a 4 column file with 3rd column listing Sex
Note that the --allow-no-sex flag is mandatory for this line to run.
Following guidance here: https://zzz.bwh.harvard.edu/plink/perm.shtml
Note: I'm using post-association clumping rather than the permutation approach to obtain significance values (two methods are not compatible)

```
# Conduct association analysis using Sex as a covariate (the first covariate listed in covars.txt)
plink --allow-extra-chr --file pruneddata --pheno treat.txt --allow-no-sex --covar covars.txt --covar-number 1 --assoc --out pruned_treat_assoc_SexCovar
```

Clump results from GWA to account for LD
```
plink --allow-extra-chr --file pruneddata --pheno treat.txt --allow-no-sex --clump pruned_treat_assoc_SexCovar.assoc
# Note that the 'pruneddata' data file specified above is used to calculate LD between SNPs in the .assoc file
```

Pruning and clumping resulted in **112 SNPs** retained as significant (at < 0.01 after Benjamini-Hochberg FDR correction)

### Step 4. Conduct GWA with knockdown time as phenotype (sex and treatment as covariates)
Here, KD.phe.txt (uploaded here) contains the individual knockdown times. covars.txt contains information on the covariates (here, sex and treatment).

```
plink --allow-extra-chr --file pruneddata --pheno KD.phe.txt --allow-no-sex --covar covars.txt --covar-number 1,2 --assoc --out pruned_KD_assoc_Covar

```

Clump results from GWA to account for LD
```
plink --allow-extra-chr --file pruneddata --pheno KD.phe.txt --allow-no-sex --clump pruned_KD_assoc_Covar.qassoc
```
Pruning and clumping resulted in **123 SNPs** retained as significant (at < 0.01 after Benjamini-Hochberg FDR correction)

### Step 5. Estimating variance explained by SNPs using GCTA
### (*all* SNPs in Steps A & B, and *focal SNPS only* in Step C)

Following tutorial here: https://yanglab.westlake.edu.cn/software/gcta/#Tutorial

**Step A**: GCTA-GRM: calculating the genetic relationship matrix (GRM) from all the autosomal SNPs
- Note: all SNPs here are likely autosomal, so use all identified SNPs  
- Note: Given errors with reading "1_RagTag" as chromosome names, prior to running the commande below, I had to alter pruneddata.bim file to replace "1_RagTag" to "1" (same for 2_RagTag and 3_RagTag)
```
gcta64 --bfile pruneddata --autosome --make-grm --out pruneddata --autosome-num 3 --thread-num 12
```

**Step B**: GCTA-GREML analysis: estimating the variance explained by the SNPs
Note: the above command createa a genetic relationship matrix among the autosomal SNPs. This file is then used to estimate variance explained by the SNPs
Note: since the phenotype here is treatment (i.e., being in the control vs heat-selected group), GCTA considers this a case-control analysis 
```
gcta64 --grm pruneddata --pheno treat.txt --reml --out pruneddata2 --thread-num 12 --prevalence 0.457 --covar covars_gcta.txt
# here prevalence of 0.457 refers to baseline larval mortality (i.e, mortality rate in control), and covariates = sex
```

**Step C** Repeating for adult thermal knockdown resistance:
```
gcta64 --grm pruneddata --pheno KD.phe.txt --reml --out pruneddata2 --thread-num 12 --prevalence 0.457 --covar covars.txt
# covariates here are sex and treatment
```

**Step D** Re-run A-B using only focal SNPs
Here, focal SNPs = those identified through Fst & GWA based on treatment.  
Plink files for this named 'pruneddata_subset'
```
gcta64 --bfile pruneddata_subset --autosome --make-grm --out pruneddata --autosome-num 3 --thread-num 12
gcta64 --grm pruneddata_subset --pheno treat.txt --reml --out pruneddata_subset --thread-num 12 --prevalence 0.457 --covar covars_gcta.txt
```




