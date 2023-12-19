![image](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/c7a285dc-3249-4c32-8a3c-22521a08fc6a)# Downstream sequence analysis and results on full dataset

## 1. PCA showing overall genomic differences between treatment groups

####  Scree plot 
Shows the cumulative proportion of variance explained by PC axes   
**Note the majority of variation is NOT explained in first 10 axes**

![ScreePlot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/dd32a321-5cef-46e2-9157-514b2fd98d47)

#### Control vs heat-selected group 
##### Principal components 1 & 2
![PCA_12](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/8b113734-0c96-49cb-8111-e5de4d15818f)

##### Principal components 3&4
![PCA_34](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/f6d5a685-2846-4904-9d04-02e4a46af17e)

#### Female vs male mosquitoes
![PCA_MF](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/53bf07c0-4626-4a25-8b93-44a4d35ca248)

#### By experimental round 1-3
![PCA_ExpRound](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/4f0750a7-949a-4028-9f24-f543e943f5b2)

## 2. Detect candidate loci using Fst
*Note: Fst estimates and outlier detection done using the OutFLANK R package**
This method detected 351 outlier SNPs at q < 0.05 (and Fst >= 0.05).

**Manhattan plot showing candidate loci (in bold) across the 3 chromosomes based on Fst**
![Fst_Manhattan](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/5c5f696b-056e-4c7a-858c-7566efaf4b98)

## 3. Detect candidate loci using GWA with treatment as the phenotype
*Note: analysis conducted in plink. Results were 'clumped' based on LD*

**Manhattan plot showing candidate loci (in bold) across the 3 chromosomes based on Fst**
![Rplot04_Treat](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/e0a2490a-6a0b-4b24-afa5-b6ed1066d07f)

## 4. Detect candidate loci using GWA with knockdown time as the phenotype
*Note: analysis conducted in plink. Results were 'clumped' based on LD*

**Manhattan plot showing candidate loci (in bold) across the 3 chromosomes based on Fst**
![KD](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/54280dd8-8c61-4be6-ab04-72365d458b4c)

## 5. Examine allele frequency differences for candidate loci 
Methods:
1. Calculate frequency of alternative allele for each SNP separately for control and heat-selected.  
2. Calculate difference in allele frequency between these groups.   
3. Generated matched controls for each SNP:
To generate matched controls: 
- Identify SNPs on same chromsome and at least 100k bp away from focal SNP
- Keep those within 5% of starting (control) frequency as focal SNP
- From this list, take random sample of 10 SNPs
- Caculate af differences (i.e. between control and heat-treated group) for these 10 SNPs
- Compare af differences between focal SNP and matched controls   

**Difference in AF for focal SNPs relative to matched controls**
*Note: here, candidate loci are those identified through methods 2 & 3 above*  
*green = focal SNP, gray = matched controls 1-10*
![AF_shifts](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/33c21553-2e17-4fe3-bdd5-f8ae9343939a)

K-S test indiciate the distribution of allele frequency differences in focal SNPs differs from that of matched controls.  
But distributions of all matched controls do not differ

## 6. Estimate Tajima's D across loci
![TajimaD](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/50789688-81c9-4a91-84e7-38c1630cd577)

