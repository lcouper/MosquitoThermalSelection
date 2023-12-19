# Downstream sequence analysis and results on full dataset

## 1. PCA showing overall genomic differences between treatment groups

#### 1a. Scree plot 
Shows the cumulative proportion of variance explained by PC axes   
**Note the majority of variation is NOT explained in first 10 axes**

![ScreePlot_Final](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/d5873555-3542-4799-93de-c0e27663c020)


#### 1a. Control vs heat-selected group
##### Principal components 1 & 2
![PCA_ControlHeat](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/249a4008-684d-4bef-8d6b-00366d133ec5)

##### Principal components 3&4
![PCA_3 4](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/3b2f71a5-4a5d-4349-a994-eb2944c2bd8b)

#### 1b. Female vs male mosquitoes
![PCA_MF](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/c3592e83-4801-41f0-966b-7472ea8a1de4)

#### 1c. By experimental round 1-3
![PCA_ExpRound](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/9cedeea1-8858-4683-b982-06a3ba0868b4)

## 2. Detect outlier SNPs using various methods

### Method 1: PCadapt
*Note: not doing this approach since most variation not explained in first ~ 10 axes

### Method 2: Detect outlier SNPs using Fst
*Note: Fst estimates and outlier detection done using the OutFLANK R package**
This method detected 351 outlier SNPs at q < 0.05 (and Fst >= 0.05).

#### Manhattan plot showing SNP outliers based on Fst 
*Note: outlier SNPs colored below in red*

![Fst_Manhattan](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/5c5f696b-056e-4c7a-858c-7566efaf4b98)


#### Examine allele frequency differences between control and heat-selected group
Step 1: calculated frequency of alternative allele for each SNP separately for control and heat-selected.  
Step 2: calculated difference in allele frequency between these groups.   
Step 3: generated matched controls for each SNP:
To generate matched controls: 
- Identify SNPs on same chromsome and at least 100k bp away from focal SNP
- Keep those within 5% of starting (control) frequency as focal SNP
- From this list, take random sample of 10 SNPs
- Caculate af differences (i.e. between control and heat-treated group) for these 10 SNPs
- Compare af differences between focal SNP and matched controls   

#### Visualize difference in AF for focal SNPs relative to matched controls


### Method 3: GLM using plink 

*Note using: plink v 2.0 through SCG



