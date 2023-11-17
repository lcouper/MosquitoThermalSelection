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

### 2a. Detect using PCadapt
*Note: not doing this approach since most variation not explained in first ~ 10 axes

### 2b. Detect outlier SNPs using Fst
*Note: Fst estimates and outlier detection done using the OutFLANK R package**
This method detected 1,850 outlier SNPs based on Fst, at a q < 0.01 threshold.

### Manhattan plot showing SNP outliers based on Fst 
*Note: outlier SNPs colored below in red. Outlier threshold: q < 0.01*

![FstOutliers_OutFLANK](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a275db80-cdf9-4294-a358-f1bbc8394aa9)


### 2c. Use glm to identify SNPs with significant treatment effects
*i.e.* 
```
model = glm(y ~ df$Treatment, family = "binomial") # where y = (0,2), (1,1), or (2,0)
```
*Note: Typical p-value adjustment methods were overly stringent (i.e. 0 SNPs were maintained at p < 0.05 after correction, given number of tests. Instead, used un-adjusted p-values and considered alpha of 0.0001 when identifying outlier SNPs (shown in red on plot below)*

![Rplot06](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/79fc35eb-0c04-4693-9dfc-e9431adf0c07)

**This method identified 73 outlier SNPs. These appear to cluster around similar areas of the genome as with the Fst approach, but No SNPs were identified as significant/outliers via both methods**


## 3. Examine allele frequency differences between control and heat-selected group
First, calculated frequency of alternative allele for each SNP separately for control and heat-selected.  
Next, calculated difference in allele frequency between these groups.   

### Visualize difference in allele frequency for SNPs identified as outliers through Fst or glm approaches above
![FstGlmOutlierSNPs](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/47318857-8fd3-4ef7-a9a5-7d14d16749de)

### 3b. For each outlier SNP, compare AF differences relative to matched controls
To generate matched controls: 
1. Identify SNPs within 1000 bp of focal SNP
2. Keep those within 5% of starting (control) frequency as focal SNP
3. From this list, take random sample of 7 SNPs
4. Caculate af differences (i.e. between control and heat-treated group) for these 7 SNPs
5. Compare af differences between focal SNP and matched controls

#### Visualize difference in AF for focal SNPs relative to matched controls
Here, pointed indicate focal SNPs with larger (red) or smaller (blue) AF differences than their matched controls

![AFrelativeDiff](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/b24caa4f-820b-438c-bc3f-8d859dbcd8a3)

