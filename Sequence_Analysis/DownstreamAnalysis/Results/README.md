# Downstream sequence analysis and results on full dataset

## 1. PCA showing overall genomic differences between treatment groups

#### 1a. Scree plot 
Shows the cumulative proportion of variance explained by PC axes   
**Note the majority of variation is NOT explained in first 10 axes**

![Rplot01](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/1b7bca27-35ab-4e21-9b55-7da657531834)

#### 1a. Control vs heat-selected group
##### Principal components 1 & 2
![Rplot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/43c324de-412e-48a6-9365-4482694b3d96)

##### Principal components 3&4
![PC3and4](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/f60182b2-0530-496e-90fc-2ee3cb6a2d89)

#### 1b. Female vs male mosquitoes
![Rplot02](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a23e1fdc-d3b3-4930-8f3e-6b4100d71c52)

#### 1c. By experimental round 1-3
![Rplot03](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/00fc1b11-f8ae-4001-8cef-564fa6aace2c)

## 2. Detect outlier SNPs using various methods

### 2a. Detect using PCadapt
*Note: not doing this approach since most variation not explained in first ~ 10 axes

### 2b. Detect outlier SNPs using Fst
*Note: Fst estimates and outlier detection done using the OutFLANK R package**
This method detected 56 outlier SNPs based on Fst, at a q < 0.01 threshold.

### Manhattan plot showing SNP outliers based on Fst 
*Note: outlier SNPs colored below in red. Outlier threshold: q < 0.01*

![Rplot07](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/7b5b6b0c-96c5-4e83-a052-e94d6a39f8ee)


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

