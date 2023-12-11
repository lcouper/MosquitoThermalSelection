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
This method detected 11,446 outlier SNPs based on Fst, a q <0.05 threshold (and 1,850 with q < 0.01).

### Manhattan plot showing SNP outliers based on Fst 
*Note: outlier SNPs colored below in red. Outlier threshold: q < 0.01*

![OutFLANK_Fst_Outliers_q0 05andFst0 05](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/745f9d15-8cb2-4a39-8745-94a91ac3a3e7)

## 3. Examine allele frequency differences between control and heat-selected group
Step 1: calculated frequency of alternative allele for each SNP separately for control and heat-selected.  
Step 2: calculated difference in allele frequency between these groups.   
Step 3: generated matched controls for each SNP:
To generate matched controls: 
- Identify SNPs within 1000 bp of focal SNP
- Keep those within 5% of starting (control) frequency as focal SNP
- From this list, take random sample of 7 SNPs
- Caculate af differences (i.e. between control and heat-treated group) for these 7 SNPs
- Compare af differences between focal SNP and matched controls   

Step 4: Compare AF differences relative to matched controls

#### Visualize difference in AF for focal SNPs relative to matched controls

###### Plot 1. Here, points indicate focal SNPs with larger (red) or smaller (blue) AF differences than their matched controls
![AFdifference_relativetoMatchedControls](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/405c2641-eebf-45d7-ab5d-4554c81e8b8a)

###### Plot 2. Density plot of AF shfits in focal SNPs vs their matched controls
![DensityPlotOfAFshift](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/0cd1aaaf-bd94-461c-8634-6358e89b06cb)

###### Plot 3. Visualize AF shifts in SNPs with largest (95th quantile) differences relative to their matched control
![AFdifferences_95thquartile](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/1830ad8d-6b6f-4dd9-b52b-8a01704a174d)



