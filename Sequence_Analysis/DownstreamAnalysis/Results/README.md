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
Note: then comparing these observed p-values to distribution of permuated pvals

### 2d. Use glm to identify SNPs significantly associated with knockdown time
*i.e.* 
```
model = glm(y ~ df$KDtime + df$Sex, family = "binomial") # where y = (0,2), (1,1), or (2,0)
```
Note: sex added to model as there were clear differences between Males and Females (e.g., due to body size)


