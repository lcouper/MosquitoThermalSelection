# Downstream sequence analysis and results on full dataset

## 1. PCA showing overall genomic differences between treatment groups

#### 1a. Scree plot 
Shows the cumulative proportion of variance explained by PC axes   
**Note the majority of variation is NOT explained in first 10 axes**

![Rplot01](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/1b7bca27-35ab-4e21-9b55-7da657531834)

#### 1a. Control vs heat-selected group
![Rplot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/43c324de-412e-48a6-9365-4482694b3d96)

#### 1b. Female vs male mosquitoes
![Rplot02](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a23e1fdc-d3b3-4930-8f3e-6b4100d71c52)

#### 1c. By experimental round 1-3
![Rplot03](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/00fc1b11-f8ae-4001-8cef-564fa6aace2c)

## 2. Detect outlier SNPs using Fst
*Note: Fst estimates and outlier detection done using the OutFLANK R package**
This method detected 56 outlier SNPs based on Fst, at a q < 0.01 threshold.

### Manhattan plot showing SNP outliers based on Fst 
*Note: outlier SNPs colored below in red. Outlier threshold: q < 0.01*
![ManhattanPlot_FstOutlier](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/9c4af6a7-7f79-45f4-b986-4f3725edcce6)

## 3. Examine allele frequency differences between control and heat-selected group
First, calculated frequency of alternative allele for each SNP separately for control and heat-selected.
Next, calculated difference in allele frequency between these groups 
Identified SNPs with allele frequency differences in 99.99th percentile (n = 119)

### SNPs with greatest *increases* (red) or *decreases* (blue) in frequency between control and heat-selected
![af_diff](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/f3bb94ca-29a6-4f30-955e-728b61d07837)

  
