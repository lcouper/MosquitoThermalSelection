# Downstream sequence analysis
*Pertains to analysis performed in R, after genotype matrix was created*

## 1. Visualize overall differences between treatment groups using PCA 

```
datasub = read.csv("GenotypeMatrix_Subset.csv”)[-1]
# center and scale genotype matrix
datasubM = as.matrix(datasub)
datasubMS = scale(datasubM, center = TRUE, scale = TRUE)
# get data into pcadapt-readable format
pcadata = read.pcadapt(t(datasubMS)) 
screetest <- pcadapt(input=pcadata,K=20) # K = # of principal components to retain
plot(screetest,option="screeplot") # Looks like ~5 is the correct number to retian
x = pcadapt(input= pcadata, K=5)

# Import file with metadata
data = read.csv("FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09")
data = data[!(data$Sample %in% Removes),] 
Treatment = data$Treatment
Sex = data$Sex
ExpRound = data$Exp.Round

```

![Screeplot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/ec6cc449-75e7-40a5-8d54-63c513fef785n)

  

### PCA: Heat-selected vs control 
```
plot(x,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot139](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/8dc59047-a108-4873-bb05-2e33ae073723)


### PCA: Female vs Male
```
plot(x,option="scores", pop=Sex, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot140](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/f84fdfe1-a313-4c5d-933f-b92edc9b49ff)

### PCA: Experimental Round (1-3)
```
plot(x,option="scores", pop=ExpRound, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot141](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/877eb0e6-991d-4afa-a734-83a2625859f7)

## 2. Detect outlier SNPs

Manhattan plot to visualize potentially important SNPs 

```
plot(x,option="manhattan") + theme_bw()
```
![Rplot142](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a71373a3-812e-4a26-96af-031a80046e56)


Adjust p-values for multiple testing
```
padjusted = p.adjust(x$pvalues, method = "fdr")
outliers = which(padjusted < 0.01) 
# Here, identified 6 significant SNPs at alpha = 0.01
colnames(datasub)[outliers]
[1] "V1038779" "V167436"  "V457915"  "V1096118" "V69628"   "V460604"
```

## 3. Ideally: Use Fst to detect outliers

*The OutFLANK package appears to be the common method for estimating Fst for outlier detection
However it requires the q-value package, which I can not seem to install for the life of me
Will try doing this on new laptop next week*

## 4. LFMM anaysis to detect SNPs associated with longer knockdown times 

Goal = identify SNPs associated with longer knockdown time (i.e., greater heat tolerance) *within* a given treatment, while also controlling for variation in sex and body size
i.e.:

![Screen Shot 2023-08-14 at 3 31 48 PM](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/8a296d4a-e41a-475e-ae14-84dae1df4bd8)

*I'm not sure of the proper way to do this. My attempt so far is below (but I don’t think this is correct)*

```
# LFMM analysis to detect SNPs associated with longer KD time

# tutorial from here: https://cran.r-project.org/web/packages/lfmm/vignettes/lfmm.html
library(lfmm)
library(scales)
library(dplyr)

# create genotype and phenotype matrices
Y2 = datasub
X2kd_unscaled = data$Kdtime
X2kd = as.numeric(scale(X2kd_unscaled, center = TRUE, scale = TRUE))
X2treat = as.numeric(recode(data$Treatment,  "control" = '1',"high temp" = '2'))
X2sex = as.numeric(recode(data$Sex,  "M" = '1',"F" = '2'))
X2 = cbind.data.frame(X2kd, X2treat, X2sex)

# Set up LFMM model
mod.lfmm <- lfmm_ridge(Y = Y2,  X = X2, K = 5) 

# Identify significant SNPs (Unsure if I'm doing this correctly)
pv <- lfmm_test(Y = Y2, X = X2,  lfmm = mod.lfmm, calibrate = "gif")
# Pulls out p-values associated for SNP - KDtime associations
pvaluesKD <- pv[["calibrated.pvalue"]][,1]
outliers2 = which(pvaluesKD < 0.01)
# Here, identified 2 SNPs:
colnames(datasub)[outliers2]
# [1] "V3439"   "V794831"
```
Next, I would identify if there are any overlaps in SNPs identified by LFMM and the outlier SNPs between treatment and control (no overlap so far, using data subset)

Manhattan-ish plot, visualizing significant SNPs in this model
```
plot(-log10(pvalues), pch = 19, cex = .8, 
     xlab = "SNP", ylab = "-Log P", col = "grey")
```
![Rplot143](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/cf82afcd-9ac2-4406-a7e2-250ee6e83ccf)





