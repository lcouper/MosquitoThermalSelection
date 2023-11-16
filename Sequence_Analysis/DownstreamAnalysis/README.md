# Downstream sequence analysis
*Pertains to analysis performed in R, after genotype matrix was created*  
*The analysis/figures shown here pertain to analysis of a **100 SNP subset** of the genotype matrix*

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
pcobj = pcadapt(input= pcadata, K=5)

# Import file with metadata
data = read.csv("FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09")
data = data[!(data$Sample %in% Removes),] 
Treatment = data$Treatment
Sex = data$Sex
ExpRound = data$Exp.Round

```
![ScreePlot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/fa113cf0-f92a-4c4b-a507-d3df67882de3)

### PCA: Heat-selected vs control 
```
plot(pcobj,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
![PCA_ControlTreat_Subset](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/5292cc8c-7920-45ab-909b-ef5009d1fdcc)

### PCA: Female vs Male
```
plot(pcobj,option="scores", pop=Sex, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
![PCA_MF_Subset](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/127393b8-ae36-4f80-bd0e-96be28eca1da)


### PCA: Experimental Round (1-3)
```
plot(pcobj,option="scores", pop=ExpRound, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
![PCA_ExpRound_Subset](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/c45bf7f7-0344-46e5-ae0a-722123421d7b)


## 2. Detect outlier SNPs using various methods

### 2a. Detect using PCadapt

Manhattan plot to visualize potentially important SNPs 
```
plot(x,option="manhattan") + theme_bw()
```
![Manhattan_Subset](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/39fd363b-7eca-4062-b621-c3361ab20d33)

Adjust p-values for multiple testing
```
padjusted = p.adjust(pcobj$pvalues, method = "fdr")
```
Identify outliers
```
outliers = which(padjusted < 0.01)
```
Here, identified 17 significant SNPs at alpha = 0.01
```
colnames(datasub)[outliers]
[1] "X3_RagTag_317222084" "X2_RagTag_285266406" "X2_RagTag_358523747" "X2_RagTag_254025183" "X1_RagTag_226675236" "X3_RagTag_152027826" "X2_RagTag_283502580" "X1_RagTag_49700296"  "X3_RagTag_108965707" "X2_RagTag_295529825" "X3_RagTag_94856485"  "X3_RagTag_157997951" "X2_RagTag_1920050"   "X2_RagTag_21477677"  "X2_RagTag_309792817" "X3_RagTag_140211440" "X3_RagTag_98148127" 
```

### 2b. Detect using Fst and OutFlankR

```
fst = MakeDiploidFSTMat(datasub, locusNames = 1:ncol(datasub), popNames = data$Treatment)
hist(fst$FST, breaks = 50)
```
![Fstvaluedist_Method2](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/b90e66fc-4850-4aca-a605-34a5270ba46a)

Identify statistical outliers based on chi-squared distribution
```
out1 <- OutFLANK(fst,NumberOfSamples = 2) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Note: with this datasubset, pvalues are not uniformly distributed
# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                       FALSE, RightZoomFraction = 0.05, titletext = NULL)
```
![Fst_NoSampleSizeCorrection](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/ec717cc8-aafa-44cc-984b-4877826f9c52)

```
# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=OF$FSTNoCorrbar,
      dfInferred=OF$dfInferred,qthreshold=0.05,Hmin=0.1)

# Which have q-values < 0.05
outliers <- which(P1$OutlierFlag==TRUE)
colnames(datasub)[outliers] # 1 outlier SNP in the data subset ("X3_RagTag_280034261")

# Create manhattan plot with outlier SNPs
plot(P1$LocusName,P1$FST,xlab="SNP Position",ylab="FST",col=rgb(0,0,0,alpha=0.3), 
     pch = 16)
points(P1$LocusName[outliers],P1$FST[outliers],col="red", pch = 16)
```
![Manhattan_FstOutliers](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/3d55f889-9c6d-4ccc-8e6f-fc02aa75ba0d)




```

