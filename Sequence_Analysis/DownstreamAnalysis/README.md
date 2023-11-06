# Downstream sequence analysis
*Pertains to analysis performed in R, after genotype matrix was created*  
*The analysis/figures shown here pertain to analysis of a **120 SNP subset** of the genotype matrix*

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
![PC_controlHightemp](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/2400bd1b-fd4e-4406-8d09-398543ed5489)

### PCA: Female vs Male
```
plot(pcobj,option="scores", pop=Sex, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
![PC_malefemale](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/dc5818f5-fe25-4ac1-9a8e-5a27d9daac54)

### PCA: Experimental Round (1-3)
```
plot(pcobj,option="scores", pop=ExpRound, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
![PC_expround](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/14771264-4411-4865-a168-da3fc019bcfe)

## 2. Detect outlier SNPs

Manhattan plot to visualize potentially important SNPs 

```
plot(x,option="manhattan") + theme_bw()
```
![manhattan](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/a4de6047-7f4d-4e9b-9a55-f641087909cb)


Adjust p-values for multiple testing
```
padjusted = p.adjust(pcobj$pvalues, method = "fdr")
```
Identify outliers
```
outliers = which(padjusted < 0.01)
```
Here, identified 8 significant SNPs at alpha = 0.01
```
colnames(datasub)[outliers]
[1] "V622304"  "V2681962" "V1942801" "V1238265" "V2338596" "V2386972" "V3379779" "V3220287"
```

## 3. Use Fst to detect outliers

Note: Method 1 currently includes all SNPs, while method 2 uses only the 100 SNPs in the data subset so they are not directly comparable. Perhaps use both approaches in final analysis and retain SNPs identified through both methods?

### Method 1: Using Fst values generated through vcftools 

See step 23 of Seqeuence Analysis for vcftools Fst call

```
Fstdata = read.delim("Fst_estimates_VCFtools/Fst_estimates_controlvsheat_vcftools.txt")
colnames(Fstdata) = c("Chrom", "SNP", "WeirFst")
hist(Fstdata$WeirFst)
```
![Fstvaldist](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/9a670d8d-f6f8-4f8b-a4ae-40b644f8bed2)

Detect outliers as those exceeding 99.9th% percentile.  
Note: following methods here: https://speciationgenomics.github.io/per_site_Fst/

```
my_threshold <- quantile(Fstdata$WeirFst, 0.999, na.rm = T)
fst <- Fstdata %>% mutate(outlier = ifelse(Fstdata$WeirFst > my_threshold, "outlier", "background"))
length(which(fst$outlier == "outlier")) # identifies 1204 SNPs
```
![Fstplot_man](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/bfe65f29-6724-44a1-80f1-3bc636e654cc)

### Method 2: Using OutFLANK R package 

```
fst = MakeDiploidFSTMat(datasub, locusNames = 1:ncol(datasub), popNames = data$Treatment)
hist(fst$FST, breaks = 50)
```

![Fstdist_method2](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/2327fa5f-fc2c-43bc-96b2-8e8de6a49d60)

Identify statistical outliers based on chi-squared distribution
```
out1 <- OutFLANK(fst,NumberOfSamples = 2) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Note: with this datasubset, pvalues are not uniformly distributed
# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                       FALSE, RightZoomFraction = 0.05, titletext = NULL)
```
![Rplot_fst](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/8b85d612-3f5c-4fa1-b2eb-bf51dfba6b5f)

```
# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=OF$FSTNoCorrbar,
      dfInferred=OF$dfInferred,qthreshold=0.05,Hmin=0.1)

# Which have q-values < 0.05
outliers <- which(P1$OutlierFlag==TRUE)
colnames(datasub)[outliers] # identified one SNP: V200615

# Create manhattan plot with outlier SNPs
plot(P1$LocusName,P1$FST,xlab="SNP Position",ylab="FST",col=rgb(0,0,0,alpha=0.3), 
     pch = 16)
points(P1$LocusName[outliers],P1$FST[outliers],col="red", pch = 16)
```

![Rplot_Fstoutiler](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/ab4e819a-52be-490a-8658-43f3e2962af0)



## 4. Detect SNPs associated with longer knockdown times 

Identifies SNPs associated with longer knockdown time (i.e., greater heat tolerance) while controlling for treatment group and sex
```
# Combine relevant metadata and genotype matrix into dataframe
df = cbind(metadata[,c(8,2,4)], datasub)
df$Treatment = as.factor(df$Treatment)
df$Sex = as.factor(df$Sex)

# Run model KD time ~ genotype matrix + Treatment + Sex (not using body size for now) 
model = lm(Kdtime ~ ., data = df)

# Pull out p-values and adjust for multiple testing
pvals = as.numeric(summary(model)$coefficients[,4])
# padjusted = p.adjust(pvals, method = "fdr") # This appears to be over-correcting
outliers = which(pvals < 0.05) 

# Pull out ID of the SNP outliers
colnames(df)[outliers] # in this example: "V790007" "V254278" "V790890" "V248541" "V489649"
```

