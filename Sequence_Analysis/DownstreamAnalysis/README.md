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
```
Identify outliers
```
outliers = which(padjusted < 0.01)
```
Here, identified 6 significant SNPs at alpha = 0.01
```
colnames(datasub)[outliers]
[1] "V1038779" "V167436"  "V457915"  "V1096118" "V69628"   "V460604"
```

## 3. Use Fst to detect outliers

Note: This includes all SNPs (not just the subset)

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

