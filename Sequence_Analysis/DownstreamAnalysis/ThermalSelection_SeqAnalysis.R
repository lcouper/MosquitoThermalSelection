#### Thermal Selection Experiment Sequence Analysis ####

#### Set wd and load libraries ####
library(pcadapt)
library(dplyr)
library(ggplot2)
library(OutFLANK)
library(dplyr)
library(tidyverse)
library(data.table)

setwd("~/Documents/Current Projects/Thermal Selection Experiment")

#### Create data subset (only need to run once)  #####
library(data.table)

setwd("~/Documents/Current Projects/Thermal Selection Experiment")
#x = fread("GenotypeMatrix/GenotypeMatrix", header = FALSE)
header = read.csv("GenotypeMatrix/GenotypeMatrix_SampleNames.csv", header = FALSE)
SampleNames = as.character(header$V1)

SampleInfo = read.csv("FullAlignment&ExperimentStats.csv", header = T)
rownames(x) = SampleNames

# Remove sample 'E-014' (did not have KD time recorded)
# Remove samples 'E-5, E-6, E-7, and E-9' <- were not properly sequenced
x = x[-c(5, 18:21),]

# subset to make analyses easier
set.seed(1)
datasub =  dplyr::select(x, sample(seq_len(ncol(x)), size = 100))
write.csv(datasub, "GenotypeMatrix/GenotypeMatrix_Subset.csv")

#### Perform/QC analysis on data subset
##### Subset: Load data files #####

# Load genotype matrix subset
datasub = read.csv("GenotypeMatrix/GenotypeMatrix_Subset.csv")[,-1]

# Load metadata
metadata = read.csv("FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09")
metadata = metadata[!(metadata$Sample %in% Removes),] 
Treatment = metadata$Treatment
Sex = metadata$Sex
ExpRound = metadata$Exp.Round


##### Subset. 1. Visualize differences using PCA #####

# Step 1: Center and scale genotype matrix
datasubM = as.matrix(datasub)
datasubMS = scale(datasubM, center = TRUE, scale = TRUE)

# Step 2: Get data into pcadapt-readable format. Input = transposed genotype matrix
# (For data subset, will get warning message about having more individuals than SNPs)
pcadata = read.pcadapt(t(datasubMS)) 
screetest <- pcadapt(input=pcadata,K=20) # K = # of principal components to retain
plot(screetest,option="screeplot") # Looks like ~5 is the correct number to retian
x = pcadapt(input= pcadata, K=5)

# Step 3: Plot PCA 
# Control vs heat-selected
plot(x,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
# Female vs Male
plot(x,option="scores", pop=Sex, 
     col = c("#e75480", "#89CFF0")) + theme_bw() + ggtitle(label = NULL)
# By Experimental Round (1-3)
plot(x,option="scores", pop=ExpRound, 
     col = c("#A50026", "#F99858", "#364B9A")) + theme_bw() + ggtitle(label = NULL)


##### Subset. 2. Detect outlier SNPs between heat and control group #####
##### Subset. 2a. Method 1: pcadapt #### 

# Note this uses the df 'x' created in step 1

# Step 1: Create Manhattan plot to visualize potentially important SNPs
plot(x,option="manhattan") + theme_bw() + ggtitle(label = NULL) 

# Step 2. Adjust p-values for multiple testing
# To look at histogram of p-values 
#hist(x$pvalues) # Slight U-shape <- is this a problem?
padjusted = p.adjust(x$pvalues, method = "fdr")
outliers = which(padjusted < 0.01) 
# Here, identified 6 significant SNPs at alpha = 0.01
colnames(datasub)[outliers]


##### Subset. 2b. Method 2: Use Fst to find outliers #####

# Following methods here: http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html

fst = MakeDiploidFSTMat(datasub, locusNames = 1:ncol(datasub), popNames = data$Treatment)
hist(fst$FST, breaks = 50, xlab = "Fst", ylab = "", main = "Fst value distribution: Method 2")

out1 <- OutFLANK(fst,NumberOfSamples = 2) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Not uniformly distributed

# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                       FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=out1$FSTNoCorrbar,
                                dfInferred=out1$dfInferred,qthreshold=0.05,Hmin=0.1)
# Which have q-values < 0.05
outliers <- which(P1$OutlierFlag==TRUE)
colnames(datasub)[outliers] # identified one SNP: V200615

# Manhattan plot with outlier SNPs
plot(P1$LocusName,P1$FST,xlab="SNP Position",ylab="FST",col=rgb(0,0,0,alpha=0.3), 
     pch = 16)
points(P1$LocusName[outliers],P1$FST[outliers],col="red", pch = 16)


##### Subset. 2c. Method 3: Model af differences ####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfaf = cbind(metadata$Treatment, datasub)
colnames(dfaf)[1] = "Treatment"
dfaf$Treatment = as.factor(dfaf$Treatment)
rownames(dfaf) = metadata$Sample

# Separate control and heat-selected treatments
control = dfaf[dfaf$Treatment == "control",-1]
heat = dfaf[dfaf$Treatment == "high temp",-1]

# calculate af
af_control = colSums(control)/ (2*nrow(control))
af_heat = colSums(heat)/ (2*nrow(heat))

# combine and reshape for model
af_table = cbind.data.frame(af_control, af_heat)
colnames(af_table) = c("control", "heat")
af_table$SNP = rownames(af_table)
afl = melt(af_table)
colnames(afl) = c("SNP","treatment", "af")
afmod = glm(af ~ treatment, data = afl, family = quasibinomial(link = "logit"))
summary(afmod)




##### Subset. 2d. Plot AF differences #####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfaf = cbind(metadata$Treatment, datasub)
colnames(dfaf)[1] = "Treatment"
dfaf$Treatment = as.factor(dfaf$Treatment)
rownames(dfaf) = metadata$Sample

# Separate control and heat-selected treatments
control = dfaf[dfaf$Treatment == "control",-1]
heat = dfaf[dfaf$Treatment == "high temp",-1]

# calculate af
af_control = colSums(control)/ (2*nrow(control))
af_heat = colSums(heat)/ (2*nrow(heat))

# combine and calculate differences
af_table = cbind.data.frame(af_control, af_heat)
colnames(af_table) = c("control", "heat")
af_table$diff = af_table$heat - af_table$control

# Identify outliers (i.e. SNPs with significant diff frequencies between control and heat)
threshold = quantile(af_table$diff, 0.99, na.rm = T) # Adjust this later
af_outliers =  af_table %>% mutate(outlier = ifelse(af_table$diff > threshold, "outlier", "background"))
rownames(af_table)[which(af_outliers$outlier == "outlier")] # 1 outlier detected: "V183021"

# Reshape and plot
aft = cbind(rownames(af_table), af_outliers[,-3])
colnames(aft)[1] = c("SNP")
aft2 = aft %>% pivot_longer(-c(SNP, outlier), names_to = "treatment", values_to = "af") 

ggplot(aft2,aes(x=treatment,y=af)) + theme_bw() +
   geom_point(aes(x = treatment, y = af, group = SNP, col = outlier)) + 
  scale_fill_manual(values = c("grey80", "red")) + 
   geom_line(aes(x = treatment, y = af, group = SNP, col = outlier)) +
  scale_color_manual(values = c("grey80", "red"))





##### Subset. 3. Use lm to identify genotypic predictors of KD time #####

# Combine relevant metadata and genotype matrix into dataframe
df = cbind(metadata[,c(8,2,4)], datasub)
df$Treatment = as.factor(df$Treatment)
df$Sex = as.factor(df$Sex)

# Model KD time ~ genotype matrix + Treatment + Sex (not using body size for now) 
model = lm(Kdtime ~ ., data = df)

# Pull out p-values and adjust for multiple testing
pvals = as.numeric(summary(model)$coefficients[,4])
padj = p.adjust(pvals, method = "fdr") 
outliers = which(padj < 0.05) 

# Pull out ID of the SNP outliers
colnames(df)[outliers] # in this example: none

# Second model excluding the Treatment and Sex predictors
df2 = cbind(metadata$Kdtime, datasub)
model2 = lm(df2$`metadata$Kdtime` ~ ., data = df2)
pvals2 = as.numeric(summary(model2)$coefficients[,4])
padj2 = p.adjust(pvals2, method = "fdr") 
outliers2 = which(padj2 < 0.05) 

colnames(df2)[outliers2] # still no significant SNPs identified in subset




#### Conduct above analysis on full dataset ####

#### Full. Load data files ####

datafull = fread("GenotypeMatrix/GenotypeMatrix.txt", header = FALSE)[,-1] # First row is sample number
# Load sample headers and locus names
header = read.csv("GenotypeMatrix/GenotypeMatrix_SampleNames.csv", header = FALSE)
SampleNames = as.character(header$V1)
rownames(datafull) = SampleNames

locus = fread("GenotypeMatrix/GenotypeMatrix_LocusNames.csv", header = FALSE)
locusNames = as.character(locus$V2)
colnames(datafull) = locusNames

# Remove samples without metadata from dataset
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09")
datafull = datafull[!(rownames(datafull) %in% Removes),] 
 
# Load metadata and sample names
metadata = read.csv("FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
metadata = metadata[!(metadata$Sample %in% Removes),] 
Treatment = metadata$Treatment
Sex = metadata$Sex
ExpRound = metadata$Exp.Round

##### Full 1. Visualize differences using PCA #####

# Step 1: Center and scale genotype matrix
datafullM = as.matrix(datafull)
datafullMS = scale(datafullM, center = TRUE, scale = TRUE)

# Step 2: Get data into pcadapt-readable format. Input = transposed genotype matrix
# (For data subset, will get warning message about having more individuals than SNPs)
pcadata = read.pcadapt(t(datafullMS)) 
screetest <- pcadapt(input=pcadata,K=20) # K = # of principal components to retain
plot(screetest,option="screeplot") + theme_bw()# Still appears that majority of variation NOT explained by first 10 pc axes
x = pcadapt(input= pcadata, K=5)

# Step 3: Plot PCA 
# Control vs heat-selected
plot(x,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
# Female vs Male
plot(x,option="scores", pop=Sex, 
     col = c("#e75480", "#89CFF0")) + theme_bw() + ggtitle(label = NULL)
# By Experimental Round (1-3)
plot(x,option="scores", pop=ExpRound, 
     col = c("#A50026", "#F99858", "#364B9A")) + theme_bw() + ggtitle(label = NULL)

##### Full. 2. Detect outlier SNPs between heat and control group #####
##### Full. 2a. Method 1: pcadapt (not using) #### 
# scree plot above indicates that majority of variation not explained by
# first ~10-20 PC axes, suggests pcadapt method not appropriate for 
# detecting SNP diffs between treatment and control

##### Full 2b. Method 2: Use Fst to find outliers #####

# Note: The call below was run in stages to avoid vector memory exhausted errors. 
# (This does not impact Fst estimates)
# fst = MakeDiploidFSTMat(datafull, locusNames = 1:ncol(datafull), popNames = metadata$Treatment)
# resulting Fst calculations were combined and output as dataframe
# which can be loaded here to skip this step in the future


fst1 = MakeDiploidFSTMat(datafull[,1:100000], locusNames = 1:100000, popNames = metadata$Treatment)
fst2 = MakeDiploidFSTMat(datafull[,100001:200000], locusNames = 100001:200000, popNames = metadata$Treatment)
fst3 = MakeDiploidFSTMat(datafull[,200001:300000], locusNames = 200001:300000, popNames = metadata$Treatment)
fst4 = MakeDiploidFSTMat(datafull[,300001:400000], locusNames = 300001:400000, popNames = metadata$Treatment)
fst5 = MakeDiploidFSTMat(datafull[,400001:500000], locusNames = 400001:500000, popNames = metadata$Treatment)
fst6 = MakeDiploidFSTMat(datafull[,500001:600000], locusNames = 500001:600000, popNames = metadata$Treatment)
fst7 = MakeDiploidFSTMat(datafull[,600001:700000], locusNames = 600001:700000, popNames = metadata$Treatment)
fst8 = MakeDiploidFSTMat(datafull[,700001:800000], locusNames = 700001:800000, popNames = metadata$Treatment)
fst9 = MakeDiploidFSTMat(datafull[,800001:900000], locusNames = 800001:900000, popNames = metadata$Treatment)
fst10 = MakeDiploidFSTMat(datafull[,900001:1000000], locusNames = 900001:1000000, popNames = metadata$Treatment)
fst11 = MakeDiploidFSTMat(datafull[,1000001:1100000], locusNames = 1000001:1100000, popNames = metadata$Treatment)
fst12 = MakeDiploidFSTMat(datafull[,1100001:ncol(datafull)], locusNames = 1100001:ncol(datafull), popNames = metadata$Treatment)


fstall = rbind(fst1, fst2, fst3, fst4, fst5, fst6, fst7, fst8, fst9, fst10, fst11, fst12)
write.csv(fstall, "~/Downloads/fstall.csv")








hist(fst$FST, breaks = 50, xlab = "Fst", ylab = "", main = "Fst value distribution: Method 2")

out1 <- OutFLANK(fst,NumberOfSamples = 2) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Not uniformly distributed

# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=out1$FSTNoCorrbar,
                                dfInferred=out1$dfInferred,qthreshold=0.05,Hmin=0.1)
# Which have q-values < 0.05
outliers <- which(P1$OutlierFlag==TRUE)
colnames(datasub)[outliers] # identified one SNP: V200615

# Manhattan plot with outlier SNPs
plot(P1$LocusName,P1$FST,xlab="SNP Position",ylab="FST",col=rgb(0,0,0,alpha=0.3), 
     pch = 16)
points(P1$LocusName[outliers],P1$FST[outliers],col="red", pch = 16)



