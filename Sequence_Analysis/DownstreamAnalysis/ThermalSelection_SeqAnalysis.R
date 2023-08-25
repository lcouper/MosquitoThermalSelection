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


##### Subset. 2c. Method 3: Model af differences (note: update this part with Mark) ####

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


##### Subset. 2d. Plot AF differences #####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfafs = cbind(metadata$Treatment, datasub)
colnames(dfafs)[1] = "Treatment"
dfafs$Treatment = as.factor(dfafs$Treatment)
rownames(dfafs) = metadata$Sample

# Separate control and heat-selected treatments
controls = dfafs[dfafs$Treatment == "control",-1]
heats = dfafs[dfafs$Treatment == "high temp",-1]

# calculate af
afs_control = colSums(controls)/ (2*nrow(controls))
afs_heat = colSums(heats)/ (2*nrow(heats))

# combine and calculate differences
afs_table = cbind.data.frame(afs_control, afs_heat)
colnames(afs_table) = c("control", "heat")
afs_table$diff = afs_table$heat - afs_table$control

# Identify outliers (i.e. SNPs with significant diff frequencies between control and heat)
thresholds = quantile(afs_table$diff, 0.99, na.rm = T) # Adjust this later
afs_outliers =  afs_table %>% mutate(outlier = ifelse(afs_table$diff > thresholds, "outlier", "background"))
rownames(afs_table)[which(afs_outliers$outlier == "outlier")] # 1 outlier detected: "V183021"

# Reshape and plot
afts = cbind(rownames(afs_table), afs_outliers[,-3])
colnames(afts)[1] = c("SNP")
afts2 = afts %>% pivot_longer(-c(SNP, outlier), names_to = "treatment", values_to = "af") 

ggplot(afts2,aes(x=treatment,y=af)) + theme_bw() +
   geom_point(aes(x = treatment, y = af, group = SNP, col = outlier)) + 
  scale_fill_manual(values = c("grey80", "red")) + 
   geom_line(aes(x = treatment, y = af, group = SNP, col = outlier)) +
  scale_color_manual(values = c("grey80", "red"))


##### Subset. 3. Use lm to identify genotypic predictors of KD time #####

# Combine relevant metadata and genotype matrix into dataframe
dfs = cbind(metadata[,c(8,2,4)], datasub)
dfs$Treatment = as.factor(dfs$Treatment)
dfs$Sex = as.factor(dfs$Sex)

# Model KD time ~ genotype matrix + Treatment + Sex (not using body size for now) 
models = lm(Kdtime ~ ., data = dfs)

# Pull out p-values and adjust for multiple testing
pvalsS = as.numeric(summary(models)$coefficients[,4])
padjS = p.adjust(pvalsS, method = "fdr") 
outliersS = which(padjS < 0.05) 

# Pull out ID of the SNP outliers
colnames(dfs)[outliersS] # in this example: none

# Second model excluding the Treatment and Sex predictors
dfs2 = cbind(metadata$Kdtime, datasub)
models2 = lm(dfs2$`metadata$Kdtime` ~ ., data = dfs2)
pvalsS2 = as.numeric(summary(models2)$coefficients[,4])
padjS2 = p.adjust(pvalsS2, method = "fdr") 
outliersS2 = which(padjS2 < 0.05) 

colnames(dfs2)[outliersS2] # still no significant SNPs identified in subset




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

##### Full. 2b. Method 2: Use Fst to find outliers #####

# Note: The call below was run in stages to avoid vector memory exhausted errors. 
# (This does not impact Fst estimates)
# fst = MakeDiploidFSTMat(datafull, locusNames = 1:ncol(datafull), popNames = metadata$Treatment)
# e.g of subset: fst1 = MakeDiploidFSTMat(datafull[,1:100000], locusNames = 1:100000, popNames = metadata$Treatment)
# resulting Fst calculations were combined and output as dataframe
# which can be loaded here to skip this step in the future

fstall = fread("fstall.csv", header = T)

hist(fstall$FST, breaks = 50, xlab = "Fst", ylab = "", main = "Fst value distribution")

out1 <- OutFLANK(fstall,NumberOfSamples = 2) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Not uniformly distributed

# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fstall,Fstbar=out1$FSTNoCorrbar,
                                dfInferred=out1$dfInferred,qthreshold=0.01,Hmin=0.1)
# Which have q-values < 0.05
outliers <- which(P1$OutlierFlag==TRUE)
length(outliers) # identified 56 SNP outliers based on Fst and q < 0.01
OutlierSNPs = colnames(datafull)[outlierIndex] 

# Manhattan plot with outlier SNPs
plot(P1$LocusName,P1$FST,xlab="SNP Position",ylab="FST",col=rgb(0,0,0,alpha=0.2), pch = 16)
points(P1$LocusName[outliers],P1$FST[outliers],col="red", pch = 16)




##### Full. 2d. Identify SNPs with largest AF differences ####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfaf = cbind(metadata$Treatment, datafull)
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
threshold = quantile(abs(af_table$diff), 0.9999, na.rm = T) # Adjust this later
af_outlier =  af_table %>% mutate(outlier = ifelse(abs(af_table$diff) > threshold, "outlier", "background"))
length(which(af_outlier$outlier == "outlier")) # identifies 119 SNPs in the top 99.99th percentile
outliers2 = rownames(af_table)[which(af_outlier$outlier == "outlier")]

# Reshape and plot
aft = cbind(colnames(dfaf)[-1], af_outlier)
colnames(aft)[1] = c("SNP")
aft2 = aft %>% pivot_longer(-c(SNP, outlier, diff), names_to = "treatment", values_to = "af") 
# Add column indicating if af increased (red) or decreased (blue) from contorl to heat-selected 
aft2$incdec = as.factor(as.numeric(aft2$diff >= 0))

# Keep only outliers
aft2o = aft2[aft2$outlier == "outlier",]

# Plot
ggplot(aft2o,aes(x=treatment,y=af)) + theme_bw() + 
  geom_point(aes(x = treatment, y = af, group = SNP, col = incdec)) + 
  scale_fill_manual(values = c(alpha("darkblue", 0.6), "darkred")) + 
  geom_line(aes(x = treatment, y = af, group = SNP, col = incdec)) +
  scale_color_manual(values = c(alpha("darkblue", 0.6), "darkred")) +
  ylab("allele frequency") + 
  theme(legend.position = "none", axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

##### Full. 2e. Plot AF differences of outlier SNPs identified through Fst ####

# Note: aft2 dataframe created above
# Pull out SNPs from this df identified as outliers based on Fsf
# Names of SNPs stored in OutlierSNPs (created in 2d)

afOut = aft2[aft2$SNP %in% OutlierSNPs , ]

ggplot(afOut, aes(x=treatment,y=af)) + theme_bw() + 
  geom_point(aes(x = treatment, y = af, group = SNP, col = incdec)) + 
  scale_fill_manual(values = c(alpha("darkblue", 0.6), "darkred")) + 
  geom_line(aes(x = treatment, y = af, group = SNP, col = incdec)) +
  scale_color_manual(values = c(alpha("darkblue", 0.6), "darkred")) +
  ylab("allele frequency") + 
  theme(legend.position = "none", axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

##### Full. 3. Use lm to identify genotypic predictors of KD time (troubleshoot with Mark) #####

# Combine relevant metadata and genotype matrix into dataframe
df = cbind(metadata[,c(8,2,4)], datafull)
colnames(df) = c("Kdtime", "Treatment", "Sex", colnames(datafull))
df$Treatment = as.factor(df$Treatment)
df$Sex = as.factor(df$Sex)

dftest = df[,-c(1000:1204990)]
# Model KD time ~ genotype matrix + Treatment + Sex (not using body size for now) 
model = lm(Kdtime ~ ., data = dftest)

# Pull out p-values and adjust for multiple testing
pvals = as.numeric(summary(model)$coefficients[,4])
padj = p.adjust(pvals, method = "bh") 
outliers = which(padj < 0.01) 

# Pull out ID of the SNP outliers
colnames(df)[outliers] # 

# Second model excluding the Treatment and Sex predictors
df2 = cbind(metadata$Kdtime, datafull)
model2 = lm(df2$`metadata$Kdtime` ~ ., data = df2)
pvals2 = as.numeric(summary(model2)$coefficients[,4])
padj2 = p.adjust(pvals2, method = "fdr") 
outliers2 = which(padj2 < 0.05) 

colnames(df2)[outliers2] 

#### Full. 4. Compare SNPs identified through Step 2b, 2d, and 3 ####

