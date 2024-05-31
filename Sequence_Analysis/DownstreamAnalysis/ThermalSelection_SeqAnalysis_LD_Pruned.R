#############################################################################
#### Thermal Selection Experiment Sequence Analysis on LD-pruned data set ### 
#### Written by Lisa Couper #################################################

# Code below does the following:
# Loads libraries/files and creates data subset for sharing/reproducibility
# Visualizes genomic differences using PCA
# Detects outlier SNPs using 3 appraoches (Fst, GWA)
# Visualizes candidate SNPs positions across the genome
# Examines allele frequency diffs in candidate SNPs relative to matched controls
# Calculates and plots selection coefficients
# Re-shapes and plots population diversity metrics

##### Set wd and load libraries ####

library(pcadapt)
library(dplyr)
library(ggplot2)
library(OutFLANK)
library(tidyverse)
library(data.table)
library(Hmisc)
library(ggVennDiagram)
library(stringi)

setwd("~/Documents/Current Projects/Thermal Selection Experiment")


##### Load data files for analysis #####

# Note: for data subset, use 'GenotypeMatrix_Pruned_Subset.csv' and
# 'GenotypeMatrix_LocusNames_Pruned_Subset.csv' below
# all other files should be the same

datapruned = fread("GenotypeMatrix/GenotypeMatrix_Pruned.csv", header = TRUE)[,-1]

# Load sample headers and locus names
header = read.csv("GenotypeMatrix/GenotypeMatrix_SampleNames.csv", header = FALSE)
SampleNames = as.character(header$V1)[-c(5, 155)] # remove E-014 and O-014 for consistent with above
rownames(datapruned) = SampleNames
lociint = read.csv("GenotypeMatrix/GenotypeMatrix_LocusNames_Pruned.csv", header = TRUE)
loci = lociint$x
colnames(datapruned) = loci

# Load metadata and sample names
metadata = read.csv("FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09", "O-014")
metadata = metadata[!(metadata$Sample %in% Removes),] 
Treatment = metadata$Treatment
Sex = metadata$Sex
ExpRound = metadata$Exp.Round
HighLow = metadata$HighLow

##### Create data subset for sharing/reproducibility #####
index <- sample(1:ncol(datapruned), 100)
datasub <- as.data.frame(datapruned)[,sort(index)]
#write.csv(datasub, "GenotypeMatrix/GenotypeMatrix_Pruned_Subset.csv")
locisub <- loci[index]
#write.csv(locisub, "GenotypeMatrix/GenotypeMatrix_LocusNames_Pruned_Subset.csv")

##### Visualize differences using PCA #####

# Step 1: Center and scale genotype matrix
dataprunedM = as.matrix(datapruned)
dataprunedMS = scale(dataprunedM, center = TRUE, scale = TRUE)

# Step 2: Get data into pcadapt-readable format. Input = transposed genotype matrix
# (For data subset, will get warning message about having more individuals than SNPs)
pcadata = read.pcadapt(t(dataprunedMS)) 
screetest <- pcadapt(input=pcadata,K=20) # K = # of principal components to retain
plot(screetest,option="screeplot") + theme_bw()# Still appears that majority of variation NOT explained by first 10 pc axes
x = pcadapt(input= pcadata, K=5)

# Step 3: Plot PCA 
# Control vs heat-selected
plot(x,option="scores", pop=Treatment, 
     col = c( "#5482AB", "firebrick3")) + theme_bw() + ggtitle(label = NULL)
# Control vs heat-selected (axes 3 & 4)
plot(x,option="scores", pop=Treatment, i = 3, j = 4,
     col = c( "#5482AB", "firebrick3")) + theme_bw() + ggtitle(label = NULL)
# Female vs Male
plot(x,option="scores", pop=Sex, 
     col = c("#e75480", "#89CFF0")) + theme_bw() + ggtitle(label = NULL)
# By Experimental Round (1-3)
plot(x,option="scores", pop=ExpRound, 
     col = c("#A50026", "#F99858", "#364B9A")) + theme_bw() + ggtitle(label = NULL)
# By "High" and "Low" knockdown times
plot(x,option="scores", pop=HighLow, 
     col = c("#FC4E07", "#00AFBB")) + theme_bw() + ggtitle(label = NULL)

##### Detect outlier SNPs between heat and control group #####
##### Method 1 Set-up: Use Fst to find outliers based on treatment #####

# Note: The call below was run in stages to avoid vector memory exhausted errors. 
# (This does not impact Fst estimates)

fst1 = MakeDiploidFSTMat(datapruned[,1:100000], locusNames = loci[1:100000], popNames = metadata$Treatment)

# resulting Fst calculations were combined and output as dataframe
# which can be loaded here to skip this step in the future

fstall = fread("OutlierSNPs/pruned_Fst_Treat_all.csv", header = T)

hist(fstall$FST, breaks = 50, xlab = "Fst", ylab = "", main = "Fst value distribution")
out1 <- OutFLANK(fstall,NumberOfSamples = 2, RightTrimFraction = 0.1) # NumberOfSamples = number of populations
hist(out1$results$pvaluesRightTail) # Not uniformly distributed

# Plot observed Fst distribution with chi-squared fit
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

# Find statistical outliers
P1 <- pOutlierFinderChiSqNoCorr(fstall,Fstbar=out1$FSTNoCorrbar,
                                dfInferred=out1$dfInferred,qthreshold=0.05,Hmin=0.1)

# Output full OutFLANK results
fwrite(P1, "OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv") 

# Obtain average Fst across all SNPs
mean(abs(P1$FST), na.rm = T) # 0.00378

# Identify SNPs with q-values <= 0.05 and Fst >= 0.05
outliers = P1[P1$qvalues <= 0.05 & P1$FST >= 0.05,] # should be 351

# Identify and output just the SNP names
outlierIndex = which(P1$qvalues <= 0.05 & P1$FST >= 0.05) 
OutlierSNPs = colnames(datapruned)[outlierIndex] 
#write.csv(OutlierSNPs, "OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned") # output this SNP list

##### Method 1 Manhattan Plot ######

# Pull in OutFlank results 
P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
# note need to add back in chromosome/SNP tag
P1$LocusName <- read.csv("GenotypeMatrix/GenotypeMatrix_LocusNames_Pruned.csv", header = TRUE)[,-1]
P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results

# Combine chromosome position and outlier status to create new column for coloring
P1$outlier = as.character(P1$FST >= 0.05)
P1$chrom = substr(P1$LocusName, 1, 8)    
P1$combo = paste(P1$chrom, P1$outlier, sep = "_")

# For manhattan plots to be comparable across methods, want to pull out bp position:
# for each SNP in P1, find its corresponding row number in loci
SNPindex = which(loci %in% P1$LocusName)

# Keep only Fst > 0.025 for visual clarity, color by outlier status
P1FP = P1[P1$FST > 0.025,]
SNPindex2 = which(loci %in% P1FP$LocusName)

# Manhattan plot
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex2, y = P1FP$FST, xlab="SNP Position",ylab=expression("F"[st]), 
     pch = 16, cex = 0.8, las = 1,
     col= as.factor(P1FP$combo), frame.plot = FALSE, xaxt = "n")
abline(h = 0.05)


##### Method 2 Set-up: GWAS approach with treatment as predictor: Conducted in plink ####
# import full list of SNPs from plink. See here for details of plink file creation and analysis: https://github.com/lcouper/MosquitoThermalSelection/blob/main/Sequence_Analysis/DownstreamAnalysis_Plink/README.md

plinkTP <- read.csv("OutlierSNPs/plink_Treat_SexCovar_pruned", sep="")

# Remove any SNPs not in the original genotype matrix 
plinkTP$SNP = gsub(":", "_", plinkTP$SNP) # first, correct SNP name for consistency
plinkTP = plinkTP[-which(plinkTP$SNP %nin% colnames(datapruned)),]

# Pull in clumped dataset
plinkTPC = read.csv("OutlierSNPs/plink_Treat_SexCovar_pruned_clumped.txt", sep="")
plinkTPC$SNP = gsub(":", "_", plinkTPC$SNP) # first, correct SNP name for consistency
# Correct GWA p-values for multiple testing
plinkTPC$P.adjust = p.adjust(plinkTPC$P, method = "BH")

# Merge dataframes (from clumped, take SNP, adjusted p-values, and linked SNP list)
plinkmerge = left_join(plinkTP, plinkTPC[,c(3,13, 12)], by = "SNP")

# Create column in merged file to indicate whether each SNP found in clumped dataframe
# (clumped here contains all the SNPs identified as significant in the GWA)
# Later: add column to identify the SNPs in LD with the significant SNP
plinkmerge$sig = as.character(!is.na(plinkmerge$P.adjust))

# output merged plink file, as well list of SNPs identified as significant from clumping
plinkTPC_SNPlist = plinkTPC$SNP
#write.csv(plinkmerge, "OutlierSNPs/Treat_SexCovar_pruned_clumped_Plinkoutput.csv")
#write.csv(plinkTPC_SNPlist, "OutlierSNPs/OutlierSNPnames_plink_Treat_SexCovar_pruned_clumped.csv")

##### Method 2 Manhattan plot ######

plinkDF <- read.csv("OutlierSNPs/Treat_SexCovar_pruned_clumped_Plinkoutput.csv")[,-1]

# For consistency with Fst analysis, remove any SNPs not in Fst calculation dataframe
#P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
#P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results
plinkDF = plinkDF[-which(plinkDF$SNP %nin% P1$LocusName),]

SNPindex2 = which(loci %in% plinkDF$SNP)

# Create column to indicate chromosome position and significance
plinkDF$combo = paste(plinkDF$CHR, plinkDF$sig, sep = "_")

# For plotting purposes only: create new column of p-values so all data points shown
plinkDF$pforplot = -log(plinkDF$P)

# Color based on chromosome & outlier status
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex2, y = plinkDF$pforplot, xlab="SNP Position",ylab=expression("-log"[10]*italic("(p)")), 
     pch = 16, cex = 0.8, las = 1,
     col= as.factor(plinkDF$combo), frame.plot = FALSE, xaxt = "n")
abline(h = -log(max(plinkDF$P.adjust, na.rm = T)))

##### Method 3 Setup: GWAS approach with KDtime as predictor: Conducted in plink #####

plinkKDP = read.csv("OutlierSNPs/plink_KD_Covar_pruned", sep="")

# Remove any SNPs not in the original genotype matrix 
plinkKDP$SNP = gsub(":", "_", plinkKDP$SNP) # first, correct SNP name for consistency
plinkKDP = plinkKDP[-which(plinkKDP$SNP %nin% colnames(datapruned)),]

# Pull in clumped dataset (123 SNPs)
plinkKDPC = read.csv("OutlierSNPs/plink_KD_Covar_pruned_clumped", sep="")
plinkKDPC$SNP = gsub(":", "_", plinkKDPC$SNP) # first, correct SNP name for consistency
# Correct GWA p-values for multiple testing
plinkKDPC$P.adjust = p.adjust(plinkKDPC$P, method = "BH")


# Merge dataframes (from clumped, take SNP, adjusted p-values, and linked SNP list)

plinkKDmerge = left_join(plinkKDP, plinkKDPC[,c(3,13, 12)], by = "SNP")

# Create column in plink KDmerge to indicate whether each SNP found in clumped dataframe
# (clumped here contains all the SNPs identified as significant in the GWA)
# Later: add column to identify the SNPs in LD with the significant SNP
plinkKDmerge$sig = as.character(!is.na(plinkKDmerge$P.adjust))

# output merged plink file, as well list of SNPs identified as significant from clumping
plinkKDPC_SNPlist = plinkKDPC$SNP
#write.csv(plinkKDmerge, "OutlierSNPs/KD_Covar_pruned_clumped_Plinkoutput.csv")
#write.csv(plinkKDPC_SNPlist, "OutlierSNPs/OutlierSNPnames_plink_KD_Covar_pruned_clumped.csv")


##### Method 3 Manhattan Plot ######

plinkDF2 <- read.csv("OutlierSNPs/KD_Covar_pruned_clumped_Plinkoutput.csv")[,-1]

# For consistency with Fst analysis, remove any SNPs not in Fst calculation dataframe
#P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
#P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results
#plinkDF2 = plinkDF2[-which(plinkDF2$SNP %nin% P1$LocusName),]

SNPindex3 = which(loci %in% plinkDF2$SNP)

# Create column to indicate chromosome position and significance
plinkDF2$combo = paste(plinkDF2$CHR, plinkDF2$sig, sep = "_")

# For plotting purposes only: create new column of p-values so all data points shown
plinkDF2$pforplot = -log(plinkDF2$P)

# Color based on chromosome & outlier status
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex3, y = plinkDF2$pforplot, xlab="SNP Position",ylab=expression("-log"[10]*italic("(p)")), 
     pch = 16, cex = 0.8, las = 1,
     col= as.factor(plinkDF2$combo), frame.plot = FALSE, xaxt = "n")
abline(h = -log(max(plinkDF2$P.adjust, na.rm = T)))



##### Examine AF differences of SNPs identified above ####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfaf = cbind(metadata$Treatment, datapruned)
colnames(dfaf)[1] = "Treatment"
dfaf$Treatment = as.factor(dfaf$Treatment)
rownames(dfaf) = metadata$Sample

# Separate control and heat-selected treatments
control = dfaf[dfaf$Treatment == "control",-1]
heat = dfaf[dfaf$Treatment == "heat-selected",-1]

# calculate af
af_control = colSums(control)/ (2*nrow(control))
af_heat = colSums(heat)/ (2*nrow(heat))

# combine and calculate differences
AFtable = cbind.data.frame(colnames(dfaf)[-1], af_control, af_heat)
colnames(AFtable) = c("SNP", "control", "heat")
AFtable$diff = AFtable$heat - AFtable$control
# Add column indicating if af increased (red) or decreased (blue) from contorl to heat-selected 
AFtable$incdec = as.factor(as.numeric(AFtable$diff >= 0))
# reshape wide to long for plotting
AFtable2 = AFtable %>% pivot_longer(-c(SNP, diff, incdec), names_to = "treatment", values_to = "af") 

##### Generate matched controls for Method 1 outliers ######

# First, bring in data files containing list of outlier SNPs 
# including those from both OutFLANK Fst approach & plink Treatment
FstTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned.csv")[,2]

# remove any not in AFtable (does not remove any)
FstTreatSNPlist = FstTreatSNPlist[which(FstTreatSNPlist %in% AFtable$SNP)]

AFoutliers2 = AFtable2[which(AFtable2$SNP %in% FstTreatSNPlist), ] 

# Matched control requirements:
# a) same chromosome but >100kbp away from focal SNP 
# b) within 5% of starting af as focal SNP

MatchedControlFstDiff = matrix(nrow = length(FstTreatSNPlist), ncol = 11)

for (i in 1:length(FstTreatSNPlist)) {
  
  # pull out focal SNP
  focalSNP = as.character(FstTreatSNPlist[i]) 
  focalSNPindex = which(AFtable$SNP == focalSNP)
  
  # Pull out SNPs from same chromosome
  focalSNPchrom = substr(focalSNP, start = 1, stop = 8)
  focalchrom = AFtable[grepl(focalSNPchrom, AFtable$SNP),]
  
  # From those, pull out SNPs at least >100k bp away
  focalSNPpos = as.numeric(sub(".*_", "", focalSNP))
  focalchrom$pos = as.numeric(sub(".*_", "", focalchrom$SNP))
  focalchrom = focalchrom[(focalchrom$pos - focalSNPpos >= 100000),]
  
  # From those, pull out SNPs with af in control within 5% of focal SNP
  focalSNPaf = AFtable[focalSNPindex, "control"]
  focalchrom$prox = focalchrom$control >= (focalSNPaf * 0.975) & focalchrom$control <= (focalSNPaf * 1.025) 
  
  # Remove any not within 5%, as well as focal SNP itself
  potentialMatches = focalchrom[focalchrom$prox == "TRUE",]
  potentialMatches = potentialMatches[potentialMatches$SNP != focalSNP,]
  
  # if there are too few SNPs within this starting frequency - remove from list
  if (nrow(potentialMatches) <= 10){MatchedControlFstDiff[i] = NA}
  
  # Otherwise, take a random sample of 5 SNPs & pull af diffs for each of these matched controls
  else  {
    MatchedControls = potentialMatches[sample(1:nrow(potentialMatches), 10, replace = FALSE),]
    MatchedControlFstDiff[i,2:11] = abs(MatchedControls$diff)
    MatchedControlFstDiff[i,1] = focalSNP # add focal SNP name
  }
}

colnames(MatchedControlFstDiff) = c("SNP", "mc1", "mc2", "mc3", "mc4", "mc5",
                                    "mc6", "mc7", "mc8", "mc9", "mc10")
MatchedControlFstDiff = as.data.frame(MatchedControlFstDiff)
MatchedControlFstDiff[,2:11] <- sapply(MatchedControlFstDiff[,2:11],as.numeric)
AFmerge = left_join(AFoutliers2, MatchedControlFstDiff, by = "SNP")
AFmerge$absdiff = abs(AFmerge$diff)

# Plot relative af change (showing af change for each focal SNP, relative to matched controls)
# First, keep only unique SNPs (as each SNP listed twice in AFoutliers dataframe)
Unique2 = AFmerge %>% distinct(SNP, .keep_all = TRUE)

# Find quantile of each focal SNP relative to its matched controls
percentile = NA
for (i in 1:nrow(Unique2)){
  percentile[i] <- ecdf(as.numeric(Unique2[i,6:15]))(Unique2[i,16])}
Unique2$percentile <- percentile

# output 
write.csv(Unique2, "AlleleFrequencies/AF_Fstonly.csv")

##### Examine differences between focal SNPs and matched controls for Method 1 outliers #####

UniqueFst = read.csv("AlleleFrequencies/AF_Fstonly.csv", header = T)[,-1]

# k-s test to determine if distributions significantly different
ks.test(UniqueFst$absdiff, UniqueFst$mc10) # significantly different from most individually
ks.test(UniqueFst$absdiff, c(UniqueFst$mc1, UniqueFst$mc2, UniqueFst$mc3, UniqueFst$mc4, UniqueFst$mc5,
        UniqueFst$mc6, UniqueFst$mc7, UniqueFst$mc8, UniqueFst$mc9, UniqueFst$mc10)) # and from all MCs as a whole
ks.test(UniqueFst$mc1, UniqueFst$mc2) # try for a few pairs

# Plot
UniqueFstLong = gather(UniqueFst, RealorMatched, AbsDiff, c(16,6:15))
mu <- plyr::ddply(UniqueFstLong, "RealorMatched", summarise, grp.mean=median(AbsDiff, na.rm = T)) # median per group

# To consolidate matched controls:
UniqueFstLong$RealorMatched[grep("mc", UniqueFstLong$RealorMatched)] <- "mc"
mu <- plyr::ddply(UniqueFstLong, "RealorMatched", summarise, group.median=median(AbsDiff, na.rm = T)) # median per group
# median difference for focal SNPs: 0.02533177
# median difference for matched controls: 0.02318501

ggplot(UniqueFstLong, aes(x=AbsDiff, fill = RealorMatched)) + 
  geom_boxplot()+  xlab("Allele frequency shift") + theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=c("gray40", rep("grey80", 10))) + 
  scale_color_manual(values=c("gray40", rep("grey80", 10))) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=19,face="bold"))


##### Generate matched controls for Method 2 outliers ######

# First, bring in data files containing list of outlier SNPs 
plinkTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_Treat_SexCovar_pruned_clumped_Intermediate.csv")[,2]

# remove any not in AFtable (removed 9)
plinkTreatSNPlist = plinkTreatSNPlist[which(plinkTreatSNPlist %in% AFtable$SNP)]
#write.csv(plinkTreatSNPlist, "OutlierSNPs/OutlierSNPnames_plink_Treat_SexCovar_pruned_clumped.csv")

AFoutliers1 = AFtable2[which(AFtable2$SNP %in% plinkTreatSNPlist), ] 

# Matched control requirements:
# a) same chromosome but >100kbp away from focal SNP 
# b) within 5% of starting af as focal SNP

MatchedControlFstDiff = matrix(nrow = length(plinkTreatSNPlist), ncol = 11)

for (i in 1:length(plinkTreatSNPlist)) {
  
  # pull out focal SNP
  focalSNP = as.character(plinkTreatSNPlist[i]) 
  focalSNPindex = which(AFtable$SNP == focalSNP)
  
  # Pull out SNPs from same chromosome
  focalSNPchrom = substr(focalSNP, start = 1, stop = 8)
  focalchrom = AFtable[grepl(focalSNPchrom, AFtable$SNP),]
  
  # From those, pull out SNPs at least >100k bp away
  focalSNPpos = as.numeric(sub(".*_", "", focalSNP))
  focalchrom$pos = as.numeric(sub(".*_", "", focalchrom$SNP))
  focalchrom = focalchrom[(focalchrom$pos - focalSNPpos >= 100000),]
  
  # From those, pull out SNPs with af in control within 5% of focal SNP
  focalSNPaf = AFtable[focalSNPindex, "control"]
  focalchrom$prox = focalchrom$control >= (focalSNPaf * 0.975) & focalchrom$control <= (focalSNPaf * 1.025) 
  
  # Remove any not within 5%, as well as focal SNP itself
  potentialMatches = focalchrom[focalchrom$prox == "TRUE",]
  potentialMatches = potentialMatches[potentialMatches$SNP != focalSNP,]
  
  # if there are too few SNPs within this starting frequency - remove from list
  if (nrow(potentialMatches) <= 10){MatchedControlFstDiff[i] = NA}
  
  # Otherwise, take a random sample of 5 SNPs & pull af diffs for each of these matched controls
  else  {
    MatchedControls = potentialMatches[sample(1:nrow(potentialMatches), 10, replace = FALSE),]
    MatchedControlFstDiff[i,2:11] = abs(MatchedControls$diff)
    MatchedControlFstDiff[i,1] = focalSNP # add focal SNP name
  }
}

colnames(MatchedControlFstDiff) = c("SNP", "mc1", "mc2", "mc3", "mc4", "mc5",
                                    "mc6", "mc7", "mc8", "mc9", "mc10")
MatchedControlFstDiff = as.data.frame(MatchedControlFstDiff)
MatchedControlFstDiff[,2:11] <- sapply(MatchedControlFstDiff[,2:11],as.numeric)
AFmerge = left_join(AFoutliers1, MatchedControlFstDiff, by = "SNP")
AFmerge$absdiff = abs(AFmerge$diff)

# Plot relative af change (showing af change for each focal SNP, relative to matched controls)
# First, keep only unique SNPs (as each SNP listed twice in AFoutliers dataframe)
Unique1 = AFmerge %>% distinct(SNP, .keep_all = TRUE)

# Find quantile of each focal SNP relative to its matched controls
percentile = NA
for (i in 1:nrow(Unique1)){
  percentile[i] <- ecdf(as.numeric(Unique1[i,6:15]))(Unique1[i,16])}
Unique1$percentile <- percentile

# output 
# write.csv(Unique1, "AlleleFrequencies/AF_GWA_Treatonly.csv")


##### Examine differences between focal SNPs and matched controls for Method 2 outliers #####

UniqueGWAt = read.csv("AlleleFrequencies/AF_GWA_Treatonly.csv", header = T)[,-1]

# k-s test to determine if distributions significantly different
ks.test(UniqueGWAt$absdiff, UniqueGWAt$mc10) # significantly different
ks.test(UniqueGWAt$absdiff, c(UniqueGWAt$mc1, UniqueGWAt$mc2, UniqueGWAt$mc3, UniqueGWAt$mc4, UniqueGWAt$mc5,
                             UniqueGWAt$mc6, UniqueGWAt$mc7, UniqueGWAt$mc8, UniqueGWAt$mc9, UniqueGWAt$mc10)) # and from all MCs as a whole
ks.test(UniqueGWAt$mc1, UniqueGWAt$mc2) # try for a few pairs

# Plot
UniqueGWAtLong = gather(UniqueGWAt, RealorMatched, AbsDiff, c(16,6:15))
mu <- plyr::ddply(UniqueGWAtLong, "RealorMatched", summarise, grp.mean=median(AbsDiff, na.rm = T)) 

# To consolidate matched controls:
UniqueGWAtLong$RealorMatched[grep("mc", UniqueGWAtLong$RealorMatched)] <- "mc"
mu <- plyr::ddply(UniqueGWAtLong, "RealorMatched", summarise, grp.median=median(AbsDiff, na.rm = T)) # median per group
# median difference for focal SNPs: 0.15788447
# median difference for matched controls: 0.02589774

ggplot(UniqueGWAtLong, aes(x=AbsDiff, fill = RealorMatched)) + 
  geom_boxplot()  +  xlab("Allele frequency shift") + theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=c("gray40", rep("grey80", 10))) + 
  scale_color_manual(values=c("gray40", rep("grey80", 10))) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=19,face="bold"))

  




##### Generate matched controls for Method 3 outliers ######

# Here, we want to compare AF differences in the top and bottom 50% of the phenotypic distribution
# These "high" and "low" designations are based on the average for each treatment and sex
# and are identified in the 'ThermalSelectionExp_DataPreProcessing' script

# First, calculate allele frequencies for the high and low groups separately
dfaf2 = cbind(metadata$HighLow, datapruned)
colnames(dfaf2)[1] = "HighLow"
dfaf2$HighLow = as.factor(dfaf2$HighLow)
rownames(dfaf2) = metadata$Sample

low = dfaf2[dfaf2$HighLow == "Low",-1]
high = dfaf2[dfaf2$HighLow == "High",-1]

# calculate af
af_low = colSums(low)/ (2*nrow(low))
af_high = colSums(high)/ (2*nrow(high))

# combine and calculate differences
AFtable2 = cbind.data.frame(colnames(dfaf2)[-1], af_low, af_high)
colnames(AFtable2) = c("SNP", "low", "high")
AFtable2$diff = AFtable2$high - AFtable2$low
# Add column indicating if af increased (red) or decreased (blue) from low to high
AFtable2$incdec = as.factor(as.numeric(AFtable2$diff >= 0))
# reshape wide to long for plotting
AFtable3 = AFtable2 %>% pivot_longer(-c(SNP, diff, incdec), names_to = "treatment", values_to = "af") 


# Next, bring in data files containing list of outlier SNPs 
plinkKDSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_KD_Covar_pruned_clumped_Intermediate.csv")[,2]

# remove any not in AFtable (removes 3)
plinkKDSNPlist = plinkKDSNPlist[which(plinkKDSNPlist%in% AFtable2$SNP)]
#write.csv(plinkKDSNPlist, "OutlierSNPs/OutlierSNPnames_plink_KD_Covar_pruned_clumped.csv")

AFoutliers3 = AFtable3[which(AFtable3$SNP %in% plinkKDSNPlist), ] 

# Matched control requirements:
# same chromosome but >100kbp away from focal SNP 

MatchedControlFstDiff = matrix(nrow = length(plinkKDSNPlist), ncol = 11)

for (i in 1:length(plinkKDSNPlist)) {
  
  # pull out focal SNP
  focalSNP = as.character(plinkKDSNPlist[i]) 
  focalSNPindex = which(AFtable2$SNP == focalSNP)
  
  # Pull out SNPs from same chromosome
  focalSNPchrom = substr(focalSNP, start = 1, stop = 8)
  focalchrom = AFtable2[grepl(focalSNPchrom, AFtable2$SNP),]
  
  # From those, pull out SNPs at least >100k bp away
  focalSNPpos = as.numeric(sub(".*_", "", focalSNP))
  focalchrom$pos = as.numeric(sub(".*_", "", focalchrom$SNP))
  focalchrom = focalchrom[(focalchrom$pos - focalSNPpos >= 100000),]
  
  # From those, pull out SNPs with af in low within 5% of focal SNP
  potentialMatches = focalchrom[focalchrom$SNP != focalSNP,]
  
  # if there are too few SNPs within this starting frequency - remove from list
  if (nrow(potentialMatches) <= 10){MatchedControlFstDiff[i] = NA}
  
  # Otherwise, take a random sample of 10 SNPs & pull af diffs for each of these matched controls
  else  {
    MatchedControls = potentialMatches[sample(1:nrow(potentialMatches), 10, replace = FALSE),]
    MatchedControlFstDiff[i,2:11] = abs(MatchedControls$diff)
    MatchedControlFstDiff[i,1] = focalSNP # add focal SNP name
  }
}

colnames(MatchedControlFstDiff) = c("SNP", "mc1", "mc2", "mc3", "mc4", "mc5",
                                    "mc6", "mc7", "mc8", "mc9", "mc10")
MatchedControlFstDiff = as.data.frame(MatchedControlFstDiff)
MatchedControlFstDiff[,2:11] <- sapply(MatchedControlFstDiff[,2:11],as.numeric)
AFmerge = left_join(AFoutliers3, MatchedControlFstDiff, by = "SNP")
AFmerge$absdiff = abs(AFmerge$diff)

# Plot relative af change (showing af change for each focal SNP, relative to matched controls)
# First, keep only unique SNPs (as each SNP listed twice in AFoutliers dataframe)
Unique3 = AFmerge %>% distinct(SNP, .keep_all = TRUE)

# Find quantile of each focal SNP relative to its matched controls
percentile = NA
for (i in 1:nrow(Unique3)){
  percentile[i] <- ecdf(as.numeric(Unique3[i,6:15]))(Unique3[i,16])}
Unique3$percentile <- percentile

# output 
#write.csv(Unique3, "AlleleFrequencies/AF_GWA_KDonly.csv")

##### Examine differences between focal SNPs and matched controls for Method 3 outliers #####

UniqueGWAkd = read.csv("AlleleFrequencies/AF_GWA_KDonly.csv", header = T)[,-1]

# k-s test to determine if distributions significantly different
ks.test(UniqueGWAkd$absdiff, UniqueGWAkd$mc10) # significantly different
ks.test(UniqueGWAkd$absdiff, c(UniqueGWAkd$mc1, UniqueGWAkd$mc2, UniqueGWAkd$mc3, UniqueGWAkd$mc4, UniqueGWAkd$mc5,
                               UniqueGWAkd$mc6, UniqueGWAkd$mc7, UniqueGWAkd$mc8, UniqueGWAkd$mc9, UniqueGWAkd$mc10)) # and from all MCs as a whole
ks.test(UniqueGWAkd$mc1, UniqueGWAkd$mc2) # try for a few pairs

# Plot
UniqueGWAkdLong = gather(UniqueGWAkd, RealorMatched, AbsDiff, c(16,6:15))
mu <- plyr::ddply(UniqueGWAkdLong, "RealorMatched", summarise, grp.mean=median(AbsDiff, na.rm = T)) # mean per group

# To consolidate matched controls:
UniqueGWAkdLong$RealorMatched[grep("mc", UniqueGWAkdLong$RealorMatched)] <- "mc"

ggplot(UniqueGWAkdLong, aes(x=AbsDiff, fill = RealorMatched)) + 
  geom_boxplot()+  xlab("Allele frequency shift") + theme_bw() + 
  coord_flip() +
  scale_fill_manual(values=c("gray40", rep("grey80", 10))) + 
  scale_color_manual(values=c("gray40", rep("grey80", 10))) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))

##### Other summary statistics #####
##### Calculate and plot SNP-level selection coefficients ######

# note: this step uses AFtable2, created above

# Fst
FstTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned.csv")[,2]
AFoutliers2 = AFtable2[which(AFtable2$SNP %in% FstTreatSNPlist), ] 
SNP1controls <- AFoutliers2 %>% subset(treatment == "control")
SNP1controls$scoef <- (SNP1controls$diff / (SNP1controls$af * (1- SNP1controls$af)))

# Plot histogram (as proportions)
h = hist(abs(SNP1controls$scoef), plot = FALSE) 
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE, col = "grey40", las = 1, cex.lab = 1.5, cex.axis = 2.5,
     xlab = "", ylab = "", main = "")
abline(v = mean(abs(SNP1controls$scoef)), col = "firebrick", lwd = 3)
# mean s: , median s: 

# GWA-treatment
plinkTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_Treat_SexCovar_pruned_clumped_Intermediate.csv")[,2]
AFoutliers3 = AFtable2[which(AFtable2$SNP %in% plinkTreatSNPlist), ] 
SNP2controls <- AFoutliers3 %>% subset(treatment == "control")
SNP2controls$scoef <- (SNP2controls$diff / (SNP2controls$af * (1- SNP2controls$af)))

# Plot histogram (as proportions)
h = hist(abs(SNP2controls$scoef), plot = FALSE) 
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE, col = "grey40", las = 1, 
     cex.lab = 1.5, cex.axis = 2.5,
     xlab = "", ylab = "", main = "")
abline(v = mean(abs(SNP2controls$scoef)), col = "firebrick", lwd = 3)


# GWA-KD time
plinkKDSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_KD_Covar_pruned_clumped.csv")[,2]
AFoutliers4 = AFtable2[which(AFtable2$SNP %in% plinkKDSNPlist), ] 
SNP3controls <- AFoutliers4 %>% subset(treatment == "control")
SNP3controls$scoef <- (SNP3controls$diff / (SNP3controls$af * (1- SNP3controls$af)))

# Plot histogram (as proportions)
h = hist(abs(SNP3controls$scoef), plot = FALSE) 
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE, col = "grey40", las = 1, 
     cex.lab = 1.5, cex.axis = 2.5,
     xlab = "", ylab = "", main = "")
abline(v = mean(abs(SNP3controls$scoef)), col = "firebrick", lwd =3)


##### Population diversity metrics #####

# 1. Inbreeding coefficient (F)

hetero <- read.delim("~/Documents/Current Projects/Thermal Selection Experiment/PopDiv/output_het.het")
# remove individuals not used

Removes = c("E-014.deduped.bam", "O-014.deduped.bam")
hetero <- hetero[!(hetero$INDV %in% Removes),]
mean(hetero$F) # -0.06167967
# Negative values = less inbreeding than expected by chance

# 2. Observed and expected heterozygosity 
hetero$O.HET <- (hetero$N_SITES - hetero$O.HOM)/hetero$N_SITES
hetero$E.HET <- (hetero$N_SITES - hetero$E.HOM)/hetero$N_SITES

mean(hetero$O.HET) # 0.3934614
mean(hetero$E.HET) # 0.3706028

hist(hetero$O.HET, xlab = "Observed heterozygosity", main = "")
abline(v = mean(hetero$E.HET), col = "red", lwd = 1.5)

# 3. Nucleotide diversity (pi)

pi <- read.delim("~/Documents/Current Projects/Thermal Selection Experiment/PopDiv/all_samples_pi.windowed.pi")

mean(pi$PI) # 0.001520873
palette(c("#018571",  "#5E3C99", "#0571B0"))
plot(pi$PI, col = as.factor(pi$CHROM), 
     ylab = "\u03C0")
abline(h = mean(pi$PI), lty = 2)

