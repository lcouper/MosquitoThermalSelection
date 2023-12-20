###########################################################################
## Thermal Selection Experiment Sequence Analysis on LD-pruned data set  
## Written by Lisa Couper
## Last updated: 12/20/23 

# Code below does the following:
# Visualizes genomic differences using PCA
# Detects outlier SNPs using 3 appraoches (Fst, GWA)
# Visualizes candidate loci positions across the genome
# Examines allele frequency differnces in candidate loci relative to matched controls
# Examines signatures of selection (Tajima's D) across the genome


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
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
# Control vs heat-selected (axes 3 & 4)
plot(x,option="scores", pop=Treatment, i = 3, j = 4,
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
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

#fst = MakeDiploidFSTMat(datafull, locusNames = snpNames[!(snpNames %in% removes)], popNames = metadata$Treatment)

# resulting Fst calculations were combined and output as dataframe
# which can be loaded here to skip this step in the future

fstall = fread("OutlierSNPs/pruned_Fst_Treat_all.csv", header = T)[-1]

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
write.csv(P1, "OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv") 

# Obtain average Fst across all SNPs
mean(abs(P1$FST), na.rm = T) # 0.00378

# Identify SNPs with q-values <= 0.05 and Fst >= 0.05
outliers = P1[P1$qvalues <= 0.05 & P1$FST >= 0.05,] # should be 351

# Identify and output just the SNP names
outlierIndex = which(P1$qvalues <= 0.05 & P1$FST >= 0.05) 
OutlierSNPs = colnames(datapruned)[outlierIndex] 
#write.csv(OutlierSNPs, "OutlierSNPs/pruned_OutlierSNPnames_Fst_Treat.csv") # output this SNP list

##### Method 1 Manhattan Plot ######

# Pull in OutFlank results 
P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results

# Combine chromosome position and outlier status to create new column for coloring
P1$outlier = as.character(P1$FST >= 0.05)
P1$chrom = substr(P1$LocusName, 1, 8)    
P1$combo = paste(P1$chrom, P1$outlier, sep = "_")

# For manhattan plots to be comparable across methods, want to pull out bp position:
# for each SNP in P1, find its corresponding row number in loci
SNPindex = which(loci %in% P1$LocusName)

# Option 1: color based on outlier status
plot(x= SNPindex, y = P1$FST, xlab="SNP Position",ylab="FST", pch = 16, cex = 0.8,
     col=ifelse(P1$OutlierFlag==TRUE, "red", alpha("black", 0.1)))

# Option 2: color based on chromosome & outlier status
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex, y = P1$FST, xlab="SNP Position",ylab="FST", pch = 16, cex = 0.8,
     col= as.factor(P1$combo))
abline(h= 0.05)

# Keep only Fst > 0.025 for visual clarity, color by outlier status
P1FP = P1[P1$FST > 0.025,]
SNPindex2 = which(loci %in% P1FP$LocusName)

# Option 3: color based on outlier status only
plot(x = SNPindex2, y = P1FP$FST,xlab="SNP Position",ylab="FST", pch = 16, cex = 0.8,
     col=ifelse(P1FP$FST >= 0.05 & P1FP$qvalues <= 0.05, alpha("red", 0.6), alpha("black", 0.1)))

# Option 4: Color based on chromosome & outlier status
plot(x = SNPindex2, y = P1FP$FST,xlab="SNP Position",ylab="FST", pch = 16, cex = 0.8,
     col=as.factor(P1FP$combo))
abline(h= 0.05)


##### Method 2 Set-up: GWAS approach with treatment as predictor: Conducted in plink ####
# import full list of SNPs from plink. See here for details of plink file creation and analysis: https://github.com/lcouper/MosquitoThermalSelection/blob/main/Sequence_Analysis/DownstreamAnalysis_Plink/README.md

plinkTP <- read.csv("OutlierSNPs/plink_Treat_pruned", sep="")

# Remove any SNPs not in the original genotype matrix 
plinkTP$SNP = gsub(":", "_", plinkTP$SNP) # first, correct SNP name for consistency
plinkTP = plinkTP[-which(plinkTP$SNP %nin% colnames(datapruned)),]

# Pull in clumped dataset
plinkTPC = read.csv("OutlierSNPs/plink_Treat_pruned_clumped", sep="")
plinkTPC$SNP = gsub(":", "_", plinkTPC$SNP) # first, correct SNP name for consistency
# Correct GWA p-values for multiple testing
plinkTPC$P.adjust = p.adjust(plinkTPC$P, method = "BH")

# Merge dataframes (from clumped, take SNP, adjusted p-values, and linked SNP list)
plinkmerge = left_join(plinkTP, plinkTPC[,c(3,13, 12)], by = "SNP")

# Create column in plink TP to indicate whether each SNP found in clumped dataframe
# (clumped here contains all the SNPs identified as significant in the GWA)
# Later: add column to identify the SNPs in LD with the significant SNP
plinkmerge$sig = as.character(!is.na(plinkmerge$P.adjust))

# output merged plink file, as well list of SNPs identified as significant from clumping
plinkTPC_SNPlist = plinkTPC$SNP
#write.csv(plinkmerge, "OutlierSNPs/Treat_pruned_clumped_Plinkoutput.csv")
#write.csv(plinkTPC_SNPlist, "OutlierSNPs/OutlierSNPnames_plink_Treat_pruned_clumped.csv")

##### Method 2 Manhattan plot ######

plinkDF <- read.csv("OutlierSNPs/Treat_pruned_clumped_Plinkoutput.csv")[,-1]

# For consistency with Fst analysis, remove any SNPs not in Fst calculation dataframe
P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results
plinkDF = plinkDF[-which(plinkDF$SNP %nin% P1$LocusName),]

SNPindex2 = which(loci %in% plinkDF$SNP)

# Create column to indicate chromosome position and significance
plinkDF$combo = paste(plinkDF$CHR, plinkDF$sig, sep = "_")

# For plotting purposes only: create new column of p-values so all data points shown
plinkDF$pforplot = -log(plinkDF$P)

# Color based on chromosome & outlier status
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex2, y = plinkDF$pforplot, xlab="SNP Position",ylab="-log(association p value)", pch = 16, cex = 0.8,
     col= as.factor(plinkDF$combo))
abline(h = -log(max(plinkDF$P.adjust, na.rm = T)))

##### Method 3 Setup: GWAS approach with KDtime as predictor: Conducted in plink #####

plinkKDP = read.csv("OutlierSNPs/plink_KD_pruned", sep="")

# Remove any SNPs not in the original genotype matrix 
plinkKDP$SNP = gsub(":", "_", plinkKDP$SNP) # first, correct SNP name for consistency
plinkKDP = plinkKDP[-which(plinkKDP$SNP %nin% colnames(datapruned)),]

# Pull in clumped dataset (123 SNPs)
plinkKDPC = read.csv("OutlierSNPs/plink_KD_pruned_clumped", sep="")
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
#write.csv(plinkKDmerge, "OutlierSNPs/KD_pruned_clumped_Plinkoutput.csv")
#write.csv(plinkKDPC_SNPlist, "OutlierSNPs/OutlierSNPnames_plink_KD_pruned_clumped.csv")


##### Method 3 Manhattan Plot ######

plinkDF2 <- read.csv("OutlierSNPs/KD_pruned_clumped_Plinkoutput.csv")[,-1]

# For consistency with Fst analysis, remove any SNPs not in Fst calculation dataframe
P1 = read.csv("OutlierSNPs/OutFLANKresults_Treatment_Pruned.csv", header = T)[,-1]
P1 = P1[!is.na(P1$pvalues),] # remove rows with NAs in outflank results
plinkDF2 = plinkDF2[-which(plinkDF2$SNP %nin% P1$LocusName),]

SNPindex3 = which(loci %in% plinkDF2$SNP)

# Create column to indicate chromosome position and significance
plinkDF2$combo = paste(plinkDF2$CHR, plinkDF2$sig, sep = "_")

# For plotting purposes only: create new column of p-values so all data points shown
plinkDF2$pforplot = -log(plinkDF2$P)

# Color based on chromosome & outlier status
palette(c("#80CDC1","#018571", "#B2ABD2", "#5E3C99", "#92C5DE","#0571B0"))
plot(x= SNPindex3, y = plinkDF2$pforplot, xlab="SNP Position",ylab="-log(association p value)", pch = 16, cex = 0.8,
     col= as.factor(plinkDF2$combo))
abline(h = -log(max(plinkDF2$P.adjust, na.rm = T)))



##### Examine AF differences of SNPs identified above ####

# Using genotype matrix, calculate allele frequencies for each SNP and treatment
dfaf = cbind(metadata$Treatment, datapruned)
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
AFtable = cbind.data.frame(colnames(dfaf)[-1], af_control, af_heat)
colnames(AFtable) = c("SNP", "control", "heat")
AFtable$diff = AFtable$heat - AFtable$control
# Add column indicating if af increased (red) or decreased (blue) from contorl to heat-selected 
AFtable$incdec = as.factor(as.numeric(AFtable$diff >= 0))
# reshape wide to long for plotting
AFtable2 = AFtable %>% pivot_longer(-c(SNP, diff, incdec), names_to = "treatment", values_to = "af") 

###### Generated matched controls ######

# First, bring in data files containing list of outlier SNPs 
# including those from both OutFLANK Fst approach & plink Treatment
FstTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned.csv")[,2]
plinkTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_Treat_pruned_clumped.csv")[,2]

combinedSNPlist = c(FstTreatSNPlist, plinkTreatSNPlist)
# remove any non in AFtable
combinedSNPlist = combinedSNPlist[-which(combinedSNPlist %nin% AFtable$SNP)]

AFoutliers = AFtable2[which(AFtable2$SNP %in% combinedSNPlist), ] 

# Matched control requirements:
# a) same chromosome but >100kbp away from focal SNP 
# b) within 5% of starting af as focal SNP

MatchedControlFstDiff = matrix(nrow = length(combinedSNPlist), ncol = 11)

for (i in 1:length(combinedSNPlist)) {
  
  # pull out focal SNP
  focalSNP = as.character(combinedSNPlist[i]) 
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
AFmerge = left_join(AFoutliers, MatchedControlFstDiff, by = "SNP")
AFmerge$absdiff = abs(AFmerge$diff)

# Plot relative af change (showing af change for each focal SNP, relative to matched controls)
# First, keep only unique SNPs (as each SNP listed twice in AFoutliers dataframe)
Unique = AFmerge %>% distinct(SNP, .keep_all = TRUE)

# Find quantile of each focal SNP relative to its matched controls
percentile = NA
for (i in 1:nrow(Unique)){
  percentile[i] <- ecdf(as.numeric(Unique[i,6:15]))(Unique[i,16])}
Unique$percentile <- percentile

# output 
#write.csv(Unique, "~/Downloads/AF_Fst.csv")

##### K-S tests to determine if dists of AF differences differ #####

ks.test(Unique$absdiff, Unique$mc10) # compare to each of mc1 - mc10
ks.test(Unique$mc1, Unique$mc2) # try for a few pairs


###### Plot AF differences for Fst outlier SNPS ######

# Plot of AF differences between focal SNP and its matched controls
#UniqueLong = gather(Unique, AbsDiff, c(diffmc1, diffmc2, diffmc3, diffmc4, diffmc5, absdiff))

UniqueLong = gather(Unique, RealorMatched, AbsDiff, c(16,6:15))
mu <- plyr::ddply(UniqueLong, "RealorMatched", summarise, grp.mean=median(AbsDiff, na.rm = T)) # mean per group

ggplot(UniqueLong, aes(x=AbsDiff, fill = RealorMatched, color = RealorMatched)) + 
  geom_boxplot()+  xlab("Allele frequency shift") + theme_bw() + 
  coord_flip() + 
  scale_fill_manual(values=c("#008080", rep("grey90", 10))) + scale_color_manual(values=c("#008080", rep("grey90", 10)))

  



###### Examine Tajimas D ######

# Note Tajima's D calculation performed using vcftools
TajD <- read.delim("TajimasD_100kbControls")
plot(TajD$TajimaD, x = 1:nrow(TajD), cex = 0.1)

##### Plot Tajima's D & focal SNPs #####

# Color points above based on whether they were in focal SNP lists
# First, break TajD df apart into 3 chromosomes
TajDc1 = TajD[1:(grep("2_RagTag", TajD$CHROM)[1]-1),]
TajDc2 = TajD[grep("2_RagTag", TajD$CHROM)[1]:(grep("3_RagTag", TajD$CHROM)[1]-1),]
TajDc3 = TajD[grep("3_RagTag", TajD$CHROM)[1]:nrow(TajD),]

# Do same for focal SNP list (starting with Fst by treatment)
FstTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned.csv")[,2]

FstSNPsc1 = FstTreatSNPlist[grep("1_RagTag", FstTreatSNPlist)]
FstSNPsc1bp = as.numeric(sub(".*_", "", FstSNPsc1)) # keep just the bp portion
FstSNPsc2 = FstTreatSNPlist[grep("2_RagTag", FstTreatSNPlist)]
FstSNPsc2bp = as.numeric(sub(".*_", "", FstSNPsc2)) 
FstSNPsc3 = FstTreatSNPlist[grep("3_RagTag", FstTreatSNPlist)]
FstSNPsc3bp = as.numeric(sub(".*_", "", FstSNPsc3)) 

C1indexes = rep(NA, length(FstSNPsc1bp))
C2indexes = rep(NA, length(FstSNPsc2bp))
C3indexes = rep(NA, length(FstSNPsc3bp))
                
for (i in 1:length(FstSNPsc1bp)) {
  C1indexes[i] = which.min(abs(FstSNPsc1bp[i] - TajDc1$BIN_START))}
for (i in 1:length(FstSNPsc2bp)) {
  C2indexes[i] = which.min(abs(FstSNPsc2bp[i] - TajDc2$BIN_START))}
for (i in 1:length(FstSNPsc3bp)) {
  C3indexes[i] = which.min(abs(FstSNPsc3bp[i] - TajDc3$BIN_START))}

TajDc1$FocalFST = 0
TajDc1$FocalFST[C1indexes] <- 1
TajDc2$FocalFST = 0
TajDc2$FocalFST[C2indexes] <- 1
TajDc3$FocalFST = 0
TajDc3$FocalFST[C3indexes] <- 1

# Combine back together
TajDcombo = rbind(TajDc1, TajDc2, TajDc3)

# Plot and color based on Focal Fst flag and chromosome
TajDcombo$combo = paste(TajDcombo$CHROM, TajDcombo$FocalFST)
palette(c(alpha("#80CDC1", 0.3),"#018571", alpha("#B2ABD2",0.3), "#5E3C99", alpha("#92C5DE",0.3),"#0571B0"))
plot(TajDcombo$TajimaD, pch = 16, cex = 0.4, col = as.factor(TajDcombo$combo),
     ylab = "Tajima's D", xlab = "")

# Repeat for GWA - Treatment SNPs
plinkTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_Treat_pruned_clumped.csv")[,2]

plinkTreatSNPsc1 = plinkTreatSNPlist[grep("1_RagTag", plinkTreatSNPlist)]
plinkTreatSNPsc1bp = as.numeric(sub(".*_", "", plinkTreatSNPsc1)) 
plinkTreatSNPsc2 = plinkTreatSNPlist[grep("2_RagTag", plinkTreatSNPlist)]
plinkTreatSNPsc2bp = as.numeric(sub(".*_", "", plinkTreatSNPsc2)) 
plinkTreatSNPsc3 = plinkTreatSNPlist[grep("3_RagTag", plinkTreatSNPlist)]
plinkTreatSNPsc3bp = as.numeric(sub(".*_", "", plinkTreatSNPsc3)) 

C1aindexes = rep(NA, length(plinkTreatSNPsc1bp))
C2aindexes = rep(NA, length(plinkTreatSNPsc2bp))
C3aindexes = rep(NA, length(plinkTreatSNPsc3bp))

for (i in 1:length(plinkTreatSNPsc1bp)) {
  C1aindexes[i] = which.min(abs(plinkTreatSNPsc1bp[i] - TajDc1$BIN_START))}
for (i in 1:length(plinkTreatSNPsc2bp)) {
  C2aindexes[i] = which.min(abs(plinkTreatSNPsc2bp[i] - TajDc2$BIN_START))}
for (i in 1:length(plinkTreatSNPsc3bp)) {
  C3aindexes[i] = which.min(abs(plinkTreatSNPsc3bp[i] - TajDc3$BIN_START))}

TajDc1$FocalplinkTreat = 0
TajDc1$FocalplinkTreat[C1aindexes] <- 1
TajDc2$FocalplinkTreat = 0
TajDc2$FocalplinkTreat[C2aindexes] <- 1
TajDc3$FocalplinkTreat = 0
TajDc3$FocalplinkTreat[C3aindexes] <- 1

# Combine back together
TajDcombo2 = rbind(TajDc1, TajDc2, TajDc3)

# Plot and color based on Focal Fst flag and chromosome
TajDcombo2$combo = paste(TajDcombo2$CHROM, TajDcombo2$FocalplinkTreat)
palette(c(alpha("#80CDC1", 0.3),"#018571", alpha("#B2ABD2",0.3), "#5E3C99", alpha("#92C5DE",0.3),"#0571B0"))
plot(TajDcombo2$TajimaD, pch = 16, cex = 0.4, col = as.factor(TajDcombo2$combo),
     ylab = "Tajima's D", xlab = "")

# Repeat for GWA - KD SNPs
plinkKDSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_KD_pruned_clumped.csv")[,2]

plinkKDSNPsc1 = plinkKDSNPlist[grep("1_RagTag", plinkKDSNPlist)]
plinkKDSNPsc1bp = as.numeric(sub(".*_", "", plinkKDSNPsc1)) 
plinkKDSNPsc2 = plinkKDSNPlist[grep("2_RagTag", plinkKDSNPlist)]
plinkKDSNPsc2bp = as.numeric(sub(".*_", "", plinkKDSNPsc2)) 
plinkKDSNPsc3 = plinkKDSNPlist[grep("3_RagTag", plinkKDSNPlist)]
plinkKDSNPsc3bp = as.numeric(sub(".*_", "", plinkKDSNPsc3)) 

C1bindexes = rep(NA, length(plinkKDSNPsc1bp))
C2bindexes = rep(NA, length(plinkKDSNPsc2bp))
C3bindexes = rep(NA, length(plinkKDSNPsc3bp))

for (i in 1:length(plinkKDSNPsc1bp)) {
  C1bindexes[i] = which.min(abs(plinkKDSNPsc1bp[i] - TajDc1$BIN_START))}
for (i in 1:length(plinkKDSNPsc2bp)) {
  C2bindexes[i] = which.min(abs(plinkKDSNPsc2bp[i] - TajDc2$BIN_START))}
for (i in 1:length(plinkKDSNPsc3bp)) {
  C3bindexes[i] = which.min(abs(plinkKDSNPsc3bp[i] - TajDc3$BIN_START))}

TajDc1$FocalplinkKD = 0
TajDc1$FocalplinkKD[C1bindexes] <- 1
TajDc2$FocalplinkKD = 0
TajDc2$FocalplinkKD[C2bindexes] <- 1
TajDc3$FocalplinkKD = 0
TajDc3$FocalplinkKD[C3bindexes] <- 1

# Combine back together
TajDcombo3 = rbind(TajDc1, TajDc2, TajDc3)

# Plot and color based on Focal Fst flag and chromosome
TajDcombo3$combo = paste(TajDcombo3$CHROM, TajDcombo3$FocalplinkKD)
palette(c(alpha("#80CDC1", 0.3),"#018571", alpha("#B2ABD2",0.3), "#5E3C99", alpha("#92C5DE",0.3),"#0571B0"))
plot(TajDcombo3$TajimaD, pch = 16, cex = 0.4, col = as.factor(TajDcombo3$combo),
     ylab = "Tajima's D", xlab = "")

