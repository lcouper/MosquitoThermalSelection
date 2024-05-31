################################################################################
##### Thermal Selection Experiment Sequence Analysis on LD-pruned data set #####  
##################### Written by Lisa Couper  ##################################
##################### Last updated: current  ###################################

# The code below identifies genes associated with focal SNPs #

##### 1. Load libraries and datasets #####

setwd("~/Documents/Current Projects/Thermal Selection Experiment")
library(grid)
library(dplyr)
library(ggVennDiagram)
library(gridExtra)
library(data.table)
library(ggplot2)
library(eulerr)

#### 2. [only need to do once] Identify number of genes. Remove redundant predictions ######

brakergenesAll <- read.delim("Gene Annotation/braker.gtf", header=FALSE)
colnames(brakergenesAll) <- c("Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "TranscriptID")

# How many genes?
brakergenesonly <- brakergenesAll %>% subset(Feature == "gene")

# remove redundant genes (identified by both Augustus and GeneMark)
transcriptRemove <- rep(NA, length(brakergenesonly))

for (i in 1:nrow(brakergenesonly)) {
  test <- brakergenesonly[i,]
  testChrom <- test$Chromosome
  allothers <- brakergenesonly[-i,]
  anydups <- length(which(allothers$Start <= test$Start & allothers$End >= test$Start))
  dupped <- allothers[which(allothers$Start <= test$Start & allothers$End >= test$Start),]
  dupped <- subset(dupped, Chromosome == testChrom)
if (nrow(dupped) == 0) 
  {transcriptRemove[i] <- NA}
 else if (anydups > 0 & dupped$Source == "GeneMark.hmm3") {
  transcriptRemove[i] <- dupped$TranscriptID}
else if (anydups > 0 & test$Source == "GeneMark.hmm3") { # prioritize keeping AUGUSTUS
  transcriptRemove[i] <- test$TranscriptID}
}

brakergenesonly <- brakergenesonly[-which(brakergenesonly$TranscriptID %in% transcriptRemove),]
#fwrite(brakergenesonly, "OutlierGenes/Braker_NonRedundantGenesOnly.csv")

# Identified 30,554 unique genes 

##### 2.  Pull in candidate SNP lists ######

brakergenes <- fread("OutlierGenes/Braker_NonRedundantGenesOnly.csv")

FstTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_Fst_Treat_pruned.csv")[,2]
plinkTreatSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_Treat_SexCovar_pruned_clumped.csv")[,2]
plinkKDSNPlist =  read.csv("OutlierSNPs/OutlierSNPnames_plink_KD_Covar_pruned_clumped.csv")[,2]

##### 3. Define function to return gene (and start & stop positions) associated with each SNP ####
GeneMatch <- function(x) {# here, x = a row in FstTreat
  # Subset to just the chromsome the SNP is on, and gene regions only
  chrom <- x$Chromosome
  brakergenesChromAug <- subset(brakergenes, Chromosome == chrom)  %>% subset(Feature == "gene") %>% subset(Source == "AUGUSTUS")
  brakergenesChromGM <- subset(brakergenes, Chromosome == chrom)  %>% subset(Feature == "gene") %>% subset(Source == "GeneMark.hmm3")
  
  # Identify gene the SNP belongs to based on start and end nucleotide positions
  pos <- as.numeric(x$SNPpos)
  GeneAugustus <- which(brakergenesChromAug$Start <= pos & brakergenesChromAug$End >= pos)
  GeneGM <- which(brakergenesChromGM$Start <= pos & brakergenesChromGM$End >= pos)
  
  # If SNP doesn't fall directly within a gene, assign it to the closest gene within 50 kb (if this exists)
  # Note this distance based on Selmoni et al. 2024 and Chappell et al. 2022 (used 5 kb but genomes 1/10th of size)
  # First, find closest gene
  GeneAugustus5kb <- brakergenesChromAug[which.min(abs(brakergenesChromAug$Start - pos ) + abs(brakergenesChromAug$End - pos )),]
  # Is it within 5k?
  distanceAug <- (abs(GeneAugustus5kb$Start - pos) <= 50000) | (abs(GeneAugustus5kb$End - pos) <= 50000)
  
  # Repeat for GeneMark
  GeneGM5kb <- brakergenesChromGM[which.min(abs(brakergenesChromGM$Start - pos ) + abs(brakergenesChromGM$End - pos )),]
  # Is it within 5k?
  distanceGM <- (abs(GeneGM5kb$Start - pos) <= 5000) | (abs(GeneGM5kb$End - pos) <= 5000)
  
  # Prefer the Augustus prediction if available, if not, use the GeneMark prediction
  
  if (length(GeneAugustus) > 0) 
  {return(c(brakergenesChromAug[GeneAugustus, "TranscriptID"], 
            brakergenesChromAug[GeneAugustus, "Start"],
            brakergenesChromAug[GeneAugustus, "End"]))}
  else if (length(GeneAugustus) == 0 & length(GeneGM) > 0)
  {return(c(brakergenesChromGM[GeneGM, "TranscriptID"], 
            brakergenesChromGM[GeneGM, "Start"],
            brakergenesChromGM[GeneGM, "End"]))}
  else if (length(GeneAugustus) == 0 & length(GeneGM) == 0 & distanceAug == TRUE)  
  {return(c(GeneAugustus5kb$TranscriptID, GeneAugustus5kb$Start, GeneAugustus5kb$End))}
  else if (length(GeneAugustus) == 0 & length(GeneGM) == 0 & distanceAug == FALSE & distanceGM == TRUE)
  {return(c(GeneGM5kb$TranscriptID, GeneGM5kb$Start, GeneGM5kb$End))}
  else 
  {return(c(NA, NA, NA))}
}

##### 4. For each SNP list, pull out chromosome and nucleotide position #####

##### 5. Fst Treat SNP list --> gene list #####

FstTreat = cbind.data.frame(substr(FstTreatSNPlist, 1, 8), substr(FstTreatSNPlist, 10, stop = 25)) # 25 is just to ensure the full name is captured. No SNPs that long probably
colnames(FstTreat) <- c("Chromosome", "SNPpos")

# Use GeneMatch function to find gene associated with each SNP
FstTreat$Gene <- NA
FstTreat$Start <- NA
FstTreat$End <- NA
for (i in 1:nrow(FstTreat)) {
  FstTreat$Gene[i] <- GeneMatch(FstTreat[i,])[1]
  FstTreat$Start[i] <- as.numeric(GeneMatch(FstTreat[i,])[2])
  FstTreat$End[i] <- as.numeric(GeneMatch(FstTreat[i,])[3])}

FstTreatGeneList <- FstTreat$Gene[!is.na(FstTreat$Gene)]
#fwrite(FstTreat, "OutlierGenes/FstTreatGeneList.csv")

# 311 / 351 SNPs map to a gene from braker gene annotation


##### 6. plink Treat SNP list --> gene list #####

plinkTreat = cbind.data.frame(substr(plinkTreatSNPlist, 1, 8), substr(plinkTreatSNPlist, 10, stop = 25)) # 25 is just to ensure the full name is captured. No SNPs that long probably
colnames(plinkTreat) <- c("Chromosome", "SNPpos")

# Use GeneMatch function to find gene associated with each SNP
plinkTreat$Gene <- NA
plinkTreat$Start <- NA
plinkTreat$End <- NA

for (i in 1:nrow(plinkTreat)) {
  plinkTreat$Gene[i] <- GeneMatch(plinkTreat[i,])[1]
  plinkTreat$Start[i] <- GeneMatch(plinkTreat[i,])[2]
  plinkTreat$End[i] <- GeneMatch(plinkTreat[i,])[3]}

plinkTreatGeneList <- plinkTreat$Gene[!is.na(plinkTreat$Gene)]

# 108 / 113 SNPs map to a gene from braker gene annotation
#fwrite(plinkTreat, "OutlierGenes/plinkTreatGeneList.csv")


##### 7. plink KDtime SNP list --> gene list #####

plinkKD = cbind.data.frame(substr(plinkKDSNPlist, 1, 8), substr(plinkKDSNPlist, 10, stop = 25)) # 25 is just to ensure the full name is captured. No SNPs that long probably
colnames(plinkKD) <- c("Chromosome", "SNPpos")

# Use GeneMatch function to find gene associated with each SNP
plinkKD$Gene <- NA
plinkKD$Start <- NA
plinkKD$End <- NA

for (i in 1:nrow(plinkKD)) {
  plinkKD$Gene[i] <- GeneMatch(plinkKD[i,])[1]
  plinkKD$Start[i] <- GeneMatch(plinkKD[i,])[2]
  plinkKD$End[i] <- GeneMatch(plinkKD[i,])[3]}

plinkKDGeneList <- plinkKD$Gene[!is.na(plinkKD$Gene)] 

# 114 / 120 SNPs map to a gene from braker gene annotation
#fwrite(plinkKD, "OutlierGenes/plinkKDGeneList.csv")

##### Calculate overlap in genes #####

#Overlaps <- fread("OutlierGenes/GeneList_AllApproaches.csv")

FstTreat <- fread("OutlierGenes/FstTreatGeneList.csv")
plinkTreat <- fread("OutlierGenes/plinkTreatGeneList.csv")
plinkKD <- fread("OutlierGenes/plinkKDGeneList.csv")

FstTreat <- FstTreat[!is.na(FstTreat$Start),]
plinkTreat <- plinkTreat[!is.na(plinkTreat$Start),]
plinkKD <- plinkKD[!is.na(plinkKD$Start),]

length(which(FstTreat$Gene %in% plinkTreat$Gene)) # 13 overlaps
length(which(FstTreat$Gene %in% plinkKD$Gene)) # 3 overlaps
length(which(plinkTreat$Gene %in% plinkKD$Gene)) # 3 overlaps

##### 8. Plot overlap in genes #####

SNPs = list(plinkTreatment = plinkTreat$Gene, FstTreat = FstTreat$Gene, plinkKDtime=plinkKD$Gene)
venn <- Venn(SNPs)
d <- process_data(venn)

colorGroups <- c(A = 'steelblue', B = 'navyblue', C = '#6A0DAD', D = "pink")     
colfunc <- colorRampPalette(colorGroups)
col <- colfunc(15)

ggplot() + 
  geom_sf(aes(fill = name), data = venn_region(d)) +
  geom_sf(aes(color = name), data = venn_setedge(d)) +
  geom_sf_text(aes(label = name), data = venn_setlabel(d)) +
  geom_sf_text(aes(label = count), data = venn_region(d)) +
  scale_fill_manual(values = alpha(col, .4)) +
  scale_color_manual(values = colorGroups) +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_x_continuous(expand = expansion(mult = .2))


fit <- euler(c("A" = 277, "B" = 89, "C" = 108, "A&B" = 13, "A&C" = 3, "B&C" = 3, 
               "A&B&C" = 0), shape = "ellipse")

plot(fit, labels = c("Fst-Treatment", "GWA-Treatment", "GWA-Knockdown"),
     fills = c(alpha("dodgerblue4", 0.7), alpha("lightblue", 0.7), alpha("#FF6865", 0.7)),
     fontsize = 12,  
     quantities = list(fontsize = 10))


                          
##### 9. Permutation: compare overlap to expectations under random chance ######


# for each focal gene with a given list (i.e., Fst Treat, plink Treat, plink KD time)
# subset the full gene list to only those:
# a) on the same chromosome and b) within one standard deviation of the length
# Then take one random draw from that subset
# Move on to next focal gene. Repeat for all within the given list.
# Repeat for other two lists.
# Identify how many overlapped
# Repeat for x # permutations 

brakergenes <- fread("OutlierGenes/Braker_NonRedundantGenesOnly.csv")
brakergenes$Length <- brakergenes$End - brakergenes$Start

l1data <- fread("OutlierGenes/FstTreatGeneList.csv")
l2data <- fread("OutlierGenes/plinkTreatGeneList.csv")
l3data <- fread("OutlierGenes/plinkKDGeneList.csv")

l1genes <- l1data[!is.na(l1data$Start),] # remove missing genes
l2genes <- l2data[!is.na(l2data$Start),]
l3genes <- l3data[!is.na(l3data$Start),]

l1genes$Length <- l1genes$End - l1genes$Start
l1genesSd <- sd(l1genes$Length)
l2genes$Length <- l2genes$End - l2genes$Start
l2genesSd <- sd(l2genes$Length)
l3genes$Length <- l3genes$End - l3genes$Start
l3genesSd <- sd(l3genes$Length)

##### Matched gene function #####
# For each focal gene, function pulls out one 'matching' genes 
# (i.e same chromosome and similar length) from the full gene set 

MatchedGenes <- function(x) { # x = a row from focal gene list
  Chrom <- x$Chromosome
  PossibleChrom <- brakergenes %>% subset(Chromosome == Chrom)
  LengthMin <- x$Length - l1genesSd
  LengthMax <- x$Length + l1genesSd
  PossibleGenes <- PossibleChrom %>% subset(Length >= LengthMin) %>% subset(Length <= LengthMax)
  Sample <- PossibleGenes[sample(1:nrow(PossibleGenes), size = 1),]
  return(Sample)
}

#### Create set of matching genes from each list (1x) #####

Match1 <- as.data.frame(matrix(nrow = nrow(l1genes), ncol = 10))
Match2 <- as.data.frame(matrix(nrow = nrow(l2genes), ncol = 10))
Match3 <- as.data.frame(matrix(nrow = nrow(l3genes), ncol = 10))

for (i in 1:nrow(l1genes)) {
  Match1[i,] <- MatchedGenes(l1genes[i,])
  colnames(Match1) <- colnames(brakergenes)
}
for (i in 1:nrow(l2genes)) {
  Match2[i,] <- MatchedGenes(l2genes[i,])
  colnames(Match2) <- colnames(brakergenes)
}
for (i in 1:nrow(l3genes)) {
  Match3[i,] <- MatchedGenes(l3genes[i,])
  colnames(Match3) <- colnames(brakergenes)
}

which(Match1$TranscriptID %in% Match2$TranscriptID) 
which(Match1$TranscriptID %in% Match3$TranscriptID) 
which(Match2$TranscriptID %in% Match3$TranscriptID) 

#### Permutation on above #####

# Create data frame to hold the number of overlapping genes between each set
overlaps <- as.data.frame(matrix(nrow = 500, ncol = 3))
colnames(overlaps) <- c("1&2", "1&3", "2&3")

for (j in 1:500) {
  
Match1 <- as.data.frame(matrix(nrow = nrow(l1genes), ncol = 10))
Match2 <- as.data.frame(matrix(nrow = nrow(l2genes), ncol = 10))
Match3 <- as.data.frame(matrix(nrow = nrow(l3genes), ncol = 10))

for (i in 1:nrow(l1genes)) {
  Match1[i,] <- MatchedGenes(l1genes[i,])
  colnames(Match1) <- colnames(brakergenes)}

for (i in 1:nrow(l2genes)) {
  Match2[i,] <- MatchedGenes(l2genes[i,])
  colnames(Match2) <- colnames(brakergenes)}

for (i in 1:nrow(l3genes)) {
  Match3[i,] <- MatchedGenes(l3genes[i,])
  colnames(Match3) <- colnames(brakergenes)}

overlaps[j,1] <- length(which(Match1$TranscriptID %in% Match2$TranscriptID))
overlaps[j,2] <- length(which(Match1$TranscriptID %in% Match3$TranscriptID))
overlaps[j,3] <- length(which(Match2$TranscriptID %in% Match3$TranscriptID))
}

#fwrite(overlaps, "OutlierGenes/FocalGenes_Permutation_ExpectedOverlap.csv")


##### Histograms of # overlapping genes from permutation #####

overlaps <- fread("OutlierGenes/FocalGenes_Permutation_ExpectedOverlap.csv")

# observed # of overlaps: 
# 1&2: 13; 1&3: 3, 2&3: 3

g1 <- ggplot(overlaps, aes(`1&2`)) + theme_minimal() + 
  geom_histogram(binwidth = 1, fill = "gray85", color = "black", position = "identity") +
  xlab("Fst-Treatment & GWA-Treatment") + 
  geom_vline(xintercept = 13, color = "red", linetype = "dashed", lwd = 1)

g2 <- ggplot(overlaps, aes(`1&3`)) + theme_minimal() + 
  geom_histogram(binwidth = 1, fill = "gray85", color = "black", position = "identity") +
  xlab("Fst-Treatment & GWA-Knockdown") + 
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", lwd = 1)

g3 <- ggplot(overlaps, aes(`2&3`)) + theme_minimal() + 
  geom_histogram(binwidth = 1, fill = "gray85", color = "black", position = "identity") +
  xlab("GWA-Treatment & GWA-Knockdown") + 
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", lwd = 1)

grid.arrange(g1,g2,g3, ncol=1, nrow = 3)

                                                                              