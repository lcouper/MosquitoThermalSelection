#### Thermal Selection Experiment Downstream Analysis Pre Processing ####

#### Set wd and load libraries ####
library(pcadapt)
library(dplyr)
library(ggplot2)
library(OutFLANK)
library(tidyverse)
library(data.table)

setwd("~/Documents/Current Projects/Thermal Selection Experiment/Intermediate")

# Note the initial genotype matrix contained the following samples
# that were removed prior to downstream analysis due to missing data or sequencing error

Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09", "O-014")

#### Add "high" vs "low" kd time ranking to metadata file #####
# based on averages for that individuals treatment and sex
# averages calculated here: aggregate(metadata$Kdtime ~ metadata$Treatment + metadata$Sex, FUN = mean, na.rm = T)
HighLow = rep(NA, nrow(metadata))
for (i in 1:nrow(metadata))
{if (metadata$Sex[i] == "F" & metadata$Treatment[i] == "control")
    HighLow[i] = as.numeric(metadata$Kdtime[i] >= 48.49492)
else if (metadata$Sex[i] == "M" & metadata$Treatment[i] == "control")
  HighLow[i] = as.numeric(metadata$Kdtime[i] >= 50.78466)
else if (metadata$Sex[i] == "F" & metadata$Treatment[i] == "high temp")
  HighLow[i] = as.numeric(metadata$Kdtime[i] >= 45.10503)
else if (metadata$Sex[i] == "M" & metadata$Treatment[i] == "high temp")
  HighLow[i] = as.numeric(metadata$Kdtime[i] >= 47.47660)
}

HighLow = recode(HighLow, '0' = "Low", '1' = "High")
metadata$HighLow = HighLow

#### Combine intermediate files in which -1 has been replaced with NA ####
setwd("NAs")
NAfiles = dir(pattern = "*.csv")
datafullNAs <- NAfiles %>%
  map(fread, header = T) %>%    
  reduce(cbind)       

fwrite(datafullNAs, "GenotypeMatrix_WithNAs.csv") # output file

#### Identify columns to keep (those without NAs) ####

index = which(colSums(is.na(datafullNAs))==0)
index = as.numeric(index)
dataKeep = datafullNAs[, ..index]

# output file
fwrite(dataKeep , "GenotypeMatrix_NAsRemoved.csv")


#### Identify locus name to keep ####
loci = read.delim("~/Documents/Current Projects/Thermal Selection Experiment/GenotypeMatrix/PriorToRemovingNAs/GenotypeMatrix_LocusNames.txt", header = FALSE)
locusNames = as.character(paste(loci$V1, loci$V2, sep="_"))

locusNamesKeep = locusNames[index]

write.csv(locusNamesKeep, "~/Downloads/LocusNames_WithNAsRemoved.csv")


#### Combine Fst estimates from OutFlank (Treatment) ####

setwd("~/Documents/Current Projects/Thermal Selection Experiment/Intermediate/Fst")

Fstfiles = dir(pattern = "*.csv")
Fstall <- Fstfiles %>%
  map(fread, header = T) %>%    
  reduce(rbind)       

fwrite(Fstall, "fstall.csv")

#### Combine Fst estimates from OutFlank (KDtime) ####

setwd("~/Documents/Current Projects/Thermal Selection Experiment/Intermediate/Fst2")

f1 = fread("fst1mil.csv")
f2 = fread("fst2mil.csv")
f3 = fread("fst3mil.csv")
f4 = fread("fst4mil.csv")
f4 = f4[,-1]

fall = rbind(f1, f2, f3, f4)

fwrite(fall, "fstKDall.csv")

#### Create LD groups based on pruning done in plink #####

setwd("~/Documents/Current Projects/Thermal Selection Experiment/plink files")

LDs = fread("LD.plink.tags.list", header = T)[,-1]

# Keep only those in pruned SNP list
# First, bring in loci retained for analysis
loci = read.csv("../GenotypeMatrix/GenotypeMatrix_LocusNames_Pruned.csv", header = TRUE)[,-1]

LDs = LDs[LDs$SNP %in% loci,]
LDs$TAGS = gsub(":", "_", LDs$TAGS) # again, correct SNP name for consistency
LDs$TAGS = str_split(LDs$TAGS, "\\|") # creates list of linked SNPs for each focal SNP

# output as csv 
LDs %>%  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write.csv('LinkedSNPs_Directory.csv', row.names = FALSE)



#### Subset full data set based on LD pruning done in plink ######

setwd("~/Documents/Current Projects/Thermal Selection Experiment")

datafull = fread("GenotypeMatrix/GenotypeMatrix_NAsRemoved.csv", header = TRUE) # First row is sample number

# Load sample headers and locus names
header = read.csv("GenotypeMatrix/GenotypeMatrix_SampleNames.csv", header = FALSE)
SampleNames = as.character(header$V1)[-c(5, 155)] # remove E-014 and O-014 for consistent with above
rownames(datafull) = SampleNames
loci = read.csv("GenotypeMatrix/GenotypeMatrix_LocusNames_WithNAsRemoved.csv", header = TRUE)[,-1]
colnames(datafull) = loci

# Keep only SNPs that are retained after LD pruning in plink # 

# Load pruned dataset
pruned = read.table("plink files/plink.prune.in", quote="\"", comment.char="")
pruned = gsub(":", "_", pruned$V1) # correct name

# Remove all other SNPs from data full and loci dataframes above
loci = loci[which(loci %in% pruned)]
# Note datafull must be data frame for below code to work
datafull = as.data.frame(datafull)
datafull = datafull[,which(names(datafull) %in% pruned)]
dim(datafull) # 227 samples, 583,889 SNPs

write.csv(datafull, "GenotypeMatrix_Pruned.csv")
write.csv(loci, "GenotypeMatrix_LocusNames_Pruned.csv")

#### Generate a subset of the pruned dataset to upload to github ####

dataprunedsub = sample(datapruned, 100)
write.csv(dataprunedsub, "~/Downloads/LDpruned_GenotypeMatrix_Subset.csv")
