#### Thermal Selection Experiment Sequence Analysis ####

setwd("~/Documents/Current Projects/Thermal Selection Experiment")

#### Create data subset (only need to run once)  #####
library(data.table)

setwd("~/Documents/Current Projects/Thermal Selection Experiment")
# data = fread("GenotypeMatrix/GenotypeMatrix", header = FALSE)
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

#### Load libraries and data files #####

library(pcadapt)
library(dplyr)
library(ggplot2)

# Load genotype matrix subset
datasub = read.csv("GenotypeMatrix/GenotypeMatrix_Subset.csv")[,-1]

# Load metadata
data = read.csv("~/Documents/Current Projects/Thermal Selection Experiment/FullAlignment&ExperimentStats.csv", header = T)
# Remove samples with missing KD or sequence data
Removes  = c("E-014", "E-05", "E-06", "E-07", "E-09")
data = data[!(data$Sample %in% Removes),] 
Treatment = data$Treatment
Sex = data$Sex
ExpRound = data$Exp.Round


#### 1. Visualize differences using PCA #####

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


#### 2. Detect outlier SNPs using pcadapt #### 

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


#### 3. Use Fst to find outliers #####

# Note: Fst were calculated in vcftools 
Fstdata = read.delim("Fst_estimates_VCFtools/Fst_estimates_controlvsheat_vcftools.txt")
colnames(Fstdata) = c("Chrom", "SNP", "WeirFst")
hist(Fstdata$WeirFst, xlab = "Weir Fst", ylab = "", 
     main = "Fst value distribution")

# Identify the mean Fst
summary(Fstdata$WeirFst)
meanFst = as.numeric(summary(Fstdata$WeirFst)[4])

# detect outliers as those exceeding 99th%
# following methods here: https://speciationgenomics.github.io/per_site_Fst/

my_threshold <- quantile(Fstdata$WeirFst, 0.999, na.rm = T)
fst <- Fstdata %>% mutate(outlier = ifelse(Fstdata$WeirFst > my_threshold, "outlier", "background"))
length(which(fst$outlier == "outlier")) # 1024 outliers detected

ggplot(fst, aes(SNP, WeirFst, colour = outlier)) + geom_point()


#### 4. Use lm to identify genotypic predictors of KD time #####

# Combine relevant metadata and genotype matrix into dataframe
df = cbind(metadata[,c(8,2,4)], datasub)
df$Treatment = as.factor(df$Treatment)
df$Sex = as.factor(df$Sex)

# Model KD time ~ genotype matrix + Treatment + Sex (not using body size for now) 
model = lm(Kdtime ~ ., data = df)

# Pull out p-values and adjust for multiple testing
pvals = as.numeric(summary(model)$coefficients[,4])

# padjusted = p.adjust(pvals, method = "fdr") # This appears to be over-correcting
outliers = which(pvals < 0.05) 

# Pull out ID of the SNP outliers
colnames(df)[outliers] # in this example: "V790007" "V254278" "V790890" "V248541" "V489649"




