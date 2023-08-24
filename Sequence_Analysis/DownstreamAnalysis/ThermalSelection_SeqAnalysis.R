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


#### 2. Detect outlier SNPs between heat and control group #####
#### 2a. Method 1: pcadapt #### 

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


#### 2b. Method 2: Use Fst to find outliers #####

### see code on new lab laptop for how to do this!!

#### 2c. Method 3: Use AF differences #####

dfaf = cbind(metadata$Treatment, datasub)
colnames(dfaf)[1] = "Treatment"
dfaf$Treatment = as.factor(dfaf$Treatment)
rownames(dfaf) = metadata$Sample

control = dfaf[dfaf$Treatment == "control",-1]
heat = dfaf[dfaf$Treatment == "high temp",-1]

af_control = colSums(control)/ (2*nrow(control))
af_heat = colSums(heat)/ (2*nrow(heat))

af_table = cbind.data.frame(af_control, af_heat)
colnames(af_table) = c("control", "heat")
af_table$diff = af_table$heat - af_table$control

# Identify outliers 
threshold = quantile(af_table$diff, 0.99, na.rm = T) # Adjust this later
af_outliers =  af_table %>% mutate(outlier = ifelse(af_table$diff > threshold, "outlier", "background"))
which(af_outliers$outlier == "outlier") # 1 outlier detected: "V183021"

# Reshape and plot
aft = cbind(rownames(af_table), af_table[,-3])
colnames(aft)[1] = c("SNP")
aft2 = aft %>% pivot_longer(-SNP, names_to = "treatment", values_to = "af") 

ggplot(aft2,aes(x=treatment,y=af)) +
    geom_point() + 
   geom_line(aes(x = treatment, y = af, group = SNP))





#### 4. Use lm to identify genotypic predictors of KD time #####

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

