### Thermal Selection Experiment Analysis 

library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)

setwd("~/Documents/Current Projects/Thermal Selection Experiment")

#### Pull in data #####

# All samples
data = read.csv("ExperimentData/ThermalKnockdownResults_AllSamples.csv", header = T)[-1]

# Limiting analysis to only individuals assayed 48-72h after eclosure 
# Removing assays that didn't include at least 1 control and 1 treatment individual
data = data[data$TimeSinceEclosed != "24h" & data$TimeSinceEclosed != "96h",]
data = data[data$Trial != "Y" & data$Trial != "X",]

# Sequenced samples only
dataseq = read.csv("ExperimentData/ThermalKnockdownResults_Sequenced&PassedFiltersSamplesOnly.csv", header = T)
dataseq = dataseq[!is.na(dataseq$KDtime),] # remove sample without KD time

#### Descriptive stats ####

# all samples
data %>% group_by(Sex, Treatment) %>% 
  dplyr::summarise(NumSamples = sum(Count)) %>% as.data.frame()
# 140 in control (66 F, 74 M), 111 in high-temp treatment (56 F, 55 M)

# sequenced samples only
dataseq %>% group_by(Sex, Treatment) %>% 
  dplyr::summarise(NumSamples = sum(Count)) %>% as.data.frame()
# 124 in control (61 F, 63 M), 107 in high-temp treatment (53 F, 54 M)

# Visualize variation in knockdown time
ggplot(dataseq, aes(x=Treatment, y=KDtime, fill=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(y = "Knockdown time (minutes)") + 
  scale_fill_manual(values=c("#5482AB", "firebrick3")) + 
  geom_jitter(color="black", size=1.0, alpha=0.9) +
  theme_bw() + theme(legend.position="none",
        plot.title = element_text(size=13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13)) 

# Visualize variation in body size
dataseqwl = dataseq[!is.na(dataseq$WingLength),]
ggplot(dataseqwl, aes(x=Treatment, y=WingLength, fill=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(y = "Wing length (mm)") + 
  scale_fill_manual(values=c("#5482AB", "firebrick3")) + 
  geom_jitter(color="black", size=1.0, alpha=0.9) +
  theme_bw() + theme(legend.position="none",
                     plot.title = element_text(size=13),
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 13)) 

# Range, mean, and sd knockdown times by group
data %>% group_by(Treatment) %>% 
  dplyr::summarise(RangeKD = range(KDtime, na.rm= T), 
         MeanKD = mean(KDtime, na.rm= T), sdKD = sd(KDtime, na.r = T)) %>% as.data.frame()

# Range, mean, and sd wing lengths by group
data %>% group_by(Treatment) %>% 
  dplyr::summarise(RangeWL = range(WingLength, na.rm= T), 
                   MeanWL = mean(WingLength, na.rm= T), sdKD = sd(WingLength, na.r = T)) %>% as.data.frame()


#### Difference in KD time #####

# Effect of treatment on knockdowm time 
# model includes experimental round and knockdown trial round as random effects
# treatment (variable of interest) and sex as fixed effects
mod = lmer(KDtime ~ Treatment + Sex + (1|ExpRound) + (1|Trial), data = dataseq)
summary(mod) # Significant effect of treatment (not of sex)
# On average, high temp treatment knocked down 2.85 minutes earlier than controls

# Effect of treatment on wing length
mod2 = lm(WingLength ~ Treatment + Sex, data = dataseq)
summary(mod2)
