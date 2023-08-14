## Downstream sequence analysis
*Pertains to analysis performed in R, after genotype matrix was created*

1. Visualize overall differences between treatment groups using PCA 

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

  

### Plot PCA: Heat-selected vs control 
```
plot(x,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot139](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/8dc59047-a108-4873-bb05-2e33ae073723)


### Plot PCA: Female vs Male
```
plot(x,option="scores", pop=Sex, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot140](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/f84fdfe1-a313-4c5d-933f-b92edc9b49ff)

### Plot PCA: Experimental Round (1-3)
```
plot(x,option="scores", pop=ExpRound, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```

![Rplot141](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/877eb0e6-991d-4afa-a734-83a2625859f7)


