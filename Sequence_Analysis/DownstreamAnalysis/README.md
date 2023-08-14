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
```

![Screeplot](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/ec6cc449-75e7-40a5-8d54-63c513fef785)


<p align="center">
  <img width="700"
    src="https://github.com/lcouper/MosquitoThermalSelection/tree/main/Sequence_Analysis/DownstreamAnalysis)/images/Screeplot.jpg">
  </p>   
  


```
# Plot PCA: Heat-selected vs control 
plot(x,option="scores", pop=Treatment, 
     col = c("#00AFBB",  "#FC4E07")) + theme_bw() + ggtitle(label = NULL)
```
