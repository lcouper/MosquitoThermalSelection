## Old results from permutation method on glm
![Rplot06](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/79fc35eb-0c04-4693-9dfc-e9431adf0c07)

**This method identified 73 outlier SNPs. These appear to cluster around similar areas of the genome as with the Fst approach, but No SNPs were identified as significant/outliers via both methods**


## 3. Examine allele frequency differences between control and heat-selected group
First, calculated frequency of alternative allele for each SNP separately for control and heat-selected.  
Next, calculated difference in allele frequency between these groups.   

### Visualize difference in allele frequency for SNPs identified as outliers through Fst or glm approaches above
![FstGlmOutlierSNPs](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/47318857-8fd3-4ef7-a9a5-7d14d16749de)

### 3b. For each outlier SNP, compare AF differences relative to matched controls
To generate matched controls: 
1. Identify SNPs within 1000 bp of focal SNP
2. Keep those within 5% of starting (control) frequency as focal SNP
3. From this list, take random sample of 7 SNPs
4. Caculate af differences (i.e. between control and heat-treated group) for these 7 SNPs
5. Compare af differences between focal SNP and matched controls

#### Visualize difference in AF for focal SNPs relative to matched controls
Here, pointed indicate focal SNPs with larger (red) or smaller (blue) AF differences than their matched controls

![AFrelativeDiff](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/b24caa4f-820b-438c-bc3f-8d859dbcd8a3)
