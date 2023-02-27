# Mosquito Thermal Selection
 
These are the protocols data, and scripts used to investigate mosquito thermal adaptation as described in the manuscript: https://docs.google.com/document/d/19Ubu2EEQXl_nQds2mQTjXvEUH_yTWDDFenYHbrP_YZc/edit 

methods for detecting chromosomal inversions from sequence data if needed:
BreakDancerMax (BD) ver.1.4.4 with the default setting was used to validate breakpoints of SVs, including translocations, inversions and deletions at the nucleotide level using the WGS data (Binary Alignment/Map format). A Poisson model16 was used to calculate the confidence score for each candidate variant. BD is able to identify inter-chromosomal translocation (CTX), inversion (INV) and deletion (DEL)

method 2: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2252-9
