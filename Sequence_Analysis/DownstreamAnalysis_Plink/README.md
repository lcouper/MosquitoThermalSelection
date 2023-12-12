## Per
 
Create files for use in plink 
Followed guidance from: https://www.biostars.org/p/109690/ 
Another useful tutorial here: http://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/

As our genome involves many scaffolds, we first need to create our own chromosome mapping file. Do so, by running:
```
bcftools view -H Filtered_VCF_All_sorted_0.995_bialleliconly.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > filename.chrom-map.txt
```
Then, make mapping file by running:
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --plink --chrom-map filename.chrom-map.txt --out myplink
```

*Note: this did not actually work, so instead just ran: 
```
vcftools --vcf Filtered_VCF_All_sorted_0.995_bialleliconly.vcf --plink --out myplink
```

#### 25. Run plink
Followed tutorial here:
https://zzz.bwh.harvard.edu/plink/tutorial.shtml
Note: running plink v 1.9 on SCG 

```
# Navigate to working directory with .map and .ped files (here, /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF)
plink --file myplink
```

**Conduct association analysis with treatment as phenotype**
Here, treat.phe.txt (uploaded here) is a 3 column file specificying the family ID and individual ID (here the same thing) and the phenotype-- here '1' for control and '2' for heat-selected. Note that the --allow-no-sex flag is mandatory for this line to run
```
plink --file myplink --pheno --allow-no-sex treat.phe.txt --assoc
```

