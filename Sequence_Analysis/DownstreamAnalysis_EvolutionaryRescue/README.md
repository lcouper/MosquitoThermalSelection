# Analysis of adaptive potential using evolutionary rescue model

Note: this analysis was conducted in R using the script 'ThermalSelectionExp_EvolRescue.R' uploaded here

In brief, this script:

- defines a range of plausible values for parameters: selection strength, heritability, phenotypic variance, and max population growth rate based on our experimental and genomic data and prior studies
- calculates the maximum evolutionary rate for each set of parameter values using the analytic, quantitative-genetic formulation from Lynch and Lande (1993):
- ![image](https://github.com/lcouper/MosquitoThermalSelection/assets/10873177/beec2a6f-6876-43af-9edf-cd52dce70790)

- plots evolutionary rates in comparison with rates of recently observed and projected warming (RCP 4.5 and 6.0)
  
