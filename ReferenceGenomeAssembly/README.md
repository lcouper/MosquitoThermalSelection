# Overview of Library Prep, Sequencing, and Computational Methods

# Library Prep

A single adult female Ae. sierrensis collected from Corvallis, Oregon was used for this assembly. 
DNA extraction, library prep and sequencing were conducted at the Oregon Genomics and Cell Characterization Core Facility. 

High molecular weight DNA extraction was performed using a Circulomics Nanobind Tissue Big DNA Kit (Cat No. 900-701-01). The resulting sample contained 2,088 ng of DNA with an average fragment size of 30.7kb.

Extracted DNA was used in library preparations - one following the SMRTbell Express Prep kit 2.0 With Low DNA Input protocol (Cat No. 101-730-400) and the other following the SMRTbell Express Prep kit 2.0 with Ultra-Low DNA Input protocol (Cat No. 101-987-800). For both preps, genomic DNA was first sheared to the target fragment size – 20kb for the low and 10kb for the ultra-low library – using the Megaruptor System. The low input library prep proceeded using 1,081ng of sheared genomic DNA, while the ultra-low prep proceeded using 20ng (the maximum permitted). Two parallel preps of the ultra-low library were used to increase the final yield for a total input of 40ng. 

Subsequent steps consisted of single-stranded overhangs, DNA damage repair, A-tailing, adaptor ligation. The low input library was then size-selected with BluePippin to omit <5kb fragments, followed by AMP purification. The ultra-low library was size-selected to omit <8kb fragments, followed by ProNex purification. Final library concentrations for the low input and ultra-low input libraries were 43.8 ng/µL and 86.8 ng/µL, respectively, and final library sizes were 14,378 bp and 11,133 bp. 

# Sequencing 

Libraries were each annealed to a sequencing primer V5 bound to a Sequel II DNA polymerase using a Sequel II Binding Kit 2.2. 
The low library was loaded onto a 8M SMRT Cell at on-plate concentration of 200 pM using diffusion loading. 
The ultra low library was loaded onto a 8M SMRT Cell at an on-plate concentration of  150 pM using adaptive loading. 
For both libraries, SMRT sequencing was performed on a Sequel II System with a 2 hour pre-extension time and 30 hour movie collection time.

# Computational 

A draft genome assembly was generated using Flye, using the low and ultra low input CCS (circular consensus sequence) (see /scripts/FlyeAssembly) 
Genome completeness was assessed using BUSCO, and haplotigs were identified and removed using purge_haplotigs (see /scripts/BUSCO and /scripts/purge_haplotigs)



