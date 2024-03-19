# Steps to identify genes from reference genome sequence 

### Step 1. Identify repeats
- Use RepeatModeler2 to construct a species-specific repeat library and mask the genome with RepeatMasker
- On SCG, version: repeatmodeler/2.0.1
- Note: RepeatModeler2 Can identify repeats de novo, But also requires an external database. Dfam is one that is free
- Potential help page: https://darencard.net/blog/2022-10-13-install-repeat-modeler-masker/
- Example of how to run repeatmodeler2: https://github.com/Dfam-consortium/RepeatModeler
- Some help from here: https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Basic
- Interpret results from repeat modeler: https://github.com/Dfam-consortium/RepeatModeler
- SCG script uesd: repeatmod.sbatch

```
module load repeatmodeler/2.0.1
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
BuildDatabase -name sierrensis asierrensis.scaffolded.fasta
RepeatModeler -database sierrensis -pa 16 &>run2.out
```


### Step 2. Mask repeats
- Use RepeatMasker to mask the repeats identified through RepeatModeler
- On SCG, version: repeatmasker/4.1.6
- Using example code from here: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.html#gsc.tab=0
- and here https://www.biostars.org/p/9495313/ 
- SCG script used: repeatmask.sbatch
  
```
module load repeatmasker/4.1.0
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
RepeatMasker -pa 16 -gff -lib sierrensis-families.fa asierrensis.scaffolded.fasta
```

- Note the above generates files including:
- **“asierrensis.scaffolded.fasta.masked”**, a fasta file like the original one but in which all bases of the detected repeats are replaced (masked) with the “N” letter.
- “asierrensis.scaffolded.fasta.tbl” that contains a brief of the detected repeats, providing details of the classified families and some statistics
- “asierrensis.scaffolded.fasta.out.gff” that provides the start and end position of each detected repeat
