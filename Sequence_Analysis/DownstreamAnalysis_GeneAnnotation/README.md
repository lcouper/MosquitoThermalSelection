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

### Step 3. Run BRAKER for gene prediction 
- Note: uses repeat masked fasta created in step 2
- Can be run with RNA-Seq and protein data if available (N/A here)
- On SCG, version: braker2/2.1.6
- Requirements:
- Augustus, On SCG, version: augustus/3.4.0. However, BRAKER requires the config file path to be writeable. To do this, had to install augustus into my SCG account (by installing on my local computer, then uploading)
- geneMark key in your home directory. Downloaded GeneMark-ES/ET/EP+ ver 4.72_lic  LINUX 64 kernel 2.6 - 4, 64 bit from here: http://exon.gatech.edu/GeneMark/license_download.cgi. Then uploaded to ref_genome folder as "gm_key"
- following guidance here: https://github.com/Gaius-Augustus/BRAKER?tab=readme-ov-file#running-braker and here: https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_to_Braker2.html#gsc.tab=0
- For this to run, had to previously change the path to perl environment in the GeneMark-ES pipeline:
- With: ref_genome/gmes_linux_64 as working directory:
```perl change_path_in_perl_scripts.pl "/usr/bin/env perl"```

```
module load braker2/2.1.6
cd /labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome
export AUGUSTUS_CONFIG_PATH=/labs/emordeca/ThermalSelectionExpSeqFiles/ref_genome/Augustus/config
braker.pl --genome=asierrensis.scaffolded.fasta.masked --esmode --threads 4
```
