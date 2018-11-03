# SRG Extractor

A bioinformatic program designed to create a skinny reference genome (SRG) for reduced-representation sequencing (RRS) analysis.

This work is based on an original idea of Davoud Torkamaneh, postdoc at Guelph University, Canada.


# Introduction

Reduced-representation sequencing (RRS) is a genome-wide scanning method for simultaneous
discovery and genotyping of thousands to millions of SNPs that is used across a wide range
of species. However, in this method a reproducible but very small fraction of the genome is
captured for sequencing, however sequencing reads are typically aligned against the entire 
reference genome. Here we present a skinny reference genome (SRG) approach in which a 
simplified reference genome is used to decrease computing time for data processing and
to increase SNP counts and accuracy. A SRG can be integrated into any RRS analytical pipeline.  

# Requirements

1. Linux with parallel installed (http://www.gnu.org/software/parallel/)  
2. Python 2.7 or higher (https://www.python.org/) 
3. Biopython (https://biopython.org/)
4. bwa (https://github.com/lh3/bwa)  
5. samtools (http://www.htslib.org/)  
6. bedtools (https://bedtools.readthedocs.io/en/latest/)  
7. srg_extractor.py (this distribution) 


# SRG Extractor workflow

* **Stage 1:** in silico fragmentation
* **Stage 2:** Characterization of reference genome
* **Stage 3:** Identification of GBS-irrelevant regions
* **Stage 4:** Masking GBS-irrelevant regions and creating an SRG


Below is a schematic of the workflow, with inputs and outputs indicated for each stage.�





# Running SRG Extractor

## Stage 1: in silico fragmentation:

Create a bed file containing the genomic regions that can be produce through RRS library prep from the original reference genome:


```./srg_extractor.py enzymeStart enzymeEnd minBpFragments maxBpFragments genome.fasta species``` 
	
	
	
List of options:

		* enzymeStart: First enzyme     
		* enzymeEnd: Second enzyme 
		* minbp: Minimal fragment size  
		* maxbp: Maximal fragment size  
		* genome.fasta: Original reference genome sequence in FASTA format  
		
	
This will create a file called: 

```genome_fragments.bed ```  

	
Example for soybean:		
	
```./srg_extractor.py ApeKI ApeKI 50 1000 Gmax_275_v2_0.fasta soybean```


  
## Stage 2: Characterization of reference genome	

Create a bed file for the original reference genome. This is simply a file containing the name and the length of each sequence: 

```./make_genome_file.py genome.fasta``` 

This will create a file called:

```genome.bed```  
	
Example for soybean:
	
```./make_genome_file.py Gmax_275_v2_0.fasta```



## Stage 3: Identification of GBS-irrelevant regions

Create a bed file including GBS-irrelevant regions from the original reference genome:
	
```./gbs_irrelevent.py genome_fragments.bed genome.bed```

This will create a file called:
	
```genome_fragments_complement.bed```


 
## Stage 4: Masking GBS-irrelevant regions and creating an SRG

Mask the GBS-irrelevant regions and create an SRG:
	
```./masking.py genome.fasta genome_fragments_complement.bed```
	
This will create your SRG:
	
```SRG.fasta```


Example for soybean:
	
```./masking.py Gmax_275_v2_0.fasta genome_fragments_complement.bed```
  


## **Reminder:** same as any reference genome you should index the SRG before using:


Example:
	
```bwa index -a bwtsw SRG.fasta```  
```samtools faidx SRG.fasta```  



## Now you can use SRG as a reference genome in your RRS pipeline.  
