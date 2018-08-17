# SRG Extractor#

### A bioinformatic program designed to create a skinny reference genome (SRG) for reduced-representation sequencing (RRS) analysis.###

###This work is based on an original idea of Davoud Torkamaneh, postdoc at Guelph University, Canada###


## Introduction ##

### Reduced-representation sequencing (RRS) is a genome-wide scanning method for simultaneous discovery and genotyping of thousands to millions of SNPs that is used across a wide range of species. However, in this method a reproducible but very small fraction of the genome is captured for sequencing, however sequencing reads are typically aligned against the entire reference genome. Here we present a skinny reference genome (SRG) approach in which a simplified reference genome is used to decrease computing time for data processing and to increase SNP counts and accuracy. A SRG can be integrated into any RRS analytical pipeline.###


## Dependencies##

### Linux with parallel installed (http://www.gnu.org/software/parallel/)###
### Python 2.7 or higher (https://www.python.org/)###
### bwa (https://github.com/lh3/bwa)###
### samtools (http://www.htslib.org/)###
### bedtools (https://bedtools.readthedocs.io/en/latest/) ###
### srg_extractor.py (this distribution) ###


## Running SRG Extractor ##

1. Create a bed file containing the regions of interest in the original reference genome:   

```./srg_extractor.py enzymeStart enzymeEnd minBpFragments maxBpFragments genome.fasta species```  

This will create a file called  

```genome_fragments.bed ```  

### enzymeStart: First enzyme ###   
### enzymeEnd: Second enzyme (optional))### 
### minbp: Minimal fragment size ### 
### maxbp: Maximal fragment size ### 
### genome.fasta: Original reference genome sequence in FASTA format


2. Create a bed file for the original genome. This is simply a file containing the name and the length of each sequence: 

```./make_genome_file.py genome.fasta``` 

This will create a file called 

```genome.bed```  

3. Create a bed file including uninterested regions from the original reference genome using bedtools complement. This is the complement of the file create in step 1.  

```bedtools complement -i genome_fragments.bed -g genome.bed > genome_fragments_complement.bed```  

4. Masking uninterested regions using bestools:  

```bedtools maskfasta -fi genome.fasta -bed genome_fragments_complement.bed -fo SRG.fasta```  

5.  Indexing the new reference genome file  

```bwa index -a bwtsw SRG.fasta```  
```samtools faidx SRG.fasta```  

6. Continue with fastgbs.  
