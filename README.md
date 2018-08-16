# GBS Estimation and Skinny Reference Genome#

###**This work is based on an original idea of Davoud Torkamaneh, postdoc at Guelph University, Canada.**###  


As the name suggests, this project has two parts. The first part deals with the *in silico* estimation of what kind
of results one is expect to obtain when using a particular parameter set on a reference genome.
For each sequence (chromosome, contig, scaffold) of the reference genome, the script create a frequence distribution
length for the fragment obtained by *in sillico* digestion. 
 
The second part makes it possible to generate a new reference genome in which the regions 
that do not correspond to GBS fragments identified *in sillico* are masked. The goal is to reduce
execution time of a real GBS analysis when using this new reduced reference genome.  

## Estimation of GBS statistics from real genomes ##

The script created for this section is ```gbs_estimation.py```

The parameters one can submit to this analysis are:  

1. enzymeStart (First enzyme)  
2. enzymeEnd (Second enzyme (optional))  
3. meth (Methylation fraction)  
4. minbp (Minimal fragment size)  
5. maxbp (Maximal fragment size)  
6. cov (Expected coverage)  
7. mperc (Fraction of missing sites)  
8. infile (File containing genome to analyze)  
9. species (Species name)  

Parameters for the frequence distribution of fragment length:  

10. init_a_bin: initial value to start the distribution  
11. max_b_bin: maximal value of the distribution  
12. step_d_bin: step value  

Two output files are created  

1. a table containing the results from the in sillico digestion and the number of fragments obtained for each sequence of the reference genome  
2. a distribution length of the fragments obtained for each sequence of the reference genome  

## Building a reduced reference genome for GBS ##

Here is the procedure to follow to created a reduced representation of the genome

1. Running this script to create a bed file containing the regions of interest in the original reference genome:    
```./srg_extractor.py enzymeStart enzymeEnd minBpFragments maxBpFragments genome.fasta species```  
This will create a file called  
```genome_bedfile_enzymeStart_enzymeEnd_minBpFragments_maxBpFragments.bed ```  

2. Create a bed file for the original genome. This is simply a file containing the name and the length of each sequence:  
```./make_genome_file.py genome.fasta```  
This will create a file called ```genome.bed```  

3. Running bedtools complement to create a bed file uninterested regions in the original reference genome. This is the complement of the file create in step 1.  
```bedtools complement -i genome_bedfile_enzymeStart_enzymeEnd_minBpFragments_maxBpFragments.bed -g genome.bed > genome_bedfile_enzymeStart_enzymeEnd_minBpFragments_maxBpFragments_complement.bed```  

4. Running bedtools to mask uninterested regions:  
```bedtools maskfasta -fi genome.fasta -bed genome_bedfile_enzymeStart_enzymeEnd_minBpFragments_maxBpFragments_complement.bed -fo genome_skinny.fasta```  

5.  Indexing the new reference genome file  
```bwa index -a bwtsw genome_skinny.fasta```  
```samtools faidxgenome_skinny.fasta```  

6. Continue with fastgbs.  
