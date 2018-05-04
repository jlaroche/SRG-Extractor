# GBS Estimation and Synthetic Reference Genome#

As the name suggests, this project has two parts. The first part deals with the *in silico* estimation of what kind
of results one is expect to obtain when using a particular parameter set on a reference genomes.
The second part allows one to generate a synthetic reference genome from the preceding set
of parameters. The goal is to reduce execution time of a real GBS analysis when using this
new synthetic genome.

## Estimation of GBS statistics from real genomes ##

The parameters one can submit to this analysis are:
First enzyme
Second enzyme (optional)
Methylation fraction
Minimal fragment size
Maximal fragment size
Expected coverage
Fraction of missing sites

The script created for this section is ```estimate_gbs.py```

From each run, one obtain this 2 files:

One that the name begin with ```Resume_Estimate_GBS_``` and the other with ```Estimate_GBS_```


## Building synthetic reference genome for GBS ##

The script created for this section is ```build_synthetic_GBS_genome.py```


## Indexing synthetic reference genome ##

```bwa index -a bwtsw genome.fasta```

```for i in $(ls -1 *.fasta);do bwa index -a bwtsw $i;done```

```samtools faidx genome.fasta```

```for i in $(ls -1 *.fasta);do samtools faidx $i;done```
