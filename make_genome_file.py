#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

"""
This script allows to create a bed file from the reference genome file.

It takes as input a fasta file containing all chromosome and scaffolds.

It gives as output a bed file containing the name of each sequence and the
total length:

chr1  1000
chr2  800
.     .
.     .
.     .
chr20 500

This file will be used with bedtools complement.

Usage:
           ./make_genome_file.py ref_genome.fasta

The name of the output file will be: ref_genome.bed

"""

import os, sys, fileinput
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC
from math import *
import numpy
class Record: pass
	
def length(genome):

	dot = genome.index(".")
	ext = genome[dot+1:]
	output = genome[:dot]

	t=open(output+".bed","w")

	for record in SeqIO.parse(genome, "fasta"):
		t.write(record.name+"\t"+str(len(record.seq))+"\n")
	t.close()

try:
	genome = sys.argv[1]
except:
	print("One or more option(s) are missing.")
	print(__doc__)
	sys.exit(1)
length(genome)