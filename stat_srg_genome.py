#!/usr/bin/env python3.5
# -*- coding: iso-8859-1 -*-

"""

This script compute basic statistics in a SRG genome

Usage:
          ./Stat_srg_genome.py srg_genome.fasta

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

	
def count_nuc(genome):
	dot = genome.index(".")
	ext = genome[dot+1:]
	output = genome[:dot]

	out=open(output+"_SRG_stat.txt","w")
	out.write("ID\tNb of Nuc.\tRRS_irrelevant\tSRG\tSRG%\n")
	
	tot_L,tot_N,tot_S=0,0,0	
	for record in SeqIO.parse(genome, "fasta"):
		liste=[]
		liste.append(record.name)
		liste.append(str(len(record.seq)))
		nb_N = record.seq.count('N')+record.seq.count('n')
		nb_A = record.seq.count('A')+record.seq.count('a')
		nb_T = record.seq.count('T')+record.seq.count('t')
		nb_C = record.seq.count('C')+record.seq.count('c')
		nb_G = record.seq.count('G')+record.seq.count('g')
		
		tot_SRG=nb_A+nb_T+nb_C+nb_G
		
		tot_L+=len(record.seq)
		tot_N+=nb_N
		tot_S+=tot_SRG
		
		liste.append(str(nb_N))
		liste.append(str(nb_A+nb_T+nb_C+nb_G))
		liste.append(str(round(((nb_A+nb_T+nb_C+nb_G)/len(record.seq))*100,2))+'%')
		out.write('\t'.join(liste)+"\n")

	out.write('Total'+'\t'+str(tot_L)+'\t'+str(tot_N)+'\t'+str(tot_S)+'\t'+str(round((tot_S/tot_L)*100,2))+'%'+'\n')
	out.close()
	
try:
	genome = sys.argv[1]

except:
	print("One or more option(s) are missing.")
	print(__doc__)
	sys.exit(1)
count_nuc(genome)
