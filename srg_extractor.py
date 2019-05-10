#!/usr/bin/env python

"""
To build a synthetic reference genome in silico

Usage:
	./srg_extractor.py enzymeStart enzymeEnd minBpFragments maxbp genome species
	
	If using only one restriction enzyme, put it twice in the command line
	
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
import argparse,math,datetime,time,copy,sys,natsort
#parser = argparse.ArgumentParser()

def analyse(enzymeStart,enzymeEnd,minbp,maxbp,genome,species):

	dot = genome.index(".")
	ext = genome[dot+1:]
	outfile = genome[:dot]

# Creation du fichier de sortie:
#	file = open(outfile+"_bedfile_"+enzymeStart+"_"+enzymeEnd+"_"+minbp+"_"+maxbp+".bed", 'w')

	file = open(outfile+"_fragments.bed", 'w')
	
	try:
		test = RestrictionBatch([enzymeStart])
	except:
		print("The first enzyme is not a valid restriction enzyme")
		sys.exit(1)

	if enzymeEnd != enzymeStart:
		try:
			test.add(enzymeEnd)
		except:
			print("The second enzyme is not a valid restriction enzyme")
			sys.exit(1)
	else:
		pass

	rb = RestrictionBatch([enzymeStart])
	if enzymeEnd != enzymeStart:
		rb.add(enzymeEnd)
	else:
		pass

	for record in SeqIO.parse(genome, "fasta"):
#		print("\n\nSearching restriction sites in sequence " + record.name)
		sites = rb.search(record.seq)

		cuts=0;
		for key, value in sites.items():
			cuts += len(value)

		#In the case of a 2 enzymes analysis, we merge all the positions that those 2 enzymes cuts
		#We use Python sets which are more faster for searching a value in
		if(enzymeEnd != enzymeStart):
			allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = set(sites[rb.get(enzymeEnd)])
		else:
			allCuts = sites[rb.get(enzymeStart)]
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = {}

		#Put the positions in order
		allCuts.sort()

		i = 0
		for values in allCuts[:-1]:# We don't process the last cut
			#print(values)
			if len(cutsEnd) > 0:
				if(values in cutsStart and allCuts[i+1] in cutsEnd) or (values in cutsEnd and allCuts[i+1] in cutsStart):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					length=endPosition-beginingPosition
					if int(minbp)<=length<=int(maxbp):
						file.write(record.name + "\t" + str(beginingPosition) + "\t" + str(endPosition) + "\n")
					else:
						pass					
				else:
					pass
			else:
				beginingPosition = values - 1
				endPosition = allCuts[i+1] - 1
				length=endPosition-beginingPosition
				if int(minbp)<=length<=int(maxbp):
					file.write(record.name + "\t" + str(beginingPosition) + "\t" + str(endPosition) + "\n")
				else:
					pass
			i += 1
		
	file.close()

# Parsing user input
try:
    enzymeStart = sys.argv[1]
    enzymeEnd = sys.argv[2]
    minbp = sys.argv[3]
    maxbp = sys.argv[4]
    genome = sys.argv[5]
    species = sys.argv[6]
except:
	print(__doc__)
	sys.exit(1)

analyse(enzymeStart,enzymeEnd,minbp,maxbp,genome,species)
