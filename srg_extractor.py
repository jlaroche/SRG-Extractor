#!/prg/python/3.5/bin/python3.5

"""
To build a synthetic reference genome in silico

Usage:
	./srg_extractor.py enzymeStart enzymeEnd minBpFragments maxBpFragments genome species
	
	If using only one restriction enzyme, put it twice in the command line
	
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
import argparse,math,datetime,time,copy,sys,natsort
#parser = argparse.ArgumentParser()


def analyse(enzymeStart,enzymeEnd,minBpFragments,maxBpFragments,genome,species):

	dot = genome.index(".")
	ext = genome[dot+1:]
	outfile = genome[:dot]
	
# Creation du fichier de sortie:
	file = open(outfile+"_bedfile_"+enzymeStart+"_"+enzymeEnd+"_"+minBpFragments+"_"+maxBpFragments+".bed", 'w')

	rb = RestrictionBatch([enzymeStart])
	if(enzymeEnd != ''):
		rb.add(enzymeEnd)

#	n=1
	for record in SeqIO.parse(genome, "fasta"):
#		print("\n\nSearching restriction sites in sequence " + record.name)
		sites = rb.search(record.seq)
		for x in sites:
			print(record.name + "\t" + str(x) + "\t" +str(len(sites[x])) + " sites")
#			for y in sites[x]:
#				print("\t"+str(y))
#			print(sites)
		#summary.write(sites)

		cuts=0;
		list_value=[]
		for key, value in sites.items():
			cuts += len(value)
			list_value.append(str(sites[key]))
#			print(key)
#			summary.write(str(key)+ ": " + str(cuts) + " cutting sites found\n")
#			summary.write(sites[key])
	
# 2 RESTRICTION ENZYME USED
		if(len(rb) == 2):
			# Since we have 2 enzymes, we merge all the positions that those 2 enzymes cuts
			allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
			# Put the positions in order
			allCuts.sort()

			i = 0
			num = 1
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = set(sites[rb.get(enzymeEnd)])
			for values in allCuts:
				if(values != allCuts[-1]):     # We don't process the last cut
					if(values in cutsStart and allCuts[i+1] in cutsEnd):
						beginingPosition = values - 1
						endPosition = allCuts[i+1] - 1
						length=endPosition-beginingPosition
						if int(minBpFragments)<=length<=int(maxBpFragments):
#							file.write(">"+record.name + "_" + str(num) + "_" +str(beginingPosition) + "_" + str(endPosition) + " length:" + str(length) + "\n" + str(record.seq[beginingPosition:endPosition]) + "\n")
							file.write(record.name + "\t" + str(beginingPosition) + "\t" + str(endPosition) + "\n")
							num += 1
						else:
							pass
					else:
						pass
				else:
					pass
				i+=1

# 1 RESTRICTION ENZYME USED
		else:
			allCuts = sites[rb.get(enzymeStart)]
			allCuts.sort()

			i = 0
			num = 1
			for values in allCuts:
				if(values != allCuts[-1]):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					length=endPosition-beginingPosition
					if int(minBpFragments)<=length<=int(maxBpFragments):
#							file.write(">"+record.name + "_" + str(num) + "_" +str(beginingPosition) + "_" + str(endPosition) + " length:" + str(length) + "\n" + str(record.seq[beginingPosition:endPosition]) + "\n")
							file.write(record.name + "\t" + str(beginingPosition) + "\t" + str(endPosition) + "\n")
							num += 1
					else:
						pass
				else:
					pass
				i+=1
			print("\tValue num\t"  +str(num))
			print("\tValue i\t"  + str(i))

#		n+=1
		
	file.close()

# Parsing user input
try:
    enzymeStart = sys.argv[1]
    enzymeEnd = sys.argv[2]
    minBpFragments = sys.argv[3]
    maxBpFragments = sys.argv[4]
    genome = sys.argv[5]
    species = sys.argv[6]
except:
	print(__doc__)
	sys.exit(1)

analyse(enzymeStart,enzymeEnd,minBpFragments,maxBpFragments,genome,species)
