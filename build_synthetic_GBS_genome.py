#!/prg/python/3.5/bin/python3.5

import argparse
#parser = argparse.ArgumentParser()

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
import math
import datetime
import time
import copy

def analyse(enzymeStart,enzymeEnd,minBpFragments,maxBpFragments,infile,species):

	dot = infile.index(".")
	ext = infile[dot+1:]
	outfile = infile[:dot]
	
# Creation de 2 fichiers:
	file = open(outfile+"_synthetic_genome_with_"+enzymeStart+"_"+enzymeEnd+".fasta", 'w')
	summary = open(outfile+"_synthetic_genome_with_"+enzymeStart+"_"+enzymeEnd+".summary", 'w')

#	summary.write("Loading restriction enzymes...")
	rb = RestrictionBatch([enzymeStart])
	if(enzymeEnd != ''):
		rb.add(enzymeEnd)

	n=1
	for record in SeqIO.parse(infile, "fasta"):
		summary.write("\n\nSearching restriction sites in sequence " + record.name)
		sites = rb.search(record.seq)
		for x in sites:
			summary.write("\n" + str(x) + " (" +str(len(sites[x])) + " sites) ")
			for y in sites[x]:
				summary.write("\t"+str(y))
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
	
	# If 2 restriction enzymes selected
		if(len(rb) == 2):
			# Since we have 2 enzymes, we merge all the positions that those 2 enzymes cuts
			allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
			# Put the positions in order
			allCuts.sort()

#			print(allCuts)

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
						if minBpFragments<length<maxBpFragments:
							summary.write("\nKeeping a fragment between cutting positions: "+str(values)+" - "+str(allCuts[i+1]))
							file.write(">synthfrag" + str(num) + " "  + record.name + " " + str(beginingPosition) + "_" + str(endPosition) + " length:" + str(length) + "\n" + str(record.seq[beginingPosition:endPosition]) + "\n")
							num += 1
						else:
							pass
					else:
						pass
				else:
					pass
				i+=1
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
					if minBpFragments<length<maxBpFragments:
						file.write(">synthfrag" + str(num) + " " + record.name + " " + str(beginingPosition) + "_" + str(endPosition) + " length:" + str(length) + "\n" + str(record.seq[beginingPosition:endPosition]) + "\n")
						num += 1
					else:
						pass
				else:
					pass
				i+=1
		n+=1
		
	file.close()
	summary.close()

if __name__ == '__main__':
	analyse('PstI','MspI',80,325,'budworm.fasta','Spruce budworm')
#	analyse('ApeKI','ApeKI',80,325,'budworm.fasta','Spruce budworm')
