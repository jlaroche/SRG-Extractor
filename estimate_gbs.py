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

#Enzymes choices: PstI	ApeKI

def analyse(enzymeStart,enzymeEnd,methylationPercent,minBpFragments,maxBpFragments,coverage,missingPercent,infile,species):

	dot = infile.index(".")
	ext = infile[dot+1:]
	outfile = infile[:dot]
	
# Creation de 2 fichiers: 
	stat = open("Estimate_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minBpFragments)+"_"+str(maxBpFragments)+".out", 'w')
	resume = open("Resume_Estimate_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minBpFragments)+"_"+str(maxBpFragments)+".out", 'w')

	resume.write("Species\tFile\tenzymeStart-enzymeEnd\tNb cutting sites\tNb fragments\t%Methyl\tNb fragments\tSite mutation\tSize limits\tNb fragments\tCoverage (X)\tMissing fragments\n")

	stat.write("Analysis" + "\n")
	stat.write("Name of file containing genome: " + infile + "\n")
	stat.write("Name of species: " + species + "\n")
	stat.write("Current date and time: " + str(datetime.datetime.now()) + "\n")
	stat.write("Start enzyme: " + enzymeStart + "\n")
	stat.write("End enzyme: " + enzymeEnd + "\n")
	stat.write("Percentage of methylation: "+ str(methylationPercent) + "\n")
	stat.write("Minimum base pair fragments: " + str(minBpFragments) + "\n")
	stat.write("Maximum base pair fragments: " + str(maxBpFragments) + "\n")
	stat.write("Average coverage: " + str(coverage) + "\n")
	stat.write("Percentage of missing: " + str(missingPercent) + "\n")

#	Loading restriction enzymes
	rb = RestrictionBatch([enzymeStart])
	if(enzymeEnd != ''):
		rb.add(enzymeEnd)
#	print(rb)

#	Loading sequences in memory
	genome = ''
	for record in SeqIO.parse(infile, "fasta"):
		genome = genome + str(record.seq)
	
#	Create a sequence object
	my_seq = Seq(genome)
	stat.write("\n")
	stat.write("Total bases in file (Genome's Length) : " + str(len(my_seq)) + " bp")
	stat.write("\n")

#	Searching restriction sites
	sites = rb.search(my_seq)
#	print(sites)

	cuts=0;
	for key, value in sites.items():
		cuts += len(value)

	stat.write("\nNumber of cutting sites found : " + str(cuts) + "\n")

	# If 2 restriction enzymes selected
	if(len(rb) == 2):
		# Since we have 2 enzymes, we merge all the positions that those 2 enzymes cuts
		allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
		# Put the positions in order
		allCuts.sort()

		i = 0
		frag = 0
		num = 1
		listlength = []

		# We use Python sets which are more faster for searching a value in
		cutsStart = set(sites[rb.get(enzymeStart)])
		cutsEnd = set(sites[rb.get(enzymeEnd)])
		
		stat.write("\nFrom the list of all cutting sites, scan for site pair where the first site belong to startEnzyme and the second site belong to endEnzyme\n")
		
		for values in allCuts:
			if(values != allCuts[-1]):     # We don't process the last cut
				if(values in cutsStart and allCuts[i+1] in cutsEnd):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					listlength.append(len(my_seq[beginingPosition:endPosition]))
					frag += 1
			num += 1
			i += 1
		stat.write("Number of fragments found that begins with " + enzymeStart + " and end with " + enzymeEnd + " : " + str(frag)+"\n")

	else:
		allCuts = sites[rb.get(enzymeStart)]
		allCuts.sort()

		i = 0
		num = 1
		listlength = []
		stat.write("\nFrom the list of all cutting sites, scan for site pair where the first site belong to startEnzyme and the second site belong to endEnzyme\n")
		for values in allCuts:
			if(values != allCuts[-1]):
				beginingPosition = values - 1
				endPosition = allCuts[i+1] - 1
				listlength.append(len(my_seq[beginingPosition:endPosition]))
			i += 1
			num += 1

		stat.write("Number of fragments found that begins with " + enzymeStart + " and end with " + enzymeEnd + " : " + str(cuts + 1)+"\n")
		frag=copy.deepcopy(cuts)

	cutsMeth = frag - ((frag * methylationPercent) / 100)
	stat.write("\nRemoving " + str(round((frag * methylationPercent)/100)) + " fragments according to methylation. Remaining : " + str(int(round(cutsMeth))))
	stat.write("\n")

	formula = -0.81 + 0.68 * (math.log(len(my_seq), 10));
	cutsFormula = cutsMeth - formula
	stat.write("\nRemoving " + str(round(formula)) + " fragments according to site mutation. Remaining : " + str(int(round(cutsFormula))))
	stat.write("\n")

	toRemove = 0
	for values in listlength:
		if(values <= minBpFragments or values >= maxBpFragments):
			toRemove += 1
	smallCuts = cutsFormula - toRemove;
	stat.write("\nRemoving " + str(toRemove) + " fragments less than " + str(minBpFragments) + " bp" + " and greater than " + str(maxBpFragments) + " bp" + ". Remaining : " + str(int(round(smallCuts))))
	stat.write("\n")

	cutsCoverage = smallCuts * coverage;
	stat.write("\nMultiply by depth of coverage (" + str(coverage) + "x). Remaining : " + str(int(round(cutsCoverage))))
	stat.write("\n")

	cutsMissing = cutsCoverage - ((cutsCoverage * missingPercent) / 100);
	stat.write("\nRemoving " + str(round(((cutsCoverage * missingPercent) / 100))) + " fragments due to missing (" + str(missingPercent) + "%). Remaining : " + str(int(round(cutsMissing))))
	stat.write("\n")


	resume.write(species + "\t" + infile + "\t" + enzymeStart + "-" + enzymeEnd + "\t" + str(cuts) + "\t" + str(frag) + "\t" + str(methylationPercent) + "\t" + str(round((frag * methylationPercent)/100)) + "\t" + str(round(formula)) + "\t" + str(minBpFragments)+"_" + str(maxBpFragments) + "\t"+ str(toRemove) + "\t" + str(coverage) + "\t" + str(round(((cutsCoverage * missingPercent) / 100))) + "\n")

	stat.close()
	resume.close()

if __name__ == '__main__':
#	analyse(enzymeStart,enzymeEnd,methylationPercent,minBpFragments,maxBpFragments,coverage,missingPercent,infile)

#	analyse('PstI','MspI',15,60,260,5,10,'barley_IBSC_PGSB_v2.fasta','Barley')
#	analyse('ApeKI','ApeKI',15,60,260,5,10,'barley_IBSC_PGSB_v2.fasta','Barley')
#	analyse('RsaI','RsaI',15,60,260,5,10,'barley_IBSC_PGSB_v2.fasta','Barley')
#	analyse('SbfI','SbfI',15,60,260,5,10,'barley_IBSC_PGSB_v2.fasta','Barley')

#	analyse('PstI','MspI',15,60,260,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('ApeKI','ApeKI',15,60,260,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('RsaI','RsaI',15,60,260,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('SbfI','SbfI',15,60,260,5,10,'Gmax_275_v2_0.fasta','Soybean')

#	analyse('PstI','MspI',15,60,300,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('ApeKI','ApeKI',15,60,300,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('RsaI','RsaI',15,60,300,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('SbfI','SbfI',15,60,300,5,10,'Gmax_275_v2_0.fasta','Soybean')

#	analyse('PstI','MspI',15,60,360,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('ApeKI','ApeKI',15,60,360,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('RsaI','RsaI',15,60,360,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('SbfI','SbfI',15,60,360,5,10,'Gmax_275_v2_0.fasta','Soybean')

#	analyse('PstI','MspI',15,60,400,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('ApeKI','ApeKI',15,60,400,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('RsaI','RsaI',15,60,400,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('SbfI','SbfI',15,60,400,5,10,'Gmax_275_v2_0.fasta','Soybean')

#	analyse('PstI','MspI',15,60,460,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('ApeKI','ApeKI',15,60,460,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('RsaI','RsaI',15,60,460,5,10,'Gmax_275_v2_0.fasta','Soybean')
#	analyse('SbfI','SbfI',15,60,460,5,10,'Gmax_275_v2_0.fasta','Soybean')

#	analyse('PstI','MspI',4,60,260,5,10,'bosTau8.fasta','Bovine')
#	analyse('ApeKI','ApeKI',4,60,260,5,10,'bosTau8.fasta','Bovine')
#	analyse('RsaI','RsaI',4,60,260,5,10,'bosTau8.fasta','Bovine')
#	analyse('SbfI','SbfI',4,60,260,5,10,'bosTau8.fasta','Bovine')

#	analyse('PstI','MspI',4,60,260,5,10,'bw6_03Feb15.fasta','Spruce budworm')
#	analyse('ApeKI','ApeKI',4,60,260,5,10,'bw6_03Feb15.fasta','Spruce budworm')
#	analyse('RsaI','RsaI',4,60,260,5,10,'bw6_03Feb15.fasta','Spruce budworm')
#	analyse('SbfI','SbfI',4,60,260,5,10,'bw6_03Feb15.fasta','Spruce budworm')

	analyse('PstI','MspI',4,80,325,5,10,'budworm.fasta','Spruce budworm')

#	analyse('PstI','MspI',4,125,325,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('ApeKI','ApeKI',4,125,325,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('RsaI','RsaI',4,125,325,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('SbfI','SbfI',4,125,325,5,10,'hs_ref_GRCh38_p7.fasta','Human')

#	analyse('PstI','MspI',4,50,100,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('ApeKI','ApeKI',4,50,100,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('RsaI','RsaI',4,50,100,5,10,'hs_ref_GRCh38_p7.fasta','Human')
#	analyse('SbfI','SbfI',4,50,100,5,10,'hs_ref_GRCh38_p7.fasta','Human')

#	analyse('PstI','MspI',4,0,100000000000,5,10,'ssa_ref_ICSASG_v2.fasta','Salmon')
#	analyse('ApeKI','ApeKI',4,0,100000000000,5,10,'ssa_ref_ICSASG_v2.fasta','Salmon')
#	analyse('SbfI','SbfI',4,0,100000000000,5,10,'ssa_ref_ICSASG_v2.fasta','Salmon')
#	analyse('RsaI','RsaI',4,0,100000000000,5,10,'ssa_ref_ICSASG_v2.fasta','Salmon')
