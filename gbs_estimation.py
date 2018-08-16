#!/prg/python/3.5/bin/python3.5
"""
This script allows an one to estimatate how many fragments can be obtained with a particular pair of enzyme.  

Usage:
	./gbs_estimation.py PstI MspI 15 50 1000 5 10 Gmax_275_v2_0_no-scaffolds.fasta Soybean 1 1000 50

"""



import argparse
#parser = argparse.ArgumentParser()

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
import math,datetime,time,copy,sys
import natsort

#Enzymes choices: PstI	ApeKI

def analyse(enzymeStart,enzymeEnd,meth,minbp,maxbp,cov,mperc,genome,species,init_a_bin,max_b_bin,step_d_bin):
	print("Current date and time: " + str(datetime.datetime.now()) + "\n")

	print("€∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞ Section generale d'initiation ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞€")

	methylationPercent=int(meth)
	minBpFragments=int(minbp)
	maxBpFragments=int(maxbp)
	coverage=int(cov)
	missingPercent=int(mperc)

	dot = genome.index(".")
	ext = genome[dot+1:]
	outfile = genome[:dot]
	
# Creation des fichiers de sortie: 
	resume = open("V3_Resume_Estimate_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minbp)+"_"+str(maxbp)+".txt", 'w')
	distrib = open("V3_Distrib_Length_Estimate_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+step_d_bin+"_"+max_b_bin+".txt", 'w')

	resume.write("Species\tFile\tSequence Name\tSeq Length\tenzymeStart-enzymeEnd\tNb cutting sites\tNb fragments\tnumlim\t%Methyl\tNb fragments\tRemaining fragments\tSite mutation\tRemaining fragments\tSize limits\tNb fragments\tRemaining fragments\tCoverage (X)\tRemaining fragments\tMissing fragments\tRemaining fragments\n")
	distrib.write("Sequence\tBin\tNumber of fragments\n")
				
	print("Analysis parameters")
	print("Name of file containing genome: " + genome)
	print("Name of species: " + species)
	print("Start enzyme: " + enzymeStart)
	print("End enzyme: " + enzymeEnd)
	print("Percentage of methylation: "+ str(methylationPercent))
	print("Minimum base pair fragments: " + str(minBpFragments))
	print("Maximum base pair fragments: " + str(maxBpFragments))
	print("Average coverage: " + str(coverage))
	print("Percentage of missing: " + str(missingPercent))

#	Loading restriction enzymes
	rb = RestrictionBatch([enzymeStart])
	if(enzymeEnd != ''):
		rb.add(enzymeEnd)
#	print(rb)

# 2 RESTRICTION ENZYME USED
	if(len(rb) == 2):
		print("€∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞ Analyse a 2 sites de restriction ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞€\n")
	
#	Loading sequences in memory
		for record in SeqIO.parse(genome, "fasta"):
#			print("Analyse de la séquence " + record.name)
#			print("\n")
#			print("Total bases in sequence : " + str(len(record.seq)) + " bp")
#			print("\n")

#	Searching restriction sites
			sites = rb.search(record.seq)
#			print("\nDictionnaire 'sites' contenant les positions genomiques du/des site(s) de restriction:")
#			print(sites)

			cuts=0;
			for key, value in sites.items():
				cuts += len(value)
#			print("\nNombre de sites de restriction trouvés (cuts) "+str(cuts)+"\n")
#			print("\nNumber of cutting sites found : " + str(cuts) + "\n")

		# Since we have 2 enzymes, we merge all the positions that those 2 enzymes cuts
			allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
		# Put the positions in order
			allCuts.sort()

#			print("allCuts:\n")
#			print(allCuts)

		# We use Python sets which are more faster for searching a value in
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = set(sites[rb.get(enzymeEnd)])
		
#			print("cutsStart:\n")
#			print(cutsStart)
#			print("cutsEnd:\n")
#			print(cutsEnd)
		
#			print("\nFrom the list of all cutting sites, scan for site pair where the first site belong to startEnzyme and the second site belong to endEnzyme\n")

			i = 0
			frag = 0
			num = 1
			numlim=0
			listlength = []
			for values in allCuts:
				if(values != allCuts[-1]):     # We don't process the last cut
					if(values in cutsStart and allCuts[i+1] in cutsEnd):
						beginingPosition = values - 1
						endPosition = allCuts[i+1] - 1
						listlength.append(len(record.seq[beginingPosition:endPosition]))
						
						# Nombre de fragments total
						frag += 1

						# Nombre de fragments dans les limites demandées par l'usager: 
						if int(minBpFragments)<=len(record.seq[beginingPosition:endPosition])<=int(maxBpFragments):
							numlim+=1
						else:
							pass

				num += 1
				i += 1
#			print("Number of fragments found that begins with " + enzymeStart + " and end with " + enzymeEnd + " : " + str(frag)+"\n")

#Section pour la creation d'un tableau de distribution de frequences des longueurs de fragments
			a,b,c,d=int(init_a_bin),int(max_b_bin),int(step_d_bin),int(step_d_bin)
			# 'a' est la valeur 1 du bin
			# 'b' represente la taille max de l'avant dernier bin,
			# 'c' est la valeur 2 du bin
			# 'd' est la valeur du pas
			# 'c' doit être egal a 'd'.

			# Ici entre directement le dernier bin en fonction des longueurs supérieures a max_b_bin
			DistFreq={(b+1,max(listlength)):0}
			while c<=b:
				DistFreq[(a,c)]=0
				c+=d
				a+=d

			for values in listlength:
				for k in DistFreq:
					if (values >= k[0] and values <= k[1]):
						DistFreq[k]+=1
					else:
						pass

#			print("Voici DistFreq:")
#			print(natsort.natsorted(DistFreq.items()))
			
#Remplissage du fichier pour la distribution de frequence des longueurs de fragments
			for x in natsort.natsorted(DistFreq.items()):
#				print(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")
				distrib.write(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")

#			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

			cutsMeth = frag - ((frag * methylationPercent) / 100)
#			print("\nRemoving " + str(round((frag * methylationPercent)/100)) + " fragments according to methylation. Remaining : " + str(int(round(cutsMeth))))
#			print("\n")

			formula = -0.81 + 0.68 * (math.log(len(record.seq), 10));
			cutsFormula = cutsMeth - formula
#			print("\nRemoving " + str(round(formula)) + " fragments according to site mutation. Remaining : " + str(int(round(cutsFormula))))
#			print("\n")

			toRemove = 0
			for values in listlength:
				if(values <= minBpFragments or values >= maxBpFragments):
					toRemove += 1

			smallCuts = cutsFormula - toRemove;
#			print("\nRemoving " + str(toRemove) + " fragments less than " + str(minBpFragments) + " bp" + " and greater than " + str(maxBpFragments) + " bp" + ". Remaining : " + str(int(round(smallCuts))))
#			print("\n")

			cutsCoverage = smallCuts * coverage;
#			print("\nMultiply by depth of coverage (" + str(coverage) + "x). Remaining : " + str(int(round(cutsCoverage))))
#			print("\n")
	
			cutsMissing = cutsCoverage - ((cutsCoverage * missingPercent) / 100);
#			print("\nRemoving " + str(round(((cutsCoverage * missingPercent) / 100))) + " fragments due to missing (" + str(missingPercent) + "%). Remaining : " + str(int(round(cutsMissing))))
#			print("\n")

			resume.write(species + "\t" + genome + "\t" + record.name + "\t" + str(len(record.seq)) + "\t" + enzymeStart + "-" + enzymeEnd + "\t" + str(cuts) + "\t" + str(frag) + "\t" + str(numlim) + "\t"+ str(methylationPercent) + "\t" + str(round((frag * methylationPercent)/100)) + "\t" + str(int(round(cutsMeth))) + "\t" + str(round(formula)) + "\t" + str(int(round(cutsFormula))) + "\t" + str(minBpFragments)+"_" + str(maxBpFragments) + "\t"+ str(toRemove) + "\t" +  str(int(round(smallCuts))) + "\t" + str(coverage) + "\t" + str(int(round(cutsCoverage))) + "\t" + str(round(((cutsCoverage * missingPercent) / 100))) + "\t" + str(int(round(cutsMissing))) + "\n")

#			print("€∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞ Fin de l'analyse ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞€\n")

# 1 RESTRICTION ENZYME USED
	else:
		print("€∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞ Analyse a 1 site de restriction ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞€\n")

		for record in SeqIO.parse(genome, "fasta"):#on traite une sequence a la fois 
#			print("Analyse de la séquence " + record.name)
#			print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#			print("\n"+record.name+"\n")
#			print("Total bases in sequence : " + str(len(record.seq)) + " bp")
#			print("\n")

#	Searching restriction sites
			sites = rb.search(record.seq)
#			print("\nDictionnaire 'sites' contenant les positions genomiques du/des site(s) de restriction:")
#			print(sites)

			cuts=0;
			for key, value in sites.items():
				cuts += len(value)
#			print("\nNombre de sites de restriction trouvés (cuts) "+str(cuts)+"\n")
#			print("\nNumber of cutting sites found : " + str(cuts) + "\n")

			allCuts = sites[rb.get(enzymeStart)]
			allCuts.sort()
#			print("\nListe de toutes les positions génomiques de sites de restriction (valeur de allCuts): ")
#			print(allCuts)
#			print("\nValeur de allCuts:")
#			print(allCuts)

			i = 0
			num = 1
			numlim=0
			listlength = []
#			print("\nFrom the list of all cutting sites, scan for site pair where the first site belong to startEnzyme and the second site belong to endEnzyme\n")
			for values in allCuts:
				if(values != allCuts[-1]):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					
					listlength.append(len(record.seq[beginingPosition:endPosition]))
					
					# Nombre de fragments dans les limites demandées par l'usager:
					if int(minBpFragments)<=len(record.seq[beginingPosition:endPosition])<=int(maxBpFragments):
						numlim+=1
					else:
						pass
				i += 1
				num += 1

# Nombre de fragments total
#			print("Number of fragments found that begins with " + enzymeStart + " and end with " + enzymeEnd + " : " + str(cuts - 1)+"\n")
			frag=copy.deepcopy(cuts -1)

# Section pour la creation d'un tableau de distribution de frequences des longueurs de fragments
			a,b,c,d=int(init_a_bin),int(max_b_bin),int(step_d_bin),int(step_d_bin)
			# 'a' est la valeur 1 du bin
			# 'b' represente la taille max de l'avant dernier bin
			# 'c' est la valeur 2 du bin
			# 'd' est la valeur du pas
			# 'c' doit être egal a 'd' au départ

			#Ici entre directement le dernier bin en fonction des longueurs supérieures a max_b_bin
			DistFreq={(b+1,max(listlength)):0}
			while c<=b:
				DistFreq[(a,c)]=0
				c+=d
				a+=d

			for values in listlength:
				for k in DistFreq:
					if (values >= k[0] and values <= k[1]):
						DistFreq[k]+=1
					else:
						pass

#			print("Voici DistFreq:")
#			print(natsort.natsorted(DistFreq.items()))
			
# Remplissage du fichier pour la distribution de frequence des longueurs de fragments
			for x in natsort.natsorted(DistFreq.items()):
#				print(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")
				distrib.write(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")

#			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
			cutsMeth = frag - ((frag * methylationPercent) / 100)
#			print("\nRemoving " + str(round((frag * methylationPercent)/100)) + " fragments according to methylation. Remaining : " + str(int(round(cutsMeth))))
#			print("\n")

			formula = -0.81 + 0.68 * (math.log(len(record.seq), 10));
			cutsFormula = cutsMeth - formula
#			print("\nRemoving " + str(round(formula)) + " fragments according to site mutation. Remaining : " + str(int(round(cutsFormula))))
#			print("\n")

			toRemove = 0
			for values in listlength:
				if(values <= minBpFragments or values >= maxBpFragments):
					toRemove += 1

			smallCuts = cutsFormula - toRemove;
#			print("\nRemoving " + str(toRemove) + " fragments less than " + str(minBpFragments) + " bp" + " and greater than " + str(maxBpFragments) + " bp" + ". Remaining : " + str(int(round(smallCuts))))
#			print("\n")

			cutsCoverage = round(smallCuts) * coverage;
#			print("\nMultiply by depth of coverage (" + str(coverage) + "x). Remaining : " + str(int(round(cutsCoverage))))
#			print("\n")
	
			cutsMissing = cutsCoverage - ((cutsCoverage * missingPercent) / 100);
#			print("\nRemoving " + str(round(((cutsCoverage * missingPercent) / 100))) + " fragments due to missing (" + str(missingPercent) + "%). Remaining : " + str(int(round(cutsMissing))))
#			print("\n")

			resume.write(species + "\t" + genome + "\t" + record.name + "\t" + str(len(record.seq)) + "\t" + enzymeStart + "-" + enzymeEnd + "\t" + str(cuts) + "\t" + str(frag) + "\t" + str(numlim) + "\t"+ str(methylationPercent) + "\t" + str(round((frag * methylationPercent)/100)) + "\t" + str(int(round(cutsMeth))) + "\t" + str(round(formula)) + "\t" + str(int(round(cutsFormula))) + "\t" + str(minBpFragments)+"_" + str(maxBpFragments) + "\t"+ str(toRemove) + "\t" +  str(int(round(smallCuts))) + "\t" + str(coverage) + "\t" + str(int(round(cutsCoverage))) + "\t" + str(round(((cutsCoverage * missingPercent) / 100))) + "\t" + str(int(round(cutsMissing))) + "\n")
#			print("€∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞ Fin de l'analyse ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞€\n")

	distrib.close()
	resume.close()

#methylationPercent=int(meth)
#minBpFragments=int(minbp)
#maxBpFragments=int(maxbp)
#coverage=int(cov)
#missingPercent=int(mperc)

# Parsing user input
try:
	enzymeStart = sys.argv[1]
	enzymeEnd = sys.argv[2]
	meth = int(sys.argv[3])
	minbp = int(sys.argv[4])
	maxbp = int(sys.argv[5])
	cov = int(sys.argv[6])
	mperc = int(sys.argv[7])
	genome = sys.argv[8]
	species = sys.argv[9]
	init_a_bin = sys.argv[10]
	max_b_bin = sys.argv[11]
	step_d_bin = sys.argv[12]
except:
	print("One or more option(s) are missing.")
	print(__doc__)
	sys.exit(1)

analyse(enzymeStart,enzymeEnd,meth,minbp,maxbp,cov,mperc,genome,species,init_a_bin,max_b_bin,step_d_bin)
