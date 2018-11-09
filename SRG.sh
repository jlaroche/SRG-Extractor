#!/bin/bash

red=`tput setaf 1`
blue=`tput setaf 4`
green=`tput setaf 2`
yellow=`tput setaf 3`
white=`tput setaf 7`
bold=`tput bold`
un=`tput smul`
nun=`tput rmul`

echo
echo ""$white"=============="$bold""$red"Welcome to SRG Extractor"$white"=========================="
echo
echo ""$un"Before strating with SRG Extractor, please confirm that you have already installed the following programs and they are present in PATH""$nun"
echo
echo ""$white"------------------------------------------------------------------------------------"
	echo "Linux with parallel "$blue"(http://www.gnu.org/software/parallel/)""$white"
	echo "Python 2.7 or higher "$blue"(https://www.python.org/)""$white"
	echo "Biopython "$blue"(https://biopython.org/)""$white"
	echo "samtools "$blue"(http://www.htslib.org/)""$white"
	echo "bedtools "$blue"(https://bedtools.readthedocs.io/en/latest/)""$white"
	echo "srg_extractor.py and make_genome_file.py "$blue"(these two scripts should be in current working directory)""$white"
echo ""$white"-------------------------------------------------------------------------------------"
echo
echo "Are we good to continue? y/n"
read REP

	if [ $REP = n ]
		then
			echo "You may at first need to install the required programs"
			echo "Have a good day and see you soon"
		exit 1
	elif [ $REP = y ]
		then
echo
echo
echo ""$green"Great, now please answer the following questions:""$white"
echo
echo "1- Where is your reference genome? (please type the complete path like:/home/user/refgenome/Gmax_275_v2_0.fasta"
read REF

echo ""$green"Thank you""$white"
echo
echo "2- What is your start enzyme? (Please specify like: ApeKI)"
read enzymeStart

echo ""$green"Thank you""$white"
echo
echo "3- What is your end enzyme? (This will be the same as start enzyme (ApeKI) if your RRS is with one enzyme)"
read enzymeEnd

echo ""$green"Thank you""$white"
echo
echo "4- Please define the size selection interval in bp? (e.g. 50 1000)"
read minBpFragments maxBpFragments

echo ""$green"Thank you""$white"
echo
echo "5- What is the name of your species? (e.g. Soybean)"
read species
echo
echo ""$green"Thank you for your responses. We are starting SRG Extarctor, it may take several minutes.""$white"
echo
exec &> SRG.log

./srg_extractor.py $enzymeStart $enzymeEnd $minBpFragments $maxBpFragments $REF $species

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-1"
                                exit 1
                fi

./make_genome_file.py $REF

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-2"
                                exit 1
                fi

bedtools complement -i genome_fragments.bed -g genome.bed > genome_fragments_complement.bed
	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-3"
                                exit 1
                fi

bedtools maskfasta -fi $REF -bed genome_fragments_complement.bed -fo SRG_$REF.fasta

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-4"
                                exit 1
                fi

echo "Your SRG is ready, please before using index your SRG"
echo ""$green"Good bye""$white"
fi
exit 0
