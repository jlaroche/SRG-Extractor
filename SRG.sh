#!/bin/bash

red=`tput setaf 1`
blue=`tput setaf 4`
green=`tput setaf 2`
yellow=`tput setaf 3`
white=`tput setaf 7`
bold=`tput bold`
un=`tput smul`
nun=`tput rmul`

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

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
echo "Are we good to continue? y/n""$yellow"
read REP

	if [ $REP = n ]
		then
			echo "You may at first need to install the required programs"
			echo "Have a good day and see you soon""$white"
		exit 1
	elif [ $REP = y ]
		then
echo
echo
echo ""$green"Great, now please answer the following questions:""$white"
echo
echo "1- What is the name of your reference genome? (Please type the complete name like: Gmax_275_v2_0.fasta (the file should be in current working directory))""$yellow"
read REF

echo ""$green"Thank you""$white"
echo
echo "2- What is the name of your start enzyme? (Please specify like: ApeKI)""$yellow"
read enzymeStart

echo ""$green"Thank you""$white"
echo
echo "3- What is the name of your end enzyme? (If your RRS is with one restriction enzyme this will be the same as start enzyme "$un"(ApeKI)"$nun" if your RRS is with two enzymes like "$un"PstI MspI"$nun", "$un"PstI"$nun" will be the start enzyme and "$un"MspI"$nun" the end enzyme)""$yellow"
read enzymeEnd

echo ""$green"Thank you""$white"
echo
echo "4- Please define the size selection interval in bp? (e.g. 50 1000)""$yellow"
read minBpFragments maxBpFragments

echo ""$green"Thank you""$white"
echo
echo "5- What is the name of your species? (e.g. Soybean)""$yellow"
read species
echo
echo ""$green"Thank you for your responses. We are starting SRG Extactor, it may take several minutes. "$yellow"Please look SRG.log file for error reports!""$white"
echo
exec &> SRG_$TIMESTAMP.log

./srg_extractor.py $enzymeStart $enzymeEnd $minBpFragments $maxBpFragments ${REF}.fasta $species

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-1"
                                exit 1
                fi

./make_genome_file.py ${REF}.fasta

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-2"
                                exit 1
                fi

bedtools complement -i ${REF}_fragments.bed -g ${REF}.bed > ${REF}_fragments_complement.bed
	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-3"
                                exit 1
                fi

bedtools maskfasta -fi ${REF}.fasta -bed ${REF}_fragments_complement.bed -fo SRG_${REF}.fasta

	if [ $? -ne 0 ]
                        then
                                printf "There is a problem at Stage-4"
                                exit 1
                fi

echo "Your SRG is ready. Before using it, index the file:"
echo "bwa index -a bwtsw SRG_${REF}.fasta"
echo "samtools faidx SRG_${REF}.fasta"

fi

exit 0
