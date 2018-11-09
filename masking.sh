#!/bin/bash

genome=$1
genome_frag_compl_bed=$2
output=$3

#bedtools maskfasta -fi Gmax_275_v2_0_no-scaffolds.fasta -bed Results.bed -fo Gmax_275_v2_0_no-scaffolds_skinny.fasta
bedtools maskfasta -fi $genome -bed $genome_frag_compl_bed -fo $output
