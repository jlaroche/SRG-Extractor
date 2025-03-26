#!/bin/bash

genome_frag_bed=$1
genome_bed=$2
genome_frag_compl_bed=$3

bedtools complement -i $genome_frag_bed -g $genome_bed > ${genome_frag_compl_bed}

