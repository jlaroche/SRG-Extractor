#!/bin/bash

#SBATCH -D /home/jelar5/bitbucket/gbs-synthetic-genome
#SBATCH -J BuildSyntheticGenome4GBS
#SBATCH -o BuildSyntheticGenome4GBS-%j.out
#SBATCH -c 1
#SBATCH -p soyagen
#SBATCH -A soyagen
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jerome.laroche@ibis.ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=2500

./build_synthetic_GBS_genome.py
