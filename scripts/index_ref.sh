#!/bin/bash
# November 2018

# This Script will index the given reference genome with BWA and SAMtools and generate
# a "mapped" folder in /BEARCAVE/.

# This script should be run from /BEARCAVE/scripts
# Example command line: sh index_ref.sh $1
# $1 = 	Fasta file containing the reference genome,  
# 		e.g. ../refgenomes/PolarBear_mt/polarbear_mt.fasta.gz

softdir=../software/miniconda3/envs/py35/bin

folder=$(echo $1 | rev | cut -d'/' -f2 | rev)

gunzip $1 2> /dev/null

file_name=$(echo $1 | sed 's/.gz$//' )
file_ext=$(echo $file_name | rev | cut -d'.' -f1 | rev)
file_core=$(echo $file_name | rev | cut --complement -d'.' -f1 | rev)

if [ $file_ext = "fasta" ] || [ $file_ext = "FASTA" ]; then
	fa=$(echo $file_core.fa)
	mv $file_name $fa
else
	fa=$file_name
fi


$softdir/bwa index $fa
$softdir/samtools faidx $fa
chmod 550 ../refgenomes/$folder/
chmod 440 ../refgenomes/$folder/*

mkdir -p ../mapped$folder
mkdir -p ../mapped$folder/${folder}_logs
chmod 770 -R ../mapped$folder