#!/bin/bash
# November 2018
softdir=../software/miniconda3/envs/py35/bin
# script for trimming adapter sequences from SE data and removing short reads of variable length.

# this script should be run from /BEARCAVE/scripts

# example command line: ./trim_SE_var_min_length.sh $1 $2 $3 $4
# $1 = SAMPLE*
# $2 = PREFIX*
# $3 = SEQ_RUN*
# $4 = min length threshold to be used
# * these can be found in /BEARCAVE/rawdata/metadata.txt


# make processing directory in trimdata directory
mkdir -p ../trimdata/$1_$4bp_processing
chmod u+w ../trimdata/$1_$4bp_processing/$2+$1* 2> /dev/null

# concatenate raw fastq files
zcat ../rawdata/$1/$2+$3*R1*.fastq.gz > ../trimdata/$1_$4bp_processing/$2+$1_R1.fastq

echo "Script: trim_SE_var_min_length.sh - Minimum read length: $4 bp" > ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_report.log

# trim adaptor seqs and short seqs from R1
$softdir/cutadapt -a AGATCGGAAGAGCACACGTC -O 1 -m $4 -o ../trimdata/$1_$4bp_processing/$2+$1_$4bp_mappable.fastq ../trimdata/$1_$4bp_processing/$2+$1_R1.fastq >> ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_report.log

# zip the output
gzip ../trimdata/$1_$4bp_processing/$2+$1_$4bp_mappable.fastq

# clean up of concatenated raw data file
rm ../trimdata/$1_$4bp_processing/$2+$1_R1.fastq

# and assign permissions
chmod 440 ../trimdata/$1_$4bp_processing/$2+$1*

echo "trim of $2 min length $4 bp complete"
echo "The output can be found here: ../trimdata/$1_$4bp_processing."
echo ''
