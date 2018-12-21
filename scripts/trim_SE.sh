#!/bin/bash
# November 2018

# script for trimming adapter sequences from SE data and removing short reads.

# this script should be run from /BEARCAVE/scripts

# example command line: ./trim_SE.sh $1 $2 $3
# $1 = SAMPLE*
# $2 = PREFIX*
# $3 = SEQ_RUN*

# * these can be found in /BEARCAVE/rawdata/metadata.txt

# set min length filter variable - this facilitates easy adjustment of min length threshold
minlength=30

# make processing directory in trimdata directory
mkdir -p ../trimdata/$1_processing
chmod u+w ../trimdata/$1_processing/$2+$1* 2> /dev/null

# concatenate raw fastq files
zcat ../rawdata/$1/$2+$3*R1*.fastq.gz > ../trimdata/$1_processing/$2+$1_R1.fastq

# trim adaptor seqs and short seqs from R1
../software/miniconda3/bin/cutadapt -a AGATCGGAAGAGCACACGTC -O 1 -m $minlength -o ../trimdata/$1_processing/$2+$1_mappable.fastq ../trimdata/$1_processing/$2+$1_R1.fastq > ../trimdata/$1_processing/$2+$1_trim_report.log

# zip the output
gzip ../trimdata/$1_processing/$2+$1_mappable.fastq

# clean up of concatenated raw data file
rm ../trimdata/$1_processing/$2+$1_R1.fastq

# and assign permissions
chmod 440 ../trimdata/$1_processing/$2+$1*

echo "trim of $2 complete"
echo "The output can be found here: ../trimdata/$1_processing/."
echo ''