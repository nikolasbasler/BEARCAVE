#!/bin/bash
# November 2018

# script for trimming adapter sequences from SS library PE data, removing short reads of user defined length, and merging overlapping PE reads.

# this script should be run from /BEARCAVE/scripts

# example command line: ./trim_merge_SS_PE_CL72_var_min_length.sh $1 $2 $3 $4
# $1 = SAMPLE*
# $2 = PREFIX*
# $3 = SEQ_RUN*
# $4 = min length threshold to be used
# * these can be found in /REPOSITORY/rawdata/metadata.txt

# set threads variable for number of cpus used by flash
threads=10


# make processing directory in trimdata directory
mkdir -p ../trimdata/$1_$4bp_processing
chmod u+w ../trimdata/$1_$4bp_processing/$2+$1* 2> /dev/null

# concatenate raw fastq files
zcat ../rawdata/$1/$2+$3*R1*.fastq.gz > ../trimdata/$1_$4bp_processing/$2+$1_$4bp_R1.fastq
zcat ../rawdata/$1/$2+$3*R2*.fastq.gz > ../trimdata/$1_$4bp_processing/$2+$1_$4bp_R2.fastq

# trim adaptor seqs and short seqs from R1 and R2
../software/miniconda3/bin/cutadapt -a AGATCGGAAGAGCACACGTC -A GGAAGAGCGTCGTGTAGGGA -O 1 -m $4 -o ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_R1.fastq -p ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_R2.fastq ../trimdata/$1_$4bp_processing/$2+$1_$4bp_R1.fastq ../trimdata/$1_$4bp_processing/$2+$1_$4bp_R2.fastq > ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_report.log

# merge R1 and R2, set max overlap (-M) to 75bp as most ancient frags should be mergeable
../software/miniconda3/bin/flash -M 75 -t $threads -d ../trimdata/$1_$4bp_processing -o $2+$1_$4bp ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_R1.fastq ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_R2.fastq > ../trimdata/$1_$4bp_processing/$2+$1_$4bp_merge_report.log

# clean up unnecessary files
# here I am assuming ancient DNA data and discarding the non-overlapping PE reads 
mkdir ../trimdata/$1_$4bp_processing/$2+$1_$4bp_save
mv ../trimdata/$1_$4bp_processing/$2+$1_$4bp_trim_report.log ../trimdata/$1_$4bp_processing/$2+$1_$4bp_merge_report.log ../trimdata/$1_$4bp_processing/$2+$1_$4bp.extendedFrags.fastq ../trimdata/$1_$4bp_processing/$2+$1_$4bp_save
rm ../trimdata/$1_$4bp_processing/$2+$1_$4bp* 2> /dev/null
mv ../trimdata/$1_$4bp_processing/$2+$1_$4bp_save/$2+$1_$4bp* ../trimdata/$1_$4bp_processing/
rmdir ../trimdata/$1_$4bp_processing/$2+$1_$4bp_save

# rename merged reads file 
mv ../trimdata/$1_$4bp_processing/$2+$1_$4bp.extendedFrags.fastq ../trimdata/$1_$4bp_processing/$2+$1_$4bp_mappable.fastq

# and zip it
gzip ../trimdata/$1_$4bp_processing/$2+$1_$4bp_mappable.fastq

# and assign permissions
chmod 440 ../trimdata/$1_$4bp_processing/$2+$1*


echo "trim of $2+$1 min length $4 bp complete"
echo "The output can be found here: ../trimdata/$1_$4bp_processing/."
echo ''