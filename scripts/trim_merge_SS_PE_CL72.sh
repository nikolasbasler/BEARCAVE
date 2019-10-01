#!/bin/bash
# November 2018

# script for trimming adapter sequences from SS library PE data, removing short reads, amd merging overlapping PE reads.

# this script should be run from /BEARCAVE/scripts

# example command line: ./trim_merge_SS_PE_CL72.sh $1 $2 $3
# $1 = SAMPLE*
# $2 = PREFIX*
# $3 = SEQ_RUN*

# * these can be found in /REPOSITORY/rawdata/metadata.txt
softdir=../software/miniconda3/envs/py35/bin
# set min length filter variable - this facilitates easy adjustment of min length threshold
# set threads variable for number of cpus used by flash
minlength=30
threads=10

# make processing directory in trimdata directory
mkdir -p ../trimdata/$1_processing
chmod u+w ../trimdata/$1_processing/$2+$1* 2> /dev/null

# concatenate raw fastq files
zcat ../rawdata/$1/$2+$3*R1*.fastq.gz > ../trimdata/$1_processing/$2+$1_R1.fastq
zcat ../rawdata/$1/$2+$3*R2*.fastq.gz > ../trimdata/$1_processing/$2+$1_R2.fastq

echo "Script: trim_merge_SS_PE_CL72.sh" > ../trimdata/$1_processing/$2+$1_trim_report.log
echo "Script: trim_merge_SS_PE_CL72.sh" > ../trimdata/$1_processing/$2+$1_merge_report.log


# trim adaptor seqs and short seqs from R1 and R2
$softdir/cutadapt -a AGATCGGAAGAGCACACGTC -A GGAAGAGCGTCGTGTAGGGA -O 1 -m $minlength -o ../trimdata/$1_processing/$2+$1_trim_R1.fastq -p ../trimdata/$1_processing/$2+$1_trim_R2.fastq ../trimdata/$1_processing/$2+$1_R1.fastq ../trimdata/$1_processing/$2+$1_R2.fastq >> ../trimdata/$1_processing/$2+$1_trim_report.log

# merge R1 and R2, set max overlap (-M) to 75bp as most ancient frags should be mergeable
$softdir/flash -M 75 -t $threads -d ../trimdata/$1_processing -o $2+$1 ../trimdata/$1_processing/$2+$1_trim_R1.fastq ../trimdata/$1_processing/$2+$1_trim_R2.fastq >> ../trimdata/$1_processing/$2+$1_merge_report.log

# clean up unnecessary files
# here I am assuming ancient DNA data and discarding the non-overlapping PE reads
mkdir ../trimdata/$1_processing/$2+$1_save
mv ../trimdata/$1_processing/$2+$1_trim_report.log ../trimdata/$1_processing/$2+$1_merge_report.log ../trimdata/$1_processing/$2+$1.extendedFrags.fastq ../trimdata/$1_processing/$2+$1_save
rm ../trimdata/$1_processing/$2+$1* 2> /dev/null
mv ../trimdata/$1_processing/$2+$1_save/$2+$1* ../trimdata/$1_processing/
rmdir ../trimdata/$1_processing/$2+$1_save

# rename merged reads file
mv ../trimdata/$1_processing/$2+$1.extendedFrags.fastq ../trimdata/$1_processing/$2+$1_mappable.fastq

# and zip it
gzip ../trimdata/$1_processing/$2+$1_mappable.fastq

# and assign permissions
chmod 440 ../trimdata/$1_processing/$2+$1*

echo "trim and merge of $2+$1 complete"
echo "The output can be found here: ../trimdata/$1_processing/."
echo ''
