#!/bin/bash
# November 2018


# script for mappping aDNA reads to a divergent reference genome
# It is identical to map_SE.sh, except that the bwa mismatch is relaxed to 0.01

# this script should be run from /BEARCAVE/scripts

# example command line: ./map_SE_0.01mismatch.sh $1 $2 $3 $4
# $1 = prefix or prefixes (if you are mapping combined data file, you should enter its while set of prefixes as one argument, like this: pr1_pr2_pr3
# $2 = reference - folder name in /refgenomes/
# $3 = 3 character taxon identifier, used for filenaming - e.g. arc, spe, kud, mar.
# $4 = sample name

# variables used for mapping
MAPQ=30		# Minimum read mapping quality: 
Mismatch=0.01	# Mismatch threshold
threads=10	# no. of computer cores used by bwa and samtools
ram=16G		# GB of RAM per cpu


##### main script ######

# make processing directory
mkdir -p ../mapped$2/$1+$4_$2_map_processing
chmod u+w ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2* 2> /dev/null


# unzip trimmed reads
#gunzip ../trimdata/$1+$4_mappable.fastq.gz
zcat ../trimdata/$1+$4_mappable.fastq.gz > ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq

# Mapping: bwa aln, samse, samtools view and sort
../software/miniconda3/bin/bwa aln -n $Mismatch -t $threads ../refgenomes/$2/*.fa ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq | ../software/miniconda3/bin/bwa samse ../refgenomes/$2/*.fa - ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq | ../software/miniconda3/bin/samtools view -Su -q $MAPQ -@ $threads - | ../software/miniconda3/bin/samtools sort -@ $threads -m $ram -T ../mapped$2/$1+$4_$2_map_processing/$1+$4_temp -o ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_sorted.bam -

# zip trimmed reads
#gzip ../trimdata/$1+$4_mappable.fastq
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq

# samtools index 
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_sorted.bam

# samtools rmdup
../software/miniconda3/bin/samtools rmdup -s ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_sorted.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam

# samtools index 
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam


### mapping log ###

# generate log file
log=../mapped$2/$1+$4_$2_map_processing/"$1+$4"_"$2"_mapping.log

# gather stats
date=$(date +"%d-%m-%Y")
mappable="$(zcat ../trimdata/$1+$4_mappable.fastq.gz | wc -l | awk '{sum=$1/4; print sum}')"
mapped="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_sorted.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
uniq="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
endo="$(echo "scale=5 ; $uniq / $mappable" | bc)"
propuniq="$(echo "scale=5 ; $uniq / $mapped" | bc)"
duprate="$(echo "scale=5 ; 1 - $propuniq" | bc)"
depth="$(../software/miniconda3/bin/samtools depth ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam | awk '{sum+=$3;cnt++}END{print sum/cnt}')"
map_bp="$(../software/miniconda3/bin/samtools depth ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam | awk '{sum+=$3;cnt++}END{print sum}')"
map_Gb="$(echo "scale=5 ; $map_bp / 1000000000" | bc)"

# print to log file
printf "sample  $1+$4
reference  $2
script  map_SE_0.01mismatch.sh
date  $date
mappable_reads  $mappable
mapped_reads  $mapped
uniq_mapped_reads  $uniq
endogenous  $endo
duplication  $duprate
depth  $depth
mapped_bp  $map_bp
mapped_Gb  $map_Gb
" > $log

# remove unwanted files
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_sorted.bam*

# rename final bam file
mv ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2_$map_Gb.bam
mv ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam.bai ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2_$map_Gb.bam.bai

# set permissions
chmod 440 ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2*

# report completion
echo "Job", $1+$4, "Complete"
echo "The output can be found here: ../mapped$2/$1+$4_$2_map_processing/"
echo ''