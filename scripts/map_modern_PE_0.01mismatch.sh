#!/bin/bash
# November 2018

# script for mappping modern DNA reads to a reference genome
# note this includes mapping of both merged and unmerged paired-end reads
# It is identical to map_modern_PE.sh, except that the bwa mismatch is relaxed to 0.01


# script for mappping modern DNA reads to a reference genome
# note this includes mapping of both merged and unmerged paired-end reads

# this script should be run from /BEARVACE/scripts

# example command line: ./map_modern_PE_0.01mismatch.sh $1 $2 $3 $4
# $1 = prefix or prefixes (if you are mapping combined data file, you should enter its while set of prefixes as one argument, like this: pr1_pr2_pr3
# $2 = reference - folder name in /refgenomes/
# $3 = 3 character taxon identifier, used for filenaming - e.g. arc, spe, kud, mar.
# $4 = sample name


# variable used for mapping
MAPQ=30		# Minimum read mapping quality:
Mismatch=0.01	# Mismatch threshold
threads=10	# no. of computer cores used by bwa and samtools
ram=16G		# GB of RAM per cpu
max_insert=1000	# max library insert size


##### main script ######

# make processing directory
mkdir -p ../mapped$2/$1+$4_$2_map_processing
chmod u+w ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2* 2> /dev/null

# unzip trimmed reads
zcat ../trimdata/$1+$4_mappable.fastq.gz > ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq
zcat ../trimdata/$1+$4_mappable_R1.fastq.gz > ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R1.fastq
zcat ../trimdata/$1+$4_mappable_R2.fastq.gz > ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R2.fastq

# Mapping: bwa aln, samse, samtools view and sort
../software/miniconda3/bin/bwa aln -n $Mismatch -t $threads ../refgenomes/$2/*.fa ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq | ../software/miniconda3/bin/bwa samse ../refgenomes/$2/*.fa - ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq | ../software/miniconda3/bin/samtools view -Su -q $MAPQ -@ $threads - | ../software/miniconda3/bin/samtools sort -@ $threads -m $ram -T ../mapped$2/$1+$4_$2_map_processing/$1+$4_temp -o ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_sorted.bam -

../software/miniconda3/bin/bwa aln -n $Mismatch -t $threads ../refgenomes/$2/*.fa ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R1.fastq > ../mapped$2/$1+$4_$2_map_processing/$1+$4_r1.sai
../software/miniconda3/bin/bwa aln -n $Mismatch -t $threads ../refgenomes/$2/*.fa ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R2.fastq > ../mapped$2/$1+$4_$2_map_processing/$1+$4_r2.sai
../software/miniconda3/bin/bwa sampe -a $max_insert ../refgenomes/$2/*.fa ../mapped$2/$1+$4_$2_map_processing/$1+$4_r1.sai ../mapped$2/$1+$4_$2_map_processing/$1+$4_r2.sai ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R1.fastq ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R2.fastq  | ../software/miniconda3/bin/samtools view -Su -q $MAPQ -@ $threads - | ../software/miniconda3/bin/samtools sort -@ $threads -m $ram -T ../mapped$2/$1+$4_$2_map_processing/$1+$4_temp_PE -o ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_sorted.bam -

# remove temp data
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable.fastq
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R1.fastq
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_mappable_R2.fastq

# samtools rmdup
../software/miniconda3/bin/samtools rmdup -s ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_sorted.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_rmdup.bam
../software/miniconda3/bin/samtools rmdup ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_sorted.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_rmdup.bam

# samtools merge
../software/miniconda3/bin/samtools merge -@ $threads ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_rmdup.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_rmdup.bam

# samtools index
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_sorted.bam
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_sorted.bam
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_rmdup.bam
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_rmdup.bam
../software/miniconda3/bin/samtools index ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam


### mapping log ###

# generate log file
log=../mapped$2/$1+$4_$2_map_processing/"$1+$4"_"$2"_mapping.log

# gather stats
date=$(date +"%d-%m-%Y")
merged_mappable="$(zcat ../trimdata/$1+$4_mappable.fastq.gz | wc -l | awk '{sum=$1/4; print sum}')"
PE_mappable="$(zcat ../trimdata/$1+$4_mappable_R1.fastq.gz | wc -l | awk '{sum=$1/4; print sum}')"
merged_mapped="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_sorted.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
PE_mapped="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_sorted.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
merged_uniq="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_merged_rmdup.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
PE_uniq="$(../software/miniconda3/bin/samtools idxstats ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_PE_rmdup.bam | awk 'BEGIN {a=0} {a += $3} END{print a}')"
merged_propuniq="$(echo "scale=5 ; $merged_uniq / $merged_mapped" | bc)"
merged_duprate="$(echo "scale=5 ; 1 - $merged_propuniq" | bc)"
PE_propuniq="$(echo "scale=5 ; $PE_uniq / $PE_mapped" | bc)"
PE_duprate="$(echo "scale=5 ; 1 - $PE_propuniq" | bc)"
depth="$(../software/miniconda3/bin/samtools depth ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam | awk '{sum+=$3;cnt++}END{print sum/cnt}')"
map_bp="$(../software/miniconda3/bin/samtools depth ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam | awk '{sum+=$3;cnt++}END{print sum}')"
map_Gb="$(echo "scale=5 ; $map_bp / 1000000000" | bc)"

# print to log file
printf "sample  $1+$4
reference  $2
script  map_modern_PE_0.01mismatch.sh
date  $date
merged_mappable  $merged_mappable
PE_mappable  $PE_mappable
merged_mapped  $merged_mapped
PE_mapped  $PE_mapped
merged_uniq  $merged_uniq
PE_uniq  $PE_uniq
merged_duprate  $merged_duprate
PE_duprate  $PE_duprate
depth  $depth
mapped_bp  $map_bp
mapped_Gb  $map_Gb
" > $log

# remove unwanted files
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_r*.sai
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2*PE*.bam*
rm ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2*merged*.bam*

# rename final bam file
mv ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2_$map_Gb.bam
mv ../mapped$2/$1+$4_$2_map_processing/$1+$4_$2_rmdup.bam.bai ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2_$map_Gb.bam.bai

# set permissions
chmod 440 ../mapped$2/$1+$4_$2_map_processing/$1+$4_$3_$2*

# report completion
echo "Job", $1+$4, "Complete"
echo "The output can be found here: ../mapped$2/$1+$4_$2_map_processing/"
echo ''
