#!/bin/bash
# October 2018

# This script will combine several fastq files from ../trimdata/, also combining the file
# name prefixes. File contents and prefixes will be combined in the given order This means
# the prefixes in the resulting file name will be in the same order of sequences in the
# file. I.e., if the resulting file is called abc_def_ghi+sample.fastq.gz, then the
# sequences will be in order of abc, def, ghi.

# The output file can be found in /BEARCAVE/trimdata/.

# Usage: ./combine_files.sh $1 $2 $3 [$4 ...]
# $1 		= SAMPLE
# $2		= mappable or mappable_R1 or mappable_R2
# $3, $4...	= PREFIXES*
#
# * If one of the files is already a combination of other files, use its whole set of
# prefixes as one argument. I.e. if you want to combine the file
# abc_def+SAMPLE_mappable.fastq.gz with ghi+SAMPLE_mappable.fastq.gz, call the script like
# this: ./combine_files.sh SAMPLE mappable abc_def ghi



prefixes=''
n=3
for f in "${@:3}"; do
    if (($n == $# )); then
		prefixes="${prefixes}${f}+"
	else
		prefixes="${prefixes}${f}_"
	fi
    n=$((n+1))
done

mkdir -p ../trimdata/$prefixes$1_comb_processing/
echo -n '' > ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq.gz

for f in "${@:3}"; do
	ls -l ../trimdata/${f}+$1_$2.fastq.gz | awk '{print $9}'
  #zcat ../trimdata/${f}+$1_$2.fastq.gz >> ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq
  cat ../trimdata/${f}+$1_$2.fastq.gz >> ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq.gz
done

#echo "Zipping..."
#gzip ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq
chmod 440 ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq.gz

echo -e "Done. The output file can be found here: ../trimdata/$prefixes$1_comb_processing/$prefixes$1_$2.fastq.gz\n"



echo -e "In order to avoid redundancy and excessive disk usage, you should delete the original files.\nWould you like to do this now (enter \"1\" or \"2\")?"
select yn in "Delete" "Keep"; do
    case $yn in
        Delete )
        for g in "${@:3}"; do
        	chmod u+w ../trimdata/${g}+$1_$2.fastq.gz
        	rm ../trimdata/${g}+$1_$2.fastq.gz
		done
	    echo "Files deleted."
        break;;

        Keep )
        echo "Files kept."
        break;
    esac
done
echo ''
