#!/bin/bash
# October 2018

## This script will search for files with the given prefix in their file name
## and the corresponding entry in metadata.txt and ask the user in each case
## whether the files/entry should be deleted or not.


# usage: ./exile_from_cave.sh $1
# $1 = Prefix


if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters."

else
	echo ''
	echo '---------------------------------------'

	num_mapped=$(find ../mapped* -name "*${1}*" | wc -l)
	if [ "$num_mapped" -eq "0" ]; then
		echo "No mapping files files have been found."
	else
		echo "The following mapping files have been found:"
		find ../mapped* -name "*${1}*"
		echo "Do you want to delete them (enter \"1\" or \"2\")?"
		select yn in "Delete" "Keep"; do
			case $yn in
				Delete )
				find ../mapped* -name "*${1}*" -exec chmod u+w {} +
				find ../mapped* -name "*${1}*" -exec rm -r {} +
				echo "Mapping files deleted."
				break;;

				Keep )
				echo "Mapping files kept."
				break;
			esac
		done
	fi
	echo '---------------------------------------'
	num_trimmed=$(find ../trimdata -name "*${1}*" | wc -l)
	if [ "$num_trimmed" -eq "0" ]; then
		echo "No trimming files files have been found."
	else
		echo "The following trimming files have been found:"
		find ../trimdata -name "*${1}*"

		echo ''
		echo "Do you want to delete them (enter \"1\" or \"2\")?"
		select yn in "Delete" "Keep"; do
			case $yn in
				Delete )
				find ../trimdata -name "*${1}*" -exec chmod u+w {} +
				find ../trimdata -name "*${1}*" -exec rm -r {} +
				echo "Trimming files deleted."
				break;;

				Keep )
				echo "Trimming files kept."
				break;
			esac
		done
	fi
	echo '---------------------------------------'
	num_raw=$(find ../rawdata -name "*${1}*" | wc -l)
	if [ "$num_raw" -eq "0" ]; then
		echo "No raw files have been found."
	else
		echo "The following raw files have been found:"
		find ../rawdata -name "*${1}*"

		echo ''
		echo "Do you want to delete them (enter \"1\" or \"2\")?"
		select yn in "Delete" "Keep"; do
			case $yn in
				Delete )
				find ../rawdata -name "*${1}*" -exec chmod u+w {} +
				find ../rawdata -name "*${1}*" -exec rm -r {} +
				echo "Raw files deleted."
				break;;

				Keep )
				echo "Raw files kept."
				break;
			esac
		done
	fi
	echo '---------------------------------------'

	num_meta=$(grep $1 ../rawdata/metadata.txt | wc -l)
	if [ "$num_meta" -eq "0" ]; then
		echo "No metadata has been found."
	else
		echo "The following lines have been found in ../rawdata/metadata.txt:"
		grep $1 ../rawdata/metadata.txt

		echo ''
		echo "Do you want to delete these lines from the file (enter \"1\" or \"2\")?"
		select yn in "Delete" "Keep"; do
			case $yn in
				Delete )

				chmod u+w ../rawdata/metadata.txt
				sed -i "/$1/d" ../rawdata/metadata.txt
				chmod u-w ../rawdata/metadata.txt


				### A backup copy of the edited metadata file is being created.
				dat=$(date +%Y-%m-%d_%H.%M_)
				backup=$dat$USER
				backup+=".txt"
				chmod u+w ../rawdata/old_metadata/
				cp ../rawdata/metadata.txt ../rawdata/old_metadata/$backup
				chmod 440 ../rawdata/old_metadata/$backup
				chmod u-w ../rawdata/old_metadata/



				echo "$num_meta lines have been deleted from ../rawdata/metadata.txt"
				echo "A backup was copied to ../rawdata/old_metadata/$backup"
				break;;

				Keep )
				echo "Metadata left unchanged."
				break;
			esac
		done
	fi
	echo '---------------------------------------'



fi
