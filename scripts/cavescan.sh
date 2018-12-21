#!/bin/bash
# October 2018

## This script will generate a comma separated table displaying the state of the cave.
## Columns 1 and 2 (Sample and Owner) contain the folders found in /rawdata/ and the owner 
## of these folders, respectively. 
## Column 3 (rawdata) lists the different prefixes of the files found in their respective
## folders, if any.
## Columns 4 to 8 (trimdata,mappedpanda,mappednewPanda,mappedpolar,mappedpseudoblack)
## state if this dataset exists in the respective folders. If files have been combined
## and therefore have several prefixes, only one of the prefixes has to be present to 
## count as "yes".
##
## When opening the table in a spread sheet viewer (such as Libreoffice Calc), make sure
## to 'format quoted fields as text', otherwise some identifiers might accidentally be 
## recognised as numbers and automatically change...

## The output table can be found in /BEARCAVE/state_of_the_cave.txt

## Usage: ./cavescan.sh

output="../state_of_the_cave.txt"

# Write the first columns of the header into output file (if the file already exists, it 
# will be overwritten).
echo -n "\"sample\",\"owner\",\"rawdata\",\"trimdata\"" > $output

# Save owners of folders in /rawdata into an array
rights=$(ls -l ../rawdata/ | grep '^d' | grep -v old_metadata | awk '{print $3}')
owners=()
for o in $rights ; do
	owners=("${owners[@]}" "$o")
done

# Save folder names in /rawdata into an array.
raw=$(ls -l ../rawdata/ | grep '^d' | grep -v old_metadata | awk '{print $9}')
raw_folders=()
for r in $raw; do
	raw_folders=("${raw_folders[@]}" "$r")
done

# Save all 'mapped' folder names into an array.
mapped=$(ls -l .. | grep '^d' | grep mapped | awk '{print $9}')
mapped_folders=()
for ma in $mapped; do
	mapped_folders=("${mapped_folders[@]}" "$ma")
	echo -n ",\"$ma\"" >> $output
done
echo '' >> $output

# Cycle through all folder names (samples) in $raw_folders.
n=0
for f in ${raw_folders[@]} ; do
	
	# Add folder name and owner to the beginning of the line.
	echo -n "\"$f\",\"${owners[$n]}\"" >> $output
	
	# Find all prefixes in the current folder. Duplicate occurrences of prefixes will be collapsed.
	rawdata=$(find ../rawdata/$f | cut -s -d'+' -f1 | rev | cut -d'/' -f1 | rev | tr "_" "\n" | sort | uniq )
	number_of_prefixes=$(echo $rawdata | wc -w)
	
	# If no prefixes exist, fill the rest of the line with "NA" and "no".
	if [ "$number_of_prefixes" == "0" ]; then
		echo -n ",\"NA\",\"no\"" >> $output
		for ma in ${mapped_folders[@]}; do
			echo -n ",\"no\"" >> $output
		done
		echo '' >> $output
		
	else
		# Cycle through all prefixes that belong to the currently processed folder name (sample).
		m=0
		for prefix in $rawdata ; do
			
			# Add the prefix into column 3 of the table.
			if (( $m==0 )) ; then
				echo -n ",\"$prefix\"" >> $output
			else
				echo -n ",,\"$prefix\"" >> $output
			fi

			# Check if prefix exists in any file name in /trimdata. And add "yes" or "no" into next column.
			in_trim=$(ls -l ../trimdata | grep $prefix | wc -l)
			if [ "$in_trim" == "0" ] ; then
				trimmed="no"
			else
				trimmed="yes"
			fi
			echo -n ",\"$trimmed\"" >> $output
			
			# Do the same for all 'mapped' folders.
			for mf in ${mapped_folders[@]}; do
				in_mapped=$(ls -l ../$mf | grep $prefix | wc -l)
				if [ "$in_mapped" == "0" ] ; then
					mapped="no"
				else
					mapped="yes"
				fi
				echo -n ",\"$mapped\"" >> $output
			done
			echo '' >> $output
			
			m=$((m+1))
		done
	fi
	n=$((n+1))
done

chmod 660 $output

echo "Done. The output table can be found in /BEARCAVE/state_of_the_cave.txt"
echo ''
