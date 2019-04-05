#!/bin/bash
### October 2018

###
### Usage: ./add_to_meta.sh inputfile.txt
###
### The input file has to be of a certain format: It has to be a comma separated table
### and the header line (column names) have to be the same as in the metadata.txt, except
### for the first column (PREFIX), like so:
### SEQ_RUN,SAMPLE,TAXON,LOCALITY,COUNTRY,AGE,DATABASE_NO,LIBRARY_NO,EXTRACT_METH,LIBRARY_METH,PLATFORM,READ_LENGTH,MODE,SEQ_PRIMER,RIGHTS
### Using special characters in the input file should mostly be ok, although you should
### double-check the outcome if you use anything beyond the ordinary.
### Things you should definitely not use are commas, semicolons and backslashes!
### Missing data should be entered as NA.
###
### The rawdata file names must not contain a + (plus) and are automatically renamed
### to get a random three-character prefix (base36 number), followed by a +.
###
### A backup copy of the edited metadata file is being stored under
### ../rawdata/old_metadata/YEAR-MONTH-DAY_HOUR.MINUTE_USER.txt

### This script grew in length and complexity over time which means it works very
### inefficiently. Adding many datasets at once will take time.
### Actually, this script should be rewritten from scratch...

# chmod u+w ../rawdata/metadata.txt

lines_before=$(wc -l ../rawdata/metadata.txt | cut -d ' ' -f 1)
echo -e "\nmetadata.txt contained $lines_before lines.\n"

######### This block is checking the input file for integrity. #######################
n=1
checked_lines=()
stripped_meta_file=$(cat ../rawdata/metadata.txt | sed -e 's/"/''/g' | sed -e 's/\x27/''/g' | cut -d',' -f2-16)
while IFS='' read -r line || [[ -n "$line" ]]; do # All that fancy stuff is supposed to avoid trouble when special characters occur.
	separators=$(grep -o "," <<< "$line" | wc -l)
    ## The first line of the file (column names) has to match the format of the metadata.txt file, without the first column.
    stripped_line=$(echo $line | sed -e 's/"/''/g' | sed -e 's/\x27/''/g')
    if (($n == 1)); then
    	if [ "$stripped_line" != 'SEQ_RUN,SAMPLE,TAXON,LOCALITY,COUNTRY,AGE,DATABASE_NO,LIBRARY_NO,EXTRACT_METH,LIBRARY_METH,PLATFORM,READ_LENGTH,MODE,SEQ_PRIMER,RIGHTS' ]; then
    		echo "ERROR: Column names do not match. Metadata file has not been changed."
    		broken=true
    		break
    	fi
    ## Empty lines are being skipped.
    elif [ "$line" == "" ]; then
   	 	n=$((n+1))
    	continue
	## This checks if the lines contain the appropriate amount of values, based on the number of commas and semicolons.
	elif (($separators != 14)); then
		echo "ERROR: Line $n does not contain the appropriate amount of values. Metadata file has not been changed."
		broken=true
		break
	elif echo "$stripped_meta_file" | grep -q "$stripped_line"; then
		echo "ERROR: Line $n is already contained in metadata.txt. Metadata file has not been changed."
		broken=true
		break
    else
        checked_lines=("${checked_lines[@]}" "$stripped_line")
    fi
    n=$((n+1))
done < "$1"
############################################################################################


## If no errors occurred, the input file is added to the metadata file, ignoring the header line and empty lines.
if ! [ "$broken" = true ]; then
	m=0
	while (( $m < ${#checked_lines[@]} )); do


		######### Random three-characters prefix is generated here. ####################################
		chars=0123456789abcdefghijkLmnopqrstuvwxyz
		while true; do
			## Here, a random three-character prefix is generated from numbers and letters.
			prefix=''
			unique=''
			for i in {1..3} ; do
	    		prefix=$prefix${chars:RANDOM%${#chars}:1}
			done
			## This checks if the generated prefix already exists in the first column of the metadata.txt
			## If it exists, the "unique" variable is set to false and therefore the above while loop
			## continues, generating a new prefix that is again checked.
			while IFS='' read -r metadata_line ; do
				meta_entry=$(echo $metadata_line | cut -d',' -f1 | sed -e 's/"/''/g')
				if [ "$meta_entry" == $prefix ]; then
					unique=false
					echo "..."
					break
				fi
			done < ../rawdata/metadata.txt
			if ! [ "$unique" = false ]; then
				break
			fi
		done
		###################################################################################################

		###################### File renaiming ############################
		sample_name=$(echo "${checked_lines[$m]}" | cut -d',' -f2)
		seqrun=$(echo "${checked_lines[$m]}" | cut -d',' -f1)

		if cd ../rawdata/$sample_name; then
			for f in $seqrun*.gz; do

				## If one of the file names already containes a + the process is being aborted.
				pluses=$(ls -1 | grep $seqrun | grep '+' | wc -l)
				if ! [ "$pluses" == "0" ] ; then
					echo -e "ERROR: At least one file name of $seqrun contains a +. This is not allowed. A possible reason might be that $seqrun is already in the system.\nProcess aborted in line $((m+1)) of the input file.\nPlease carefully check what went wrong and view the ../rawdata/metadata.txt to see which data has been incorporated before you add the rest (possibly with a new input file).\n"
					broken=true
					break
				fi
				## If an error occurs during the file renameing, the process is being aborted.
				if ! [ "$broken" = true ]; then
					if ! mv "$f" "$prefix+$f" ; then
						echo -e "ERROR occurred durung file renaming. Process aborted in line $((m+1)) of the input file. \nPlease carefully check what went wrong and view the ../rawdata/metadata.txt to see which data has been incorporated before you add the rest (possibly with a new input file).\n"
						broken=true
						break
					fi
					echo "$f has been renamed to $prefix+$f"
					chmod 440 $prefix+$f
				fi
			done
			cd ../../scripts/
		else
			echo "ERROR: Folder ../rawdata/$sample_name not found."
			broken=true
			break
		fi
		echo ""
		#################################################################


		if ! [ "$broken" = true ]; then
		    quoted_line=$(echo ${checked_lines[$m]} | sed 's/^/"/' | sed 's/$/"/' | sed 's/,/","/g')

		    cp ../rawdata/metadata.txt ../rawdata/temp.meta
		    chmod u+w ../rawdata/temp.meta
		    echo "\"$prefix\",$quoted_line" >> ../rawdata/temp.meta
		    mv -f ../rawdata/temp.meta ../rawdata/metadata.txt
		    chmod u-w ../rawdata/metadata.txt

    		#echo "\"$prefix\",$quoted_line" >> ../rawdata/metadata.txt
		fi

		m=$((m+1))
	done

	if ! [ "$broken" = true ]; then
		### A backup copy of the edited metadata file is being created.
		dat=$(date +%Y-%m-%d_%H.%M_)
		backup=$dat$USER
		backup+=".txt"

		cp ../rawdata/metadata.txt ../rawdata/old_metadata/$backup
		chmod 440 ../rawdata/old_metadata/$backup
	### Info message.
		echo -e "${#checked_lines[@]} lines were added to metadata.txt.\n"
		lines_after=$(wc -l ../rawdata/metadata.txt | cut -d ' ' -f 1)
		echo -e "metadata.txt now contains $lines_after lines.\n"
		echo -e "A backup copy was saved to ../rawdata/old_metadata/$backup\n"
	fi

fi

chmod u-w ../rawdata/metadata.txt
