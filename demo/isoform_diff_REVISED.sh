#!/bin/bash

#Run echino_setup.sh script
echo "Running echino_setup.sh script..."
./echino_setup_REVISED.sh

#### Read in the URLS and run wget requests from here
#### Save the protein names in an array to access later when calling glycan detection etc. 

# Input file containing URLs
input_file="isoforms.txt"
# Store all isoforms in array file names to be accessed for subsiquent use
file_names=()

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file not found!"
    exit 1
fi

# Read each line from the file and extract the value after the last equal sign to represent protein name
while IFS= read -r url; do
    # Extract the name of the protein from the URL
    URL="$url"
    name=$(echo "$url" | grep -oE '[^=]+$')
    file_names+=("$name")
    # run wget request from the url and name the fasta file the protein name extracted from the URL
    wget -O "$name".fasta "$url"
    echo "Downloading fasta file "$name""
done < "$input_file"

### Revise the rest below to instead to loop through array of isoform names and run glycan detection that way
 
#Check if echino_setup.sh script was successful
if [ $? -eq 0 ]; then
    echo "echino_setup.sh script completed successfully."

	cd glycan_detection || exit

	# loop through file names and run glycan detection on each fasta file
	for file in "$filenames"; do
		python3 glycan.py -in ../"$file".fasta -out "$file" -gap 0
	done

	
	# Loop through file names and run ubiquitination site detection, save results to file for each fasta file
	cd ..
	echo "Building AAIndex..."
	cd ESA-UbiSite
	cd src/aaindex
	make
	cd ..
	cd ..
	echo "Building LIBSVM model..."
	cd src/libsvm_320
	make
	cd ..
	cd ..
	
	for file in "$filenames"; do
		mkdir "$file"_output
		perl ESAUbiSite_main.pl ../"$file".fasta "$file"_output
	done
	
	echo "Moving results files..."
	# Create the PTM_Results directory if it doesn't exist
	mkdir -p PTM_Results
	
	# Move the output files from glycan_detection directory to PTM_Results directory
	mv glycan_detection/*_glycans_binary.out PTM_Results/
	
	# Move the fasta_ubicolor files from ESA-UbiSite directories to PTM_Results directory
	mv ESA-UbiSite/*/*.fasta_ubicolor PTM_Results/

	# Run the rmsd_plot.py on each file in file_names 

	length=${#file_names[@]}

	for ((i=1; i<"$length"; i)) do
		python3 rmsd_plot.py ${file_names[$i]} ${file_names[$i-1]}
		#### THIS STEP IS LEADING TO AN INFINITE LOOP CURRENTLY
		# mv *.png PTM_Results/
	done

	
else
    echo "Error: echino_setup.sh script failed to complete."
fi