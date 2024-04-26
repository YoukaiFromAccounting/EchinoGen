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

	### This portion still needs revision to remove hard coding
	
	#Run ubiquitination site detection, save results to file
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
	#echo "Extracting sites for p58_B"
	#mkdir B_output
	#perl ESAUbiSite_main.pl ../p58_B.fasta B_output
	
	echo "Extracting UBI sites for p58_A1..."
	mkdir A1_output
	perl ESAUbiSite_main.pl ../p58_A1.fasta A1_output
	
	echo "Extracting UBI sites for p58_A2..."
	mkdir A2_output
	perl ESAUbiSite_main.pl ../p58_A2.fasta A2_output
	
	echo "Extracting UBI sites for p58_A3..."
	mkdir A3_output
	perl ESAUbiSite_main.pl ../p58_A3.fasta A3_output
	cd ..
	
	
	echo "Moving results files..."
	# Create the PTM_Results directory if it doesn't exist
	mkdir -p PTM_Results
	
	# Move the output files from glycan_detection directory to PTM_Results directory
	mv glycan_detection/*p58_A*_glycans_binary.out PTM_Results/
	
	# Move the fasta_ubicolor files from ESA-UbiSite directories to PTM_Results directory
	mv ESA-UbiSite/*/p58_A*.fasta_ubicolor PTM_Results/
	
else
    echo "Error: echino_setup.sh script failed to complete."
fi