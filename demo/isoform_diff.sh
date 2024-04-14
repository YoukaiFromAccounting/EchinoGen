#!/bin/bash

#Run echino_setup.sh script
echo "Running echino_setup.sh script..."
./echino_setup.sh

#Check if echino_setup.sh script was successful
if [ $? -eq 0 ]; then
    echo "echino_setup.sh script completed successfully."

    #Run N-linked glycan site detection, save results to file
    echo "Running glycan_detection"
	cd glycan_detection || exit
	echo "Extracting NLG sites for p58_B..."
	python3 glycan.py -in ../p58_B.fasta -out B -gap 0
	echo "Extracting NLG sites for p58_A1..."
	python3 glycan.py -in ../p58_A1.fasta -out A1 -gap 0
	echo "Extracting NLG sites for p58_A2..."
	python3 glycan.py -in ../p58_A2.fasta -out A2 -gap 0
	echo "Extracting NLG sites for p58_A3..."
	python3 glycan.py -in ../p58_A3.fasta -out A3 -gap 0
	cd ..
	
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