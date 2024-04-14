#!/bin/bash

#Ensure available memory is above 40GB for running of protein prediction pipeline
echo "Checking free memory"
available_memory=$(free -h | awk '/^Mem:/{print $2}')
if [[ $available_memory =~ ([0-9.]+)Gi ]]; then
    memory_gb=${BASH_REMATCH[1]}
    if (( $(echo "$memory_gb < 40" | bc -l) )); then
        echo "Warning: Available memory is below 40 GB ($memory_gb GiB)"
    fi
fi

#Check available cores
echo "Checking available cores"
num_cores=$(nproc)
if ((num_cores < 4)); then
    echo "Warning: Number of CPU cores is less than 4 ($num_cores cores)"
fi

#Download glycan_detection software NGLYCOSYLATION
echo "Downloading glycan_detection software..."
git clone https://github.com/davidmatten/glycan_detection
#cd glycan_detection || exit

#Download ubiquitination detection software UBIQUITINATION
echo "Downloading ubiquitination_detection software..."
git clone https://github.com/NctuICLab/ESA-UbiSite

#Download acetylation site detection software ACETYLATION

#Download phosformer software PHOSPHORYLATION
echo "Downloading phosformer software..."
git clone https://github.com/waylandy/phosformer

#Create and activate conda environment
echo "Creating conda environment..."
#conda env create -f environment.yaml
echo "Activating conda environment..."
#conda activate phosformer

#Check for required Python packages and install if not present
echo "Checking for required Python packages..."
if ! python -c 'import pandas, argparse, os, regex, smallBixTools, numpy, scikit-learn, tensorflow, matplotlib' &> /dev/null; then
    echo "Installing required Python packages..."
    pip install pandas argparse regex smallBixTools numpy scikit-learn tensorflow matplotlib
fi

#Download sea urchin protein fasta files from NCBI 
#Rename to p58_B
wget -O p58_B.fasta 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=fasta&sort=&val=XP_799905.2'
#Rename to p58_A1
wget -O p58_A1.fasta 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=fasta&sort=&val=XP_003724499.1'
#Rename to p58_A2
wget -O p58_A2.fasta 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=fasta&sort=&val=72005711'
#Rename to p58_A3
wget -O p58_A3.fasta 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=fasta&sort=&val=XP_003724498.1'

#Download sample protein fasta file 
echo "Downloading sample protein fasta file..."
wget -O sample_protein.fasta 'https://www.uniprot.org/uniprot/P12345.fasta'

#Run glycan_detection with the sample protein fasta file
echo "Running glycan_detection with sample protein fasta file..."
./glycan_detection sample_protein.fasta

echo "Setup complete."