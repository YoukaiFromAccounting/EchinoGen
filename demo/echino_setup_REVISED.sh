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


echo "Setup complete."