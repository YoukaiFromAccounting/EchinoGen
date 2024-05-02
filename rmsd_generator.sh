#!/bin/bash

# Input file containing URLs
input_file="isoforms.txt"
# Store all isoforms in array file names to be accessed for subsiquent use
file_names=()

# Read each line from the file and extract the value after the last equal sign to represent protein name
while IFS= read -r url; do
    # Extract the name of the protein from the URL
    URL="$url"
    name=$(echo "$url" | grep -oE '[^=]+$')
    file_names+=("$name")
done < "$input_file"


for file1 in "${file_names[@]}"; do
    echo "$file1"
    for file2 in "${file_names[@]}"; do
        echo "$file2"
        if [ "$file1" != "$file2" ]
        then
            echo Running rmsd analysis on "$file1" and "$file2"
            python3 rmsd_plot.py "$file1".pdb "$file2".pdb  
        fi
    done
done
