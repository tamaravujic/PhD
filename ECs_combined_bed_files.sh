#!/bin/bash

# Define the main directory
main_dir="ECs_ATAC"

# Change directory to the main directory
cd "$main_dir" 

# Define the output file
output_file="ECs_combined_enhancer_peaks.bed"


# Iterate over each folder in the main directory
for folder in */; do
    folder_name=$(basename "$folder")

    # Define the input file in the current folder
    input_file="${folder}${folder_name}_enhancer_peaks_named.bed"
    
    # Check if the input file exists
    if [ -f "$input_file" ]; then
        # Append contents of the input file to the output file
        cat "$input_file" >> "$output_file"
    fi
done
