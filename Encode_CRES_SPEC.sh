#!/bin/bash

# Define the main directory
main_dir="ENCODE_ATAC_COPY"

# Change directory to the main directory
cd "$main_dir" 

# Iterate over each folder in the main directory
for folder in */; do
    folder_name=$(basename "$folder")

    # Combine rep1 and rep2 broadPeak files
    cat "${folder}rep1_sorted_peaks.broadPeak" "${folder}rep2_sorted_peaks.broadPeak" > "${folder}${folder_name}_combined.broadPeak"

    # Print the cell name in a new column
    awk -F'\t' -v folder_name="$folder_name" '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" folder_name}' "${folder}${folder_name}_combined.broadPeak" > "${folder}${folder_name}_named_combined.bed"

    # Intersect bed with the ERG_CREs.bed file
    intersectBed -c -wa -a "../ERG_CREs.bed" -b "${folder}${folder_name}_named_combined.bed" > "${folder}${folder_name}_enhancer_peaks.bed"

    awk -F'\t' -v folder_name="$folder_name" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" folder_name "\t" $5}' "${folder}${folder_name}_enhancer_peaks.bed" > "${folder}${folder_name}_enhancer_peaks_named.bed"

done

