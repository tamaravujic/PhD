#!/bin/bash

# ==============================================================================
# Script: search_variants_in_vcf.sh
#
# Purpose:
#   Search for a specific variant across multiple VCF files listed in a master
#   file. 
#
# Input:  
#   - file_paths.txt: text file with one VCF file path per line
#   - variant: variant string to search for (e.g., "chr1:1000-2000" or "rs12345")
#
# Output: 
#   Console output showing which files contain the variant
# ==============================================================================

# Path to file containing all VCF file paths (one per line)
master_file="file_paths.txt"

# Variant or region you're searching for
variant=""

# Loop through each file path in the master file
while IFS= read -r filepath; do
	echo "Checking $filepath"
	
	# Search for variant in compressed VCF file using zgrep
	# -q flag suppresses output, only returns exit code
	if zgrep -q "$variant" "$filepath"; then
		echo "Variant found in $filepath"
	else
		echo "Variant not found in $filepath"
	fi
done < "$master_file"
