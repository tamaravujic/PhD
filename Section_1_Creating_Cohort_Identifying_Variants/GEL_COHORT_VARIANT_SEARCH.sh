#!/bin/bash

# ==============================================================================
# Script: PL_COHORT_VARIANT_SEARCH.sh
#
# Purpose:
#   Intersect VCF variant files with genomic coordinate of Cis regulatory elements (BED format).
#   Extracts variants falling within specified enhancer or regulatory regions
#   for downstream analysis. Processes multiple VCF files in batch mode.
#
# Workflow:
#   1. Read list of VCF file paths from input file
#   2. For each VCF file:
#      - Extract unique sample identifier from file path
#      - Intersect VCF with coordinate file using bedtools
#      - Save intersected variants to output file
#      - Report completion for each file
#
# Input Files:
#   vcf_paths.txt - text file with one VCF file path per line
#   Enhancer_regions.txt - BED-format file with genomic coordinates
#
# Output:
#   {unique_id}_variants.txt - variants falling within coordinate regions
#
# Prerequisites:
#   - BEDTools installed and in PATH
#   - Input VCF files in gzip-compressed format (.vcf.gz)
#   - BED file with enhancer or regulatory region coordinates
#
# Notes:
#   - Extracts unique sample identifier from 7th directory level of file path
#   - Creates one output file per VCF file processed
#   - Skips missing files with warning message (does not stop execution)
# ==============================================================================

# ==============================================================================
# SECTION 1: DEFINE INPUT FILE PATHS
# ==============================================================================
# Specify locations of the VCF file list and coordinate regions file


# File containing paths to .vcf.gz files (one path per line)
VCF_PATHS="vcf_paths.txt"

# File containing BED-format coordinates for intersecting
# Format: chromosome, start, end
COORDINATES="Enhancer_regions.txt"

# ==============================================================================
# SECTION 2: PROCESS EACH VCF FILE
# ==============================================================================
# Loop through each line in the VCF paths file

while IFS= read -r vcf_file || [[ -n "$vcf_file" ]]; do
    
    # Check that the VCF file exists before processing
    # Skip this iteration if file is missing (with warning)
    if [[ ! -f "$vcf_file" ]]; then
        echo "Warning: File $vcf_file not found, skipping..."
        continue
    fi

    # ==============================================================================
    # SECTION 3: EXTRACT UNIQUE SAMPLE IDENTIFIER
    # ==============================================================================
    # Extract the 7th component of the file path as the unique identifier
    # Example: /path/to/samples/patient123/vcf/file.vcf.gz → patient123
    # This identifier is used to name the output file uniquely
    unique_id=$(echo "$vcf_file" | awk -F '/' '{print $7}')

    # ==============================================================================
    # SECTION 4: DEFINE OUTPUT FILE NAME
    # ==============================================================================
    # Create output filename based on the extracted unique identifier
    # All output files follow the naming convention: {unique_id}_variants.txt
    output_file="${unique_id}_variants.txt"

    # ==============================================================================
    # SECTION 5: INTERSECT VCF WITH COORDINATE REGIONS
    # ==============================================================================
    # Use bedtools intersect to find all variants within the coordinate regions
    # -a: query file (VCF file with variants)
    # -b: reference file (BED file with regions)
    # -header: preserve header lines from VCF file
    # Output: variants that overlap with enhancer regions
    bedtools intersect -a "$vcf_file" -b "$COORDINATES" -header > "$output_file"

    # Report completion of this file's processing
    echo "Variants for $vcf_file written to $output_file"

done < "$VCF_PATHS"
