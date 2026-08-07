#!/usr/bin/env bash

# ==============================================================================
# Script: GEL_find_SV_breakpoints.sh
#
# Purpose:
#   Identifies structural variants (SVs) with breakpoints located within
#   CRE regions and ERG gene body. Unlike finding overlapping SVs, this script specifically 
#   locates variants where the breakpoint coordinates
#   fall within enhancer boundaries.
#
# Workflow:
#   1. Extract all structural variants from VCF files
#   2. For each variant, check if start AND end positions overlap  regions
#   3. Output variants with breakpoints in enhancers to TSV file
#
# Input:
#   VCF files with structural variant calls (gzip-compressed)
#   BED file with CRE region coordinates (chromosome, start, end)
#
# Output:
#   ERG_breakpoints.tsv - tab-delimited file with breakpoint information
#     Columns: Sample ID, chromosome, position, end, variant type, ID
#
#
# Notes:
#   - Checks if variant start position AND end position both overlap BED regions
#   - Only retains variants where breakpoints fall within enhancer coordinates
# ==============================================================================

#BSUB -q long
#BSUB -P re_gecip_cardiovascular
#BSUB -cwd "/re_gecip/cardiovascular/tvujic/SVs/hg38"
#BSUB -o SVs_job_.stdout
#BSUB -e SVs_job_.stderr
#BSUB -n 1
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000

# ==============================================================================
# SECTION 1: SET ERROR HANDLING AND LOAD MODULES
# ==============================================================================
# Enable pipefail to catch errors in piped commands
set -euo pipefail

# Load required bioinformatics tools
module load bcftools
module load bedtools

# ==============================================================================
# SECTION 2: DEFINE INPUT AND OUTPUT FILES
# ==============================================================================
# Specify paths to input VCF list, CRE regions, and output file

# Text file containing paths to structural variant VCF files (one per line)
VCF_LIST=""

# BED file containing CRE region coordinates
# Format: chromosome, start, end
BED_FILE=""

# Output file for variants with breakpoints in enhancers
OUT=""

# ==============================================================================
# SECTION 3: CREATE OUTPUT FILE AND TEMPORARY FILE
# ==============================================================================
# Columns: Sample ID, chromosome, position (start), end, variant type, ID
# Remove any existing temporary file to start fresh

echo -e "Sample\tCHROM\tPOS\tEND\tSVTYPE\tID" > "$OUT"
rm -f tmp_all.tsv

# ==============================================================================
# SECTION 4: PROCESS EACH VCF FILE
# ==============================================================================
# Loop through each VCF file in the input list

while read -r VCF; do
    # Skip empty lines in the VCF list
    [[ -z "$VCF" ]] && continue
    
    # Check if VCF file exists and is non-empty
    if [[ ! -s "$VCF" ]]; then
        echo "SKIP empty or missing: $VCF"
        continue
    fi
    
    # Report progress
    echo "Processing $VCF"
    
    # Extract sample ID using bcftools query
    # bcftools query -l lists all sample names in the VCF
    SAMPLE=$(bcftools query -l "$VCF")
    
    # ==============================================================================
    # SECTION 5: EXTRACT ALL STRUCTURAL VARIANTS FROM VCF
    # ==============================================================================
    # Use bcftools query to extract structural variants with their coordinates
    # Output format: SAMPLE, chromosome, position (start), end, variant type, ID
    # Store in temporary file for filtering
    
    bcftools query \
        -f "${SAMPLE}\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\n" \
        "$VCF" >> tmp_all.tsv

done < "$VCF_LIST"

# ==============================================================================
# SECTION 6: FILTER BREAKPOINTS - BED REGION OVERLAP CHECK
# ==============================================================================
# Use awk to safely filter variants by checking if breakpoint positions
# fall within BED region coordinates


awk 'NR==FNR {
    # Load BED file regions: store each region indexed by line number
    a[i]=$1 "\t" $2 "\t" $3
    i++
    next
}
{
    # For each variant, check if breakpoint positions overlap any BED region
    for(k in a){
        # Split stored BED region into components
        split(a[k], x, "\t")
        # Check: chromosome match AND position within start/end coordinates
        if($2==x[1] && $3>=x[2] && $3<=x[3]) 
            print $0
    }
}' "$BED_FILE" tmp_all.tsv >> "$OUT"

# ==============================================================================
# SECTION 7: CLEAN UP TEMPORARY FILES
# ==============================================================================
# Remove intermediate files created during processing

# Delete temporary TSV file containing all extracted variants
rm -f tmp_all.tsv

echo "Done"
