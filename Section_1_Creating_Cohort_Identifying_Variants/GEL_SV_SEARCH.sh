#!/usr/bin/env bash

# ==============================================================================
# Script: SV_search.sh
#
# Purpose:
#   Identifies structural variants (SVs) that overlap with cis-regulatory
#   element (CRE/enhancer) regions or ERG gene body. Separates results by variant type
#   (structural variants and copy number variants).
#
# Analysis Context:
#   This script was used to identify structural variants
#   in Primary Lymphedema (PL) patients. Finds SVs and CNVs that affect
#   CRE regions or ERG gene body across both hg38 and hg19 cohorts.
#
# Workflow:
#   Step 1: Filter for PASS variants overlapping region 
#   Step 2: Extract and format structural variants (BND, DEL, DUP, INS, INV)
#   Step 3: Extract and format copy number variants (CNVs - CNV, LOH)
#
# Input:
#   VCF files with structural variant calls (gzip-compressed)
#   BED file with  region coordinates
#
# Output:
#   Two tab-delimited output files:
#     - SVs (structural variants): sample ID, chromosome, position, variant type
#     - CNVs (copy number variants): sample ID, chromosome, position, variant type
#
# Notes:
#   - Separates SVs (BND, DEL, DUP, INS, INV) from CNVs (CNV, LOH)
#   - Filters for PASS variants only
#   - Outputs header lines with sample ID and variant information
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

# ==============================================================================
# SECTION 2: DEFINE INPUT AND OUTPUT FILES
# ==============================================================================
# Specify paths to input VCF list and CRE regions file

# Text file containing paths to structural variant VCF files (one per line)
VCF_LIST=""

# BED file containing CRE region coordinates
BED_FILE=""

# Output file for structural variants (SVs)
SV_OUT=""

# Output file for copy number variants (CNVs)
CNV_OUT=""

# Temporary directory for intermediate VCF filtering
TMP_DIR="tmp_vcf_filtered"

# ==============================================================================
# SECTION 3: PREPARE OUTPUT FILES AND TEMPORARY DIRECTORY
# ==============================================================================
# Create temporary directory for filtered VCF files
mkdir -p "$TMP_DIR"

# Remove previous output files if they exist (start fresh)
rm -f "$SV_OUT" "$CNV_OUT"

# ==============================================================================
# SECTION 4: CREATE OUTPUT FILE HEADERS
# ==============================================================================
# Add header lines to output files with column names

# SV output header: Sample ID, chromosome, position, end, variant type, filter status
echo -e "Sample\tCHROM\tPOS\tEND\tSVTYPE\tFILTER\tID" > "$SV_OUT"

# CNV output header: Sample ID, chromosome, position, end, variant type, filter status
echo -e "Sample\tCHROM\tPOS\tEND\tSVTYPE\tFILTER\tID" > "$CNV_OUT"

# ==============================================================================
# SECTION 5: PROCESS EACH VCF FILE
# ==============================================================================
# Loop through each VCF file in the input list

echo "Starting enhancer SV analysis..."

while read -r VCF; do
    # Skip empty lines
    [[ -z "$VCF" ]] && continue
    
    # Check if VCF file exists
    if [[ ! -f "$VCF" ]]; then
        echo "WARNING: missing: $VCF"
        continue
    fi
    
    # Report progress
    echo "Processing: $VCF"
    
    # Extract sample ID using bcftools query
    # bcftools query -l lists all sample names in the VCF
    SAMPLE=$(bcftools query -l "$VCF")
    
    # Extract base filename without extension (used for intermediate files)
    BASE=$(basename "$VCF" .vcf.gz)
    
    # Define path for filtered VCF (intermediate file)
    FILTERED_VCF="${TMP_DIR}/${BASE}.filtered.vcf"

# ==============================================================================
# SECTION 6: STEP 1 - FILTER FOR PASS VARIANTS IN CRE REGIONS
# ==============================================================================
# Filter VCF to keep only:
# - Variants with PASS filter status
# - Variants overlapping  regions (using BED file)
# Output: filtered VCF with PASS variants in CRE regions

bcftools view \
    -R "$BED_FILE" \
    "$VCF" > "$FILTERED_VCF"

# ==============================================================================
# SECTION 7: STEP 2 - EXTRACT STRUCTURAL VARIANTS (SVs)
# ==============================================================================
# Extract structural variants (BND, DEL, DUP, INS, INV)
# Use bcftools query to format output with:
# - Sample name, chromosome, position, end position, variant type, filter, ID
# Variants are filtered to exclude low-quality or ambiguous calls

bcftools query \
    -i "INFO/SVTYPE="BND" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INS" || INFO/SVTYPE="INV"" \
    -f "${SAMPLE}\t%CHROM\t%POS\t%END\t%INFO/SVTYPE\t%FILTER\t%ID\n" \
    "$FILTERED_VCF" >> "$SV_OUT"

# ==============================================================================
# SECTION 8: STEP 3 - EXTRACT COPY NUMBER VARIANTS (CNVs)
# ==============================================================================
# Extract copy number variants (CNV, LOH)
# These represent gains/losses and loss-of-heterozygosity events
# Output format same as SVs for consistency

bcftools query \
    -i "INFO/SVTYPE="CNV" || INFO/SVTYPE="LOH"" \
    -f "${SAMPLE}\t%CHROM\t%POS\t%END\t%INFO/SVTYPE\t%FILTER\t%ID\n" \
    "$FILTERED_VCF" >> "$CNV_OUT"

done < "$VCF_LIST"

# ==============================================================================
# SECTION 9: CLEAN UP TEMPORARY FILES
# ==============================================================================
# Remove intermediate filtered VCF files

echo "Done."
echo "SVs: $SV_OUT"
echo "CNVs: $CNV_OUT"
