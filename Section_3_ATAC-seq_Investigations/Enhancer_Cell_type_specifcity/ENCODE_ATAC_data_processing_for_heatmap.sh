#!/bin/bash

# ==============================================================================
# Script: ENCODE_ATAC_processing.sh
#
# Purpose:
#   Processes ENCODE ATAC-seq data: combines replicates, intersects peaks with
#   ERG regulatory regions, and formats output for downstream analysis.
#
# Workflow:
#   1. Combine rep1 and rep2 broadPeak files for each cell type
#   2. Convert to BED format and add cell type annotation
#   3. Intersect peaks with ERG_CREs.bed reference file
#   4. Add cell type name to final BED output
#
# Input:
#   Directory structure: ENCODE_ATAC_COPY/
#     ├── celltype1/rep1_sorted_peaks.broadPeak
#     ├── celltype1/rep2_sorted_peaks.broadPeak
#     ├── celltype2/rep1_sorted_peaks.broadPeak
#     └── celltype2/rep2_sorted_peaks.broadPeak
#   Reference: ../ERG_CREs.bed (ERG regulatory elements)
#
# Output:
#   For each cell type folder:
#     - celltype_combined.broadPeak (merged replicates)
#     - celltype_named_combined.bed (converted to BED)
#     - celltype_enhancer_peaks.bed (intersected with ERG_CREs)
#     - celltype_enhancer_peaks_named.bed (final output with cell type column)
#
# Prerequisites:
#   - BEDTools module (bedtools intersect command)
#   - ERG_CREs.bed file in parent directory
# ==============================================================================

# Define directories
main_dir="ENCODE_ATAC_COPY"
erg_ref="../ERG_CREs.bed"

# Navigate to main directory
cd "$main_dir" || exit 1

echo "Processing ENCODE ATAC-seq data..."

# Iterate through each cell type folder
for folder in */; do
	folder_name=$(basename "$folder")
	echo "Processing $folder_name..."
	
	# Step 1: Combine rep1 and rep2 broadPeak files
	combined_broad="${folder}${folder_name}_combined.broadPeak"
	cat "${folder}rep1_sorted_peaks.broadPeak" "${folder}rep2_sorted_peaks.broadPeak" > "$combined_broad"
	
	# Step 2: Convert broadPeak to BED format and add cell type name
	# Extract cols 1-3 (chr, start, end) and col 10 (qvalue), add cell type name
	named_bed="${folder}${folder_name}_named_combined.bed"
	awk -F'\t' -v cell_type="$folder_name" '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" cell_type}' \
		"$combined_broad" > "$named_bed"
	
	# Step 3: Intersect with ERG_CREs reference file
	# -c: count overlaps, -wa: write original features from ERG_CREs
	intersected="${folder}${folder_name}_enhancer_peaks.bed"
	bedtools intersect -c -wa -a "$erg_ref" -b "$named_bed" > "$intersected"
	
	# Step 4: Format final output with cell type column
	# Rearrange columns for final BED output
	final_output="${folder}${folder_name}_enhancer_peaks_named.bed"
	awk -F'\t' -v cell_type="$folder_name" \
		'{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" cell_type "\t" $5}' \
		"$intersected" > "$final_output"
	
	echo "  Created: $final_output"
done

echo "ENCODE ATAC processing complete!"
