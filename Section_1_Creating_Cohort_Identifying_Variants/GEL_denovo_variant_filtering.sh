#!/bin/bash

# ==============================================================================
# Script: denovo_variant_filtering.sh
#
# Purpose:
#   Filters de novo variants from GEL data for chr21, intersects with ERG
#   regulatory regions (enhancers and gene body), and identifies affected
#   families in the primary lymphedema cohort.
#
# Workflow:
#   1. Convert TSV variant files to BED format
#   2. Intersect de novo variants with ERG enhancer regions
#   3. Intersect de novo variants with ERG gene body
#   4. Count unique variants and families in each dataset
#   5. Extract family IDs from variants in ERG regions
#   6. Combine family IDs across genome assemblies
#   7. Identify cohort patients with de novo variants in ERG regions
#
# Prerequisites:
#   - denovo_flagged_variants_grch38.tsv (downloaded from GEL)
#   - denovo_flagged_variants_grch37.tsv (downloaded from GEL)
#   - ERG_ENHANCERS_grch38.bed, ERG_ENHANCERS_grch37.bed
#   - ERG_BODY_grch38.bed, ERG_BODY_grch37.bed
#   - cohort.tsv (list of patients)
#
# Input Files:
#   - denovo_flagged_variants_grch38/37.tsv 
#   - ERG bed files 
#   - cohort.tsv
#
# Output Files:
#   - *_grch38.bed, *_grch37.bed (formatted variant files)
#   - ERG_ENHANCERS_DE_NOVO_INTERSECTED_grch38/37.txt
#   - ERG_BODY_DE_NOVO_INTERSECTED_grch38/37.txt
#   - family_IDs_ERG_enhancers.txt, family_IDs_ERG_body.txt
#   - patient_exp_ERG_enhancer_var.txt, patient_exp_ERG_body_var.txt
# ==============================================================================

# ==============================================================================
# SECTION 1: VARIABLES AND SETUP
# ==============================================================================

# Genome assemblies to process
ASSEMBLIES=("grch38" "grch37")

# Input files
INPUT_TSV_PREFIX="denovo_flagged_variants_"
ERG_ENHANCERS_PREFIX="ERG_ENHANCERS_"
ERG_BODY_PREFIX="ERG_BODY_"
COHORT_FILE="cohort.tsv"

# ==============================================================================
# SECTION 2: CONVERT TSV FILES TO BED FORMAT
# ==============================================================================
# Convert denovo_flagged_variants files from TSV to BED format
# BED format: chromosome, start-1, start, variant_ID
# Also sort and remove duplicates

echo "Converting TSV files to BED format..."
for assembly in "${ASSEMBLIES[@]}"; do
	input_file="${INPUT_TSV_PREFIX}${assembly}.tsv"
	output_file="${INPUT_TSV_PREFIX}${assembly}.bed"
	
	awk -v FS='\t' -v OFS='\t' '{print $4, $5-1, $5, $2}' "$input_file" | \
	sort -k 1V,1 -k2,2n | \
	uniq > "$output_file"
	
	echo "Created $output_file"
done

# ==============================================================================
# SECTION 3: LOAD BEDTOOLS MODULE
# ==============================================================================
# Check how modules are loaded into the terminal and update line
module load bio/BEDTools/2.27.1-foss-2018b

# ==============================================================================
# SECTION 4: INTERSECT DE NOVO VARIANTS WITH ERG REGIONS
# ==============================================================================
# Intersect de novo variants with ERG enhancers and gene body
# -wo: write overlapping features and overlap length
# -u: write only features from first file with overlap

echo "Intersecting variants with ERG regions..."
for assembly in "${ASSEMBLIES[@]}"; do
	variant_bed="${INPUT_TSV_PREFIX}${assembly}.bed"
	enhancer_bed="${ERG_ENHANCERS_PREFIX}${assembly}.bed"
	erg_body_bed="${ERG_BODY_PREFIX}${assembly}.bed"
	
	# Variants overlapping ERG enhancers
	bedtools intersect -wo -a "$enhancer_bed" -b "$variant_bed" \
		> "ERG_ENHANCERS_DE_NOVO_INTERSECTED_${assembly}.txt"
	
	# Variants overlapping ERG gene body
	bedtools intersect -u -a "$variant_bed" -b "$erg_body_bed" \
		> "ERG_BODY_DE_NOVO_INTERSECTED_${assembly}.txt"
	
	echo "Created intersection files for $assembly"
done

# ==============================================================================
# SECTION 5: COUNT UNIQUE VARIANTS AND FAMILIES
# ==============================================================================
# Summary statistics: how many unique variants and families in each dataset
# sort -u removes duplicates from specified column

echo "Counting unique variants and families..."

for assembly in "${ASSEMBLIES[@]}"; do
	input_tsv="${INPUT_TSV_PREFIX}${assembly}.tsv"
	enhancer_intersect="ERG_ENHANCERS_DE_NOVO_INTERSECTED_${assembly}.txt"
	body_intersect="ERG_BODY_DE_NOVO_INTERSECTED_${assembly}.txt"
	
	echo "=== $assembly ==="
	echo "Total variants:"
	sort -u -k5,5 "$input_tsv" | wc -l
	echo "Total families:"
	sort -u -k5,5 "$input_tsv" | wc -l
	
	echo "Enhancers - unique variants:"
	sort -u -k5,5 "$enhancer_intersect" | wc -l
	echo "Enhancers - unique families:"
	sort -u -k7,7 "$enhancer_intersect" | wc -l
	
	echo "Gene body - unique variants:"
	sort -u -k3,3 "$body_intersect" | wc -l
	echo "Gene body - unique families:"
	sort -u -k4,4 "$body_intersect" | wc -l
	echo ""
done

# ==============================================================================
# SECTION 6: EXTRACT FAMILY IDS FROM VARIANT-REGION INTERSECTIONS
# ==============================================================================
# Get family/patient IDs from variants in ERG regions

echo "Extracting family IDs..."
for assembly in "${ASSEMBLIES[@]}"; do
	enhancer_intersect="ERG_ENHANCERS_DE_NOVO_INTERSECTED_${assembly}.txt"
	body_intersect="ERG_BODY_DE_NOVO_INTERSECTED_${assembly}.txt"
	
	# Column 4 for body intersects, column 7 for enhancer intersects
	cut -f4 "$body_intersect" | sort | uniq > "family_IDs_ERG_body_${assembly}.txt"
	cut -f7 "$enhancer_intersect" | sort | uniq > "family_IDs_ERG_enhancers_${assembly}.txt"
done

# ==============================================================================
# SECTION 7: COMBINE FAMILY IDS ACROSS GENOME ASSEMBLIES
# ==============================================================================
# Merge family IDs from both GRCh37 and GRCh38 into single files
# (some families may have variants in both assemblies)

echo "Combining family IDs across genome assemblies..."

# Combine ERG body family IDs
> family_IDs_ERG_body.txt
for assembly in "${ASSEMBLIES[@]}"; do
	cat "family_IDs_ERG_body_${assembly}.txt" >> family_IDs_ERG_body.txt
done
sort -u family_IDs_ERG_body.txt > temp && mv temp family_IDs_ERG_body.txt

# Combine ERG enhancer family IDs
> family_IDs_ERG_enhancers.txt
for assembly in "${ASSEMBLIES[@]}"; do
	cat "family_IDs_ERG_enhancers_${assembly}.txt" >> family_IDs_ERG_enhancers.txt
done
sort -u family_IDs_ERG_enhancers.txt > temp && mv temp family_IDs_ERG_enhancers.txt

echo "Created combined family ID files"

# ==============================================================================
# SECTION 8: IDENTIFY COHORT PATIENTS WITH DE NOVO VARIANTS
# ==============================================================================
# Cross-reference family IDs with cohort file to identify affected patients

echo "Identifying patients in cohort with de novo variants..."

grep -wf family_IDs_ERG_enhancers.txt "$COHORT_FILE" > patient_exp_ERG_enhancer_var.txt
grep -wf family_IDs_ERG_body.txt "$COHORT_FILE" > patient_exp_ERG_body_var.txt

echo "Variant filtering complete!"
echo "Output files:"
echo "  - patient_exp_ERG_enhancer_var.txt"
echo "  - patient_exp_ERG_body_var.txt"
