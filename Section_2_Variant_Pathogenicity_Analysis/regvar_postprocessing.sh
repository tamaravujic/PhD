#!/bin/bash

# ==============================================================================
# Script: regvar_postprocessing.sh
#
# Purpose:
#   Post-processing for RegVar (regulatory variant prediction) output.
#   Cleans CSV formatting, filters by variant-gene list, and removes duplicates.
#
# Workflow:
#   1. Clean RegVar CSV file (remove quotes and carriage returns)
#   2. Filter results to keep only variants in whitelist
#   3. Remove duplicate lines from filtered output
#
# Input:
#   regvar-result-*.csv - RegVar output file (may have formatting issues)
#   correct_var_gene_list.csv - list of variant-gene pairs to keep
#     Format: col1=variant, col2=gene (comma-separated)
#
# Output:
#   regvar-result-*.csv.bak - backup of original file
#   filtered.csv - filtered variants (only those in whitelist)
#   filtered_nodup.csv - final output with duplicates removed
#
# Notes:
#   - Creates backup of original file before modifying
#   - Removes Windows-style carriage returns (\r)
#   - Filters to only variant-gene combinations in list
#   - When running RegVAR you need to input the genes where the variants are located.
#   This is a problem if the list contains variants from multiple genes as each variant
#   Will have a prediction for each gene. This script will keep only the output where
#   The variant and gene match. Then columns can be merged per tissue.
# ==============================================================================

# ==============================================================================
# SECTION 1: DEFINE INPUT FILES
# ==============================================================================
# Modify these paths to match your file locations

REGVAR_FILE="regvar-result-a7c55a738d3f4995b4771143e1b11769.csv"
WHITELIST="correct_var_gene_list.csv"
FILTERED_OUTPUT="filtered.csv"
FINAL_OUTPUT="filtered_nodup.csv"

# ==============================================================================
# SECTION 2: CLEAN REGVAR CSV FILE
# ==============================================================================
# Remove problematic characters from RegVar output:
# - Remove double quotes (sometimes present in CSV exports)
# - Remove carriage returns (Windows line endings \r)
# Creates backup of original file automatically

echo "Cleaning RegVar output file..."
sed -i.bak 's/"//g; s/\r//g' "$REGVAR_FILE"
echo "Backup created: ${REGVAR_FILE}.bak"

# ==============================================================================
# SECTION 3: FILTER TO WHITELIST VARIANTS
# ==============================================================================
# Keep only rows where (column1 + column2) combination exists in whitelist
# This allows you to focus on specific variant-gene pairs of interest

echo "Filtering to whitelist variants..."
awk -F',' '
NR==FNR { 
    # Read whitelist file first, store variant-gene combinations
    keep[$1 FS $2]; 
    next 
} 
# For each row in RegVar results, check if variant-gene pair is in whitelist
($1 FS $2) in keep
' "$WHITELIST" "$REGVAR_FILE" > "$FILTERED_OUTPUT"

echo "Filtered variants saved to: $FILTERED_OUTPUT"

# ==============================================================================
# SECTION 4: REMOVE DUPLICATE LINES
# ==============================================================================
# Keep only first occurrence of each unique line

echo "Removing duplicates..."
awk '!seen[$0]++' "$FILTERED_OUTPUT" > "$FINAL_OUTPUT"

echo "Final output saved to: $FINAL_OUTPUT"
echo "Post-processing complete!"
