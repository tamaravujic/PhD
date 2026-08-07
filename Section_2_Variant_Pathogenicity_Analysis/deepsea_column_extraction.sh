#!/bin/bash

# ==============================================================================
# Script: deepsea_column_extraction.sh
#
# Purpose:
#   Extract and reformat specific columns from DeepSEA output files.
#   Useful for selecting columns of interest from large DeepSEA predictions
#   and reorganizing them into a new TSV file.
#
# Workflow:
#   1. Format header text file (convert space-separated to line-separated)
#   2. Read header list and identify column indices in DeepSEA file
#   3. Extract selected columns from DeepSEA output
#   4. Save extracted columns to new TSV file
#
# Input:
#   my_headers.tsv - space-separated or tab-separated header names to extract
#   *_DEEPSEA_diffs.tsv - DeepSEA difference scores file
#
# Output:
#   extracted_headers.txt - line-separated list of headers
#   extracted_columns.tsv - TSV file with only selected columns
#
# Notes:
#   - Preserves header row and data rows
#   - Uses awk for efficient column selection
#   - Handles tab-delimited input
# ==============================================================================

# ==============================================================================
# SECTION 1: FORMAT HEADER TEXT FILE
# ==============================================================================
# Convert space-separated or tab-separated headers into one header per line
# This creates a list file that we'll use to identify which columns to extract

echo "Formatting header file..."
tr ' ' '\n' < my_headers.tsv > extracted_headers.txt

echo "Created extracted_headers.txt"

# ==============================================================================
# SECTION 2: EXTRACT SELECTED COLUMNS FROM DEEPSEA FILE
# ==============================================================================
# Use awk to:
# 1. Read header row and store column indices
# 2. Match headers from extracted_headers.txt with column indices
# 3. Extract and print only matching columns for all rows

echo "Extracting selected columns from DeepSEA file..."

awk -F'\t' '
NR==1 {
    # For each column in the header row, store its index keyed by header name
    for (i=1; i<=NF; i++) {
        header[$i] = i
    }
    
    # Read the list of headers to extract
    while ((getline line < "extracted_headers.txt") > 0) {
        if (line in header) {
            cols[header[line]] = 1  # Mark this column index for extraction
        }
    }
    
    # Print selected headers
    first = 1
    for (i=1; i<=NF; i++) {
        if (i in cols) {
            printf "%s%s", (first ? "" : "\t"), $i
            first=0
        }
    }
    printf "\n"
}
NR>1 {
    # Print values for the selected columns from data rows
    first = 1
    for (i=1; i<=NF; i++) {
        if (i in cols) {
            printf "%s%s", (first ? "" : "\t"), $i
            first=0
        }
    }
    printf "\n"
}
' 7c6a0fdb-fe04-43a1-9bc6-c558c45382d5_DEEPSEA_diffs.tsv > extracted_columns.tsv

echo "Extraction complete!"
echo "Output saved to: extracted_columns.tsv"
