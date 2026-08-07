#!/bin/bash
# ==============================================================================
# Script: variant_pipeline.sh
#
# Purpose:
#   Complete variant analysis pipeline for filtering and consolidating
#   genomic variants. Extracts variants from desired genomic regions then
#   filters for rare variants, and consolidates results into a single
#   tab-delimited output file.
#
#
#   Analysis 1 - Non-coding variants in enhancer regions:
#     - Extracted variants within cis-regulatory elements (enhancers)
#     - Filtered for rare variants (AF < 1%)
#     - Run on both hg38 and hg19 cohorts
#
#   Analysis 2 - Coding variants in PL genes:
#     - Identified variants in 40 known PL genes
#     - Filtered for rare  variants
#     - Run on both hg38 and hg19 genome assemblies
#
#   Each analysis used configuration files (variant_pipeline.conf) with
#   appropriate paths to region files or gene coordinates for each assembly.
#
# Workflow:
#   Step 1: Extract variants in  regions (PASS filter only)
#   Step 2: Filter for rare variants (allele frequency < 0.01)
#   Step 3: Consolidate results into single TSV file with variant metadata
#
# Input:
#   Configuration file (passed as command-line argument)
#   VCF paths in txt file (paths specified in configuration)
#   regions txt file (path specified in configuration)
#
# Output:
#   Directory structure with three subdirectories:
#     01_extracted/ - PASS variants overlapping CRE regions (.vcf.gz)
#     02_filtered/ - Rare variants after AF filtering (.vcf.gz)
#     03_consolidated.tsv - Final results table with all variant data
#
# Usage:
#   in terminal
# ./variant_pipeline.sh variant_pipeline.conf
#
# ==============================================================================
source "$1"

mkdir -p "$WORK_DIR/01_extracted" "$WORK_DIR/02_filtered" 

echo "Step 1: Extracting variants in regions (PASS only)..."
step1_total=0
while read vcf; do
    [[ -z "$vcf" ]] && continue
    patient=$(basename "$vcf" .vcf.gz)
    bedtools intersect -a "$vcf" -b "$CRE_REGIONS" -header | awk '$7 == "PASS" || /^#/' | bgzip -c > "$WORK_DIR/01_extracted/${patient}_extracted.vcf.gz"
    tabix -p vcf "$WORK_DIR/01_extracted/${patient}_extracted.vcf.gz"
    count=$(zcat "$WORK_DIR/01_extracted/${patient}_extracted.vcf.gz" | grep -v '^#' | wc -l)
    echo "  $patient: $count variants"
    step1_total=$((step1_total + count))
done < "$VCF_PATHS_FILE"
echo "Step 1 total: $step1_total variants"
echo ""

echo "Step 2: Filtering for rare variants (AF < $AF_THRESHOLD or no AF)..."
step2_total=0
for vcf in "$WORK_DIR/01_extracted"/*.vcf.gz; do
    patient=$(basename "$vcf" _extracted.vcf.gz)
    bcftools view "$vcf" | awk -v threshold="$AF_THRESHOLD" 'BEGIN {FS="\t"; OFS="\t"} /^#/ {print; next} {
        af="NA"
        split($8, info, ";")
        for(i in info) if(info[i] ~ /^AF1000G=/) {af=substr(info[i],9); break}
        if(af=="NA" || af=="") print
        else if(af+0 < threshold+0) print
    }' | bgzip -c > "$WORK_DIR/02_filtered/${patient}_filtered.vcf.gz"
    tabix -p vcf "$WORK_DIR/02_filtered/${patient}_filtered.vcf.gz"
    count=$(zcat "$WORK_DIR/02_filtered/${patient}_filtered.vcf.gz" | grep -v '^#' | wc -l)
    echo "  $patient: $count variants"
    step2_total=$((step2_total + count))
done
echo "Step 2 total: $step2_total variants"
echo ""

echo "Step 3: Consolidating to single file..."
step3_total=0
echo -e "patient_id\tchromosome\tposition\trsid\tref\talt\tfilter\tgenotype\tallele_frequency" > "$WORK_DIR/03_consolidated.tsv"
for vcf in "$WORK_DIR/02_filtered"/*.vcf.gz; do
    patient=$(basename "$vcf" _filtered.vcf.gz)
    bcftools view "$vcf" | awk -v patient="$patient" 'BEGIN {FS="\t"; OFS="\t"} /^#/ {next} {
        af="NA"
        split($8, info, ";")
        for(i in info) if(info[i] ~ /^AF1000G=/) {af=substr(info[i],9); break}
        split($10, sample, ":")
        print patient, $1, $2, $3, $4, $5, $7, sample[1], af
    }' >> "$WORK_DIR/03_consolidated.tsv"
    count=$(zcat "$vcf" | grep -v '^#' | wc -l)
    echo "  $patient: $count variants"
    step3_total=$((step3_total + count))
done
echo "Step 3 total: $step3_total variants"
echo ""