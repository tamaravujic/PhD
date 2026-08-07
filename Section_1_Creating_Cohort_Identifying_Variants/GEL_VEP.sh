#!/bin/bash

# ==============================================================================
# Script: VEP_INSIDE_GEL.sh
#
# Purpose:
#   Variant Effect Predictor (VEP) annotation pipeline for variants within
#   GEL data. Annotates VCF files with functional predictions, filters by
#   allele frequency and pathogenicity scores, and retains only heterozygous
#   variants that pass quality filters.
#
# Workflow:
#   1. Run VEP 111 annotation with multiple plugins (dbNSFP, LoF, SpliceAI, CADD, REVEL, MutFunc)
#   2. Filter variants by allele frequency (AF < 0.01) and pathogenicity (CADD > 5 or REVEL > 0.6)
#   3. Further filter to keep only PASS variants that are heterozygous
#   4. Remove empty files and clean up intermediate outputs
#
# Input:
#   VCF file (gzip-compressed) with variants to annotate
#   Configuration file: vep_109.conf (defines environment variables and paths)
#
# Output:
#   *_annotated.vcf - VEP-annotated variants
#   *_filtered_variants.vcf - variants passing AF and pathogenicity filters
#   *_pass_filter_heterozygous.vcf - final output: heterozygous PASS variants
#
# Cluster Requirements (LSF #BSUB parameters):
#   - Queue: short
#   - Memory: 2000 MB
#   - CPUs: 1
#   - Project: re_gecip_cardiovascular
#
# Prerequisites:
#   - Singularity container with VEP 111 installed
#   - vep_109.conf configuration file
#   - Input VCF files
#   - Reference data (FASTA, cache, annotation databases)
#
# Notes:
#   - Runs inside GEL environment (re_gecip project space)
#   - Uses singularity containers for VEP execution
#   - Filters for heterozygous variants specifically
# ==============================================================================

#BSUB -q short
#BSUB -P re_gecip_cardiovascular
#BSUB -cwd "/re_gecip/cardiovascular/tvujic/109"
#BSUB -o vep_job_.stdout
#BSUB -e vep_job_.stderr
#BSUB -n 1
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000

# ==============================================================================
# SECTION 1: LOAD MODULES AND CONFIGURATION
# ==============================================================================
# Clear existing modules and load singularity for container-based VEP execution

module purge
module load singularity/3.8.3

# Source configuration file containing environment variables and paths
# This file defines: MOUNT_WD, MOUNT_GENOMES, MOUNT_GEL_DATA_RESOURCES, IMG, etc.
source /re_gecip/cardiovascular/tvujic/109/vep_109.conf

# Set Perl library path to include VEP plugin directories
export PERLLIB=$PERLLIB:/plugins/loftee-GRCh38:/plugins/loftee-GRCh37:/plugins

# ==============================================================================
# SECTION 2: DEFINE INPUT AND OUTPUT FILES
# ==============================================================================
# Specify input VCF file and generate output filename from input basename

input='/re_gecip/cardiovascular/tvujic/109/PL_patients_vcf_intersected_with_enhancers/LP3001187-DNA_G11.vcf.gz'

# Extract filename without extension (e.g., LP3001187-DNA_G11)
output='$(basename "$input" | cut -f 1,2 -d '.')''
echo "$output"

# ==============================================================================
# SECTION 3: RUN VEP 111 ANNOTATION WITH MULTIPLE PLUGINS
# ==============================================================================
# Execute VEP using singularity container with comprehensive annotation plugins
# Plugins include: dbNSFP, LoF, SpliceAI, CADD, REVEL, MutFunc

singularity run --bind "${MOUNT_WD}","${MOUNT_GENOMES}","${MOUNT_GEL_DATA_RESOURCES}","${MOUNT_PUBLIC_DATA_RESOURCES}","${MOUNT_SCRATCH}" "${IMG}" vep \
-cache \
--species homo_sapiens \
--format vcf \
--dir_plugins /plugins/ \
--offline --vcf \
--cache \
--force_overwrite \
--assembly GRCh38 \
--cache_version 109 \
--dir_cache "${CACHE}" \
--fasta "${REFFASTA}" \
--input_file "${input}" \
--warning_file "${output}"_errors.vcf \
--output_file "${output}"_annotated.vcf \
--plugin dbNSFP,"${DB_NSFP}","${DB_NSFP_REPLACEMENT_LOGIC}",ALL \
--plugin LoF,loftee_path:/plugins/loftee-GRCh38,human_ancestor_fa:"${LOFTEE38HA}",gerp_bigwig:"${LOFTEE38GERP}",conservation_file:"${LOFTEE38SQL}" \
--plugin SpliceAI,snv="${SPLICEAIRAW38}",indel="${SPLICEAIINDEL38}" \
--plugin CADD,"${CADD16}" \
--plugin REVEL,"/public_data_resources/vep_resources/REVEL/revel_v1.3_GRCh38.tsv.gz" \
--plugin mutfunc,db="${MUTFUNC_DB}"

# ==============================================================================
# SECTION 4: FILTER ANNOTATED VARIANTS BY ALLELE FREQUENCY AND PATHOGENICITY
# ==============================================================================
# Apply two-part filter:
# Part 1: Keep variants with allele frequency < 0.01 OR variants without AF information
# Part 2: Keep variants with CADD_PHRED score > 5 OR REVEL score > 0.6

singularity exec --bind "${MOUNT_WD}","${MOUNT_GENOMES}","${MOUNT_GEL_DATA_RESOURCES}","${MOUNT_PUBLIC_DATA_RESOURCES}","${MOUNT_SCRATCH}" "${IMG}" filter_vep -i "${output}"_annotated.vcf --filter "(AF < 0.01 or not AF)" --filter "(CADD_PHRED > 5 or REVEL > 0.6)" > "${output}"_filtered_variants.vcf

# Variants were then filtered so that only those that were heterozygous and PASS the filter were retained

# ==============================================================================
# SECTION 5: LOAD TOOLS FOR VCF PROCESSING
# ==============================================================================
# Load tools for VCF manipulation and filtering

module load vcftools/0.1.16
module load bcftools/1.16
module load tabix/1.18

# ==============================================================================
# SECTION 6: PROCESS FILTERED VARIANTS - HETEROZYGOUS FILTERING
# ==============================================================================
# For each filtered VCF file:
# 1. Extract basename for output file naming
# 2. Filter to keep only PASS variants and recode to VCF format
# 3. Compress and index with tabix
# 4. Extract heterozygous variants only
# 5. Clean up intermediate files

for file in "${output}"_filtered_variants.vcf; do

    # Extract the filename without the extension
    filename=$(basename "$file" .vcf)

    # Filter for variants that pass filter and are heterozygous
    # --recode: output as VCF format
    # --recode-INFO-all: preserve all INFO fields
    # --keep-filtered PASS: keep only variants with PASS filter
    vcftools --vcf "$file" --out "$filename" --recode --recode-INFO-all --keep-filtered PASS
    
    # Compress the recoded VCF file
    bgzip "$filename".recode.vcf
    
    # Index the compressed VCF file for quick access
    tabix -p vcf "$filename".recode.vcf.gz
    
    # Extract only heterozygous genotypes (-g het)
    bcftools view -g het "$filename".recode.vcf.gz -o "$filename"_pass_filter_heterozygous.vcf
    
    # Remove temporary compressed files
    rm -f *recode.vcf.gz *.recode.vcf.gz.tbi

# ==============================================================================
# SECTION 7: REMOVE EMPTY FILES
# ==============================================================================
# Clean up intermediate annotation summary and error files

rm -f "${output}"_annotated.vcf_summary.html -f "${output}"_errors.vcf

# If the filtered variants file contains no variants (only header), remove it
if ! grep -q -v "^#" "${output}"_filtered_variants.vcf; then
    rm -f "${output}"_filtered_variants.vcf
fi

done

# ==============================================================================
# SECTION 8: CLEAN UP EMPTY HETEROZYGOUS FILTER OUTPUT
# ==============================================================================
# If the final heterozygous-filtered file is empty (contains only header lines),
# remove it to keep output directory clean

if ! grep -q -v "^#" "$filename"_filtered_variants_pass_filter_heterozygous.vcf; then
    rm -f "$filename"_filtered_variants_pass_filter_heterozygous.vcf
fi

done
