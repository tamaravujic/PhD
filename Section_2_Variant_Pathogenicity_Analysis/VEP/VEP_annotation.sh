#!/bin/bash

# ==============================================================================
# Script: VEP_annotation.sh
#
# Purpose:
#   Variant Effect Predictor (VEP) annotation pipeline. Annotates VCF file
#   with functional and population frequency data using multiple databases.
#   Runs on cluster via PBS job submission.
#
# Workflow:
#   1. Run VEP on input VCF with annotation plugins (REVEL, CADD, gnomAD, SpliceAI)
#   2. Compress and index output VCF
#   3. Convert annotated VCF to TSV format for easier analysis
#   4. (Optional) Filter and process TSV output
#
# Input:
#   variants.vcf.gz - compressed VCF file with variants to annotate
#
# Output:
#   variants_annotated_variants.vcf.gz - annotated VCF (bgzip compressed)
#   variants_annotated_variants.tsv - tab-separated annotation table
#
# Annotation Databases:
#   - REVEL: pathogenicity predictions for missense variants
#   - CADD: variant deleteriousness scores
#   - gnomAD: population allele frequencies (multiple populations)
#   - SpliceAI: splice site disruption predictions
#
# Notes:
#   - Uses Singularity container for VEP installation
#   - Assembly: GRCh38
# ==============================================================================

#PBS -l select=1:mem=60gb:ncpus=20
#PBS -l walltime=8:00:00
#PBS -N VEP_hg38_3

# ==============================================================================
# SECTION 1: DEFINE VARIABLES AND PATHS
# ==============================================================================

# Input VCF file
INROAD=/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/variants.vcf.gz
SAMPLE="variants"

# Directory structure
INPUT_DIR="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/"
OUTPUT_DIR="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/annotated"
VEP_DATA="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/vep_data/"
REF="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/reference_genomes"

# Annotation databases
REVEL="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/REVEL"
CADD="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/CADD_v1.7"
GNOMAD="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/gnomAD_data/"
SpliceAI="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/spliceAI"

# VEP container
VEP_CONTAINER="/rds/general/user/tv722/projects/tamara_vujic_phd/live/ensembl-vep/vep.sif"

# ==============================================================================
# SECTION 2: RUN VEP ANNOTATION (VCF OUTPUT)
# ==============================================================================
# First VEP run: annotate variants and output as compressed VCF

singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR},${REF},${VEP_DATA},${REVEL},${CADD},${GNOMAD},${SpliceAI} \
  ${VEP_CONTAINER} \
  vep --dir ${VEP_DATA} \
  --cache --offline --format vcf --vcf --force_overwrite \
  --fork 20 \
  --input_file ${INROAD} \
  --output_file ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.vcf.gz \
  --compress_output bgzip \
  --assembly GRCh38 \
  --plugin REVEL,${REVEL}/new_tabbed_revel_grch38.tsv.gz \
  --plugin CADD,${CADD}/whole_genome_SNVs.tsv.gz,${CADD}/gnomad.genomes.r4.0.indel.tsv.gz \
  --plugin gnomADc,${GNOMAD}/gnomad.ch.genomesv3.tabbed.tsv.gz \
  --plugin SpliceAI,snv=${SpliceAI}/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${SpliceAI}/spliceai_scores.raw.indel.hg38.vcf.gz \
  --custom file=${GNOMAD}/gnomad.genomes.r2.1.sites.grch38.chr9_noVEP.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF%AF_AFR%AF_AMR%AF_ASJ%AF_EAS%AF_FIN%AF_NFE%AF_OTH%AF_ami

# ==============================================================================
# SECTION 3: CONVERT ANNOTATED VCF TO TSV FORMAT
# ==============================================================================
# Second VEP run: convert VCF to tab-separated format for easier analysis
# Output includes all annotation fields

singularity exec \
  -B ${OUTPUT_DIR},${VEP_DATA},${REF},${REVEL},${CADD},${GNOMAD},${SpliceAI} \
  ${VEP_CONTAINER} \
  vep --input_file ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.vcf.gz \
      --output_file ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.tsv \
      --tab \
      --force_overwrite \
      --cache --offline \
      --dir ${VEP_DATA} \
      --assembly GRCh38 \
      --everything \
      --plugin REVEL,${REVEL}/new_tabbed_revel_grch38.tsv.gz \
      --plugin CADD,${CADD}/whole_genome_SNVs.tsv.gz,${CADD}/gnomad.genomes.r4.0.indel.tsv.gz \
      --plugin gnomADc,${GNOMAD}/gnomad.ch.genomesv3.tabbed.tsv.gz \
      --plugin SpliceAI,snv=${SpliceAI}/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${SpliceAI}/spliceai_scores.raw.indel.hg38.vcf.gz

# ==============================================================================
# SECTION 4: OPTIONAL TERMINAL FILTERING COMMANDS (for manual post-processing)
# ==============================================================================
# These commands can be run in the terminal after VEP completes for filtering
# TSV output. Uncomment and modify as needed.

# Select only main transcript or canonical versions
# awk -F"\t" 'NR==1 || $24 != "-" || $23 == "YES"' ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.tsv > ${OUTPUT_DIR}/${SAMPLE}_main_transcripts.tsv

# Remove duplicates that do not have a gene SYMBOL
# awk -F"\t" '!seen[$1]++ || NR==1' ${OUTPUT_DIR}/${SAMPLE}_main_transcripts.tsv > ${OUTPUT_DIR}/${SAMPLE}_main_transcripts_dup_removed.tsv

# Collapse all rows per variant into single row (preserves all information)
# awk -F"\t" '
# NR==1 {
#   header=$0
#   ncol=NF
#   next
# }
# {
#   key=$1
#   seen[key]=1
#   for (i=1; i<=NF; i++) {
#     if (vals[key,i] == "") {
#       vals[key,i]=$i
#     } else if (vals[key,i] !~ "(^|;)"$i"(;|$)") {
#       vals[key,i]=vals[key,i]";"$i
#     }
#   }
# }
# END {
#   print header
#   for (k in seen) {
#     row=""
#     for (i=1; i<=ncol; i++) {
#       row = row (i==1 ? vals[k,i] : OFS vals[k,i])
#     }
#     print row
#   }
# }' OFS="\t" ${OUTPUT_DIR}/${SAMPLE}_main_transcripts.tsv > ${OUTPUT_DIR}/${SAMPLE}_main_transcripts_collapsed.tsv

# ==============================================================================
# SECTION 5: OPTIONAL - EXTRACT SPECIFIC FIELDS FROM ANNOTATED VCF
# ==============================================================================
# Use bcftools to query and extract specific annotation fields
# Uncomment and modify as needed

# Extract specific INFO fields into TSV
# bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/AF\t%INFO/CADD\t%INFO/REVEL\t%INFO/SpliceAI_DS_AG\n' \
#   ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.vcf.gz > ${OUTPUT_DIR}/${SAMPLE}_selected_fields.tsv

# View all available INFO fields in header
# bcftools view -h ${OUTPUT_DIR}/${SAMPLE}_annotated_variants.vcf.gz | grep "^##INFO"
