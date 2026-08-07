#!/bin/bash

# ==============================================================================
# Script: Submit_dbaCount.sh
#
# Purpose:
#   PBS cluster submission script for DiffBind peak counting step.
#   Reads pre-loaded DiffBind objects and performs read count normalisation
#   across ATAC-seq samples. Used in Enhancer_specificity_heatmaps.ipynb workflow.
#
# Workflow:
#   1. Load R/DiffBind environment on cluster
#   2. Read previously saved DiffBind sample object
#   3. Count reads in peaks using summarizeOverlaps
#   4. Save normalised count object for downstream analysis
#
# Prerequisites:
#   - dbObj: saved DiffBind sample sheet object (from Enhancer_specificity_heatmaps.ipynb)
#   - peaks: saved GRanges object of merged peaks (from notebook)
#   - diffbind conda environment configured on cluster
#
# Output:
#   dbObj_with_TV_CREs.counted: DiffBind object with normalised peak counts
#
# Usage:
#   qsub Submit_dbaCount.sh
# ==============================================================================

# Navigate to working directory containing dbObj and peaks files
cd /rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/heatmap

# Load conda environment 
module load anaconda3/personal
source activate diffbind

# Run R script inline 
Rscript - <<EOF

# Load required libraries
library(DiffBind)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)

# Load pre-computed DiffBind sample object and peak regions
load(file='dbObj')
load(file='peaks')

# Count reads in peaks using DiffBind's dba.count function
# Parameters:
#   bUseSummarizeOverlaps: use GenomicAlignments::summarizeOverlaps for counting
#   summits: peak region extension (±100 bp from summit)
#   peaks: use merged peak set defined above
#   bParallel: use 16 CPUs for parallel processing (speeds up counting)
dbObj.counted <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits=100, 
                           peaks=peaks, bParallel=TRUE)

# Save normalised count object for downstream analysis (visualization, statistics)
save(dbObj.counted, file='dbObj_with_TV_CREs.counted')

EOF
