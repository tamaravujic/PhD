## used in Enhancer_specificty_heatmaps.ipynb to generate dbacount file

#PBS -l walltime=03:00:00
#PBS -l select=1:mem=120gb:ncpus=16

cd /rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/heatmap

module load anaconda3/personal
source activate diffbind

Rscript - <<EOF
library(DiffBind)
library(tidyverse)
#library(ComplexHeatmap)
library(ggplot2)
library(GenomicRanges)

load(file='dbObj')
load(file='peaks')
dbObj.counted <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,summits=100,peaks=peaks, bParallel=TRUE) #DBA_NORM_DEFAULT

save(dbObj.counted, file='dbObj_with_TV_CREs.counted')
EOF
