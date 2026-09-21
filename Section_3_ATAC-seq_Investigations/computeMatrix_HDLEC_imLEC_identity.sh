#PBS -l select=1:mem=60gb:ncpus=16
#PBS -l walltime=02:00:00
#PBS -N computeMatrix_HDLEC_imLEC_identity

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/tv722/home/anaconda3/envs/ATAC

DIR=/rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/ECs_ATAC/imlec_atac_comparison
REFDIR=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38

LEC1=$DIR/bigwigs/LEC_rep1.BPM.bw
LEC2=$DIR/bigwigs/LEC_rep2.BPM.bw
IML1=$DIR/bigwigs/imLEC_rep1.BPM.bw
IML2=$DIR/bigwigs/imLEC_rep2.BPM.bw

# Genome-wide HDLEC vs imLEC identity: regions specific to each cell type
# plus a subsampled shared set, using only the two cell types of interest
# (HUVEC dropped -- not relevant to the CRISPRi model-validity question)
computeMatrix reference-point --referencePoint center -p 16 -S $LEC1 $LEC2 $IML1 $IML2 -R $DIR/shared_HDLEC_imLEC_25000.bed $DIR/topLEC_vsIMLEC.bed $DIR/topIMLEC_vsLEC.bed -bs 100 -b 2000 -a 2000 -o $DIR/HDLEC_imLEC_identity_shared_first.matrix.gz --blackListFileName $REFDIR/hg38-blacklist.v2.bed