#PBS -l select=1:mem=40gb:ncpus=16
#PBS -l walltime=01:00:00
#PBS -N computeMatrix_ERG_functional_HDLEC_imLEC

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/tv722/home/anaconda3/envs/ATAC

DIR=/rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/ECs_ATAC/imlec_atac_comparison
REFDIR=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38


LEC1=$DIR/bigwigs/LEC_rep1.BPM.bw
LEC2=$DIR/bigwigs/LEC_rep2.BPM.bw
IML1=$DIR/bigwigs/imLEC_rep1.BPM.bw
IML2=$DIR/bigwigs/imLEC_rep2.BPM.bw

computeMatrix reference-point --referencePoint center -p 16 \
    -S $LEC1 $LEC2 $IML1 $IML2 \
    -R $DIR/ERG_CREs_functional_retained_in_imLEC.bed $DIR/ERG_CREs_functional_lost_in_imLEC.bed \
    -bs 50 -b 1000 -a 1000 \
    -o $DIR/ERG_functional_HDLEC_imLEC.matrix.gz \
    --blackListFileName $REFDIR/hg38-blacklist.v2.bed \
    --sortRegions keep