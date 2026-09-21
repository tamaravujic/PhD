#PBS -l select=1:mem=64gb:ncpus=16
#PBS -l walltime=02:00:00
#PBS -N corr_pooled_EC

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/tv722/home/anaconda3/envs/ATAC

DIR=/rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/ECs_ATAC/imlec_atac_comparison
REFDIR=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38

cd $DIR

multiBigwigSummary bins -b $DIR/bigwigs/LEC_pooled.BPM.bw $DIR/bigwigs/imlec_ATAC_pooled.BPM.bw --labels HDLEC imLEC  --binSize 10000 --blackListFileName $REFDIR/hg38-blacklist.v2.bed -p 16 -out $DIR/EC_pooled_corr.npz

plotCorrelation -in $DIR/EC_pooled_corr.npz --corMethod pearson --whatToPlot heatmap --plotNumbers --skipZeros --removeOutliers --plotTitle "Pearson correlation, 10kb bins" -o $DIR/EC_pooled_heatmap.pdf --outFileCorMatrix $DIR/EC_pooled_Pearson.tsv

plotCorrelation -in $DIR/EC_pooled_corr.npz --corMethod pearson --whatToPlot scatterplot --skipZeros --removeOutliers -o $DIR/EC_pooled_scatter.pdf