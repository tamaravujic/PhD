#PBS -l select=1:mem=64gb:ncpus=16
#PBS -l walltime=03:00:00
#PBS -N corr_replicates_EC

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/tv722/home/anaconda3/envs/ATAC

DIR=/rds/general/user/tv722/projects/erg_lymphangiogenesis_project/live/ECs_ATAC/imlec_atac_comparison
REFDIR=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38

cd $DIR

multiBigwigSummary bins -b $DIR/bigwigs/LEC_rep1.BPM.bw $DIR/bigwigs/LEC_rep2.BPM.bw $DIR/bigwigs/imLEC_rep1.BPM.bw $DIR/bigwigs/imLEC_rep2.BPM.bw --labels HDLEC1 HDLEC2 imLEC1 imLEC2 --binSize 10000 --blackListFileName $REFDIR/hg38-blacklist.v2.bed -p 16 -out $DIR/EC_reps_corr.npz --outRawCounts $DIR/EC_reps_bins.tsv

plotCorrelation -in $DIR/EC_reps_corr.npz --corMethod pearson --whatToPlot heatmap --plotNumbers --skipZeros --removeOutliers --plotTitle "Pearson correlation between replicates, 10kb bins" -o $DIR/EC_reps_heatmap.pdf --outFileCorMatrix $DIR/EC_reps_Pearson.tsv

plotCorrelation -in $DIR/EC_reps_corr.npz --corMethod pearson --whatToPlot scatterplot --skipZeros --removeOutliers -o $DIR/EC_reps_scatter.pdf