#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=4:mem=24gb:ngpus=1:gpu_type=L40S
#PBS -N variant-scorer
#PBS -J 0-4

tid=$PBS_ARRAY_INDEX

# --- Edit to your WORKDIR ---
WORKDIR=/rds/general/user/tv722/projects/tamara_vujic_phd/live/ChromBPNet/chrombpnet0.5
NAME=chr21_dbSnp155Common_lymphoedema

# --- Set paths ---
genomeDir=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38_gencode44/
sifPATH=/rds/general/project/cebola_lab_general/live/variant-scorer

# Run the variant-scorer pipeline
# variant_scoring.py
apptainer exec --nv \
	--bind "${genomeDir}":/genome \
	--bind "${WORKDIR}":/work \
	--env PYTHONPATH=/scratch/variant-scorer/src \
	"${sifPATH}"/variant-scorer.sif \
	python /scratch/variant-scorer/src/variant_scoring.py \
	-l /work/"${NAME}".tsv \
	-g /genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	-s /work/hg38.chrom.sizes \
	-m /work/results_fold_"${tid}"/models/chrombpnet_nobias.h5 \
	-o /work/results_fold_"${tid}"/"${NAME}"_variant_scores \
	-sc chrombpnet \
	-dm