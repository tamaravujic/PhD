#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=4:mem=50gb:ngpus=1:gpu_type=L40S
#PBS -N train_bias_model
#PBS -J 0-4

# Training models for five folds (0-4)
tid=$PBS_ARRAY_INDEX

# ============================
# Script: train_bias_model_sif.sh
# Purpose: Train a bias model for use with ChromBPNet, using the Imperial HPC with Singularity and GPU support.
# Usage:
#   qsub -v BAM=input.bam,PEAKS=peaks.bed,WORKDIR=/path/to/work,THRESHOLD=0.1,NAME=myproject train_bias_model_sif.sh
# Environment variables:
#   BAM      : path to input BAM file (required).
#   PEAKS    : path to input peaks file in narrowPeak format (required)
#   WORKDIR  : working directory (optional, default used if unset)
#   THRESHOLD: optional threshold value (default=0.1)
#   NAME     : name prefix for output files (optional, default="myproject")
# ============================

# Function to print help
show_help() {
    echo "Usage: qsub -v BAM=input.bam,PEAKS=peaks.bed,WORKDIR=/path/to/work,THRESHOLD=0.1,NAME=myproject train_bias_model_sif.sh"
    echo
    echo "Environment variables:"
    echo "  BAM      : path to input BAM file (required). Should be filtered and indexed."
    echo "  PEAKS    : path to input peaks file (required). Expecting .narrowPeak format, recommended MACS2 p = 0.01, with blacklist regions removed."
    echo "  WORKDIR  : working directory (optional). Default: directory created in the cebola_lab_general ephemeral space)"
    echo "  THRESHOLD: optional threshold value (default=0.1)"
    echo "  NAME     : name prefix for output files (optional, default=\"myproject\")"
    echo
    echo "Job uses the following set paths:
    genomeDir=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38_gencode44/
    sifPATH=/rds/general/project/cebola_lab_general/live/chrombpnet"
    echo
    exit 0
}

# Set the CUDA environment variables
OMPI_MCA_opal_cuda_support=true
export OMP_NUM_THREADS=4
OMPI_MCA_pml="ucx" 
OMPI_MCA_osc="ucx"
UCX_MEMTYPE_CACHE=n

# Check if user requested help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
fi

# Set paths to reference files
genomeDir=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38_gencode44/
sifPATH=/rds/general/project/cebola_lab_general/live/chrombpnet

# If WORKDIR has been set as an argument, use that. Otherwise, use the default ephemeral directory.
if [ -n "$WORKDIR" ]; then
	echo "Using WORKDIR from argument: $WORKDIR"
else
	cd /rds/general/project/cebola_lab_general/ephemeral
	# Make chromBPNet directory if it doesn't exist
	if [ ! -d chromBPNet ]; then
		mkdir -p chromBPNet
	fi
	cd chromBPNet
	# Make a unique directory for this job
	mkdir bias_model_"${PBS_JOBID}"
	cd bias_model_"${PBS_JOBID}"
	WORKDIR=$(pwd)
	echo "Using default WORKDIR: ${WORKDIR}"
fi

# Set project name
if [ ! -n "$NAME" ]; then
	NAME="myproject"
	echo "NAME not set. Using default: ${NAME}"
else
	echo "Using NAME from argument: ${NAME}"
fi

# ---------------------------------------

echo "Job started on $(date)"
echo "Output files are located at ${WORKDIR}"

# Get the ncpus from qstat -f $PBS_JOBID
NCPUS=$(qstat -f $PBS_JOBID | grep "Resource_List.ncpus" | cut -d "=" -f 2)
# Get the memory from qstat -f $PBS_JOBID
MEMORY=$(qstat -f $PBS_JOBID | grep "Resource_List.mem" | cut -d "=" -f 2)
# -------------

# Exit if no GPUs are available
# nvidia-smi
if [ $(apptainer exec --nv "${sifPATH}"/chrombpnet.sif python -c "import tensorflow as tf; print(len(tf.config.list_physical_devices('GPU')))") -eq 0 ]; then
    echo "No GPUs available. Exiting."
    exit 1
fi

#################################################

###### Create the chromosome splits for the training, validation and testing sets
# Define the genome directory that contains the reference genome (FASTA) and the chrom.sizes file
# The chrom.sizes file can be generated using the following command: awk '{print $1"\t"$2}' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > hg38.chrom.sizes
# Filter the chrom.sizes file to remove unwanted chromosomes
cat "${genomeDir}"/hg38.chrom.sizes | grep -ve chrUn -ve chrEBV -ve random -ve chrM > hg38.chrom.sizes

# --- Folds ---
# Copy the five folds json files if they do not already exist
cp --force -r /rds/general/project/cebola_lab_general/live/chrombpnet/splits .

# Define train, validation and test chromosome splits 
# Note: To create the folds we want to specify only the training, validation and testing chromosomes. For the human chromosomes we use chr1 to chr22, chrX and chrY and filter the remaining.

###### Train a custom bias model
# Copy the peaks file to the current directory
if [ ! -f peaks.narrowPeak ]; then
    echo "Copying peaks file to working directory..."
    cp "${PEAKS}" peaks.narrowPeak
else
    echo "Peaks file already exists in working directory. Skipping copy."
fi

# Generate non-peak regions which are GC matched
# -fl CHR_FOLD_PATH, --chr-fold-path CHR_FOLD_PATH, Fold information - dictionary with test, valid and train keys and values with corresponding chromosomes
# -p: Input peaks in narrowPeak file format with 10 columns

# Run chrombpnet prep nonpeaks
# Use the --nv flag for compatibility with an NVIDIA GPU
if [ ! -f output_fold_"${tid}"_negatives.bed ]; then
    echo "Generating non-peak regions..."
    apptainer exec --nv \
	--bind "${genomeDir}":/genome \
	--bind $(pwd):/work \
	"${sifPATH}"/chrombpnet.sif \
	chrombpnet prep nonpeaks \
	-g /genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	-p /work/peaks.narrowPeak \
	-c /work/hg38.chrom.sizes \
	-fl /work/splits/fold_"${tid}".json \
	-br /genome/hg38-blacklist.v2.bed \
	-o /work/output_fold_"${tid}"
else
    echo "Non-peak regions file already exists. Skipping this step."
fi

# Run chrombpnet bias pipeline
if [ ! -f input.bam ]; then
    echo "Copying input BAM file to working directory..."
    cp "${BAM}" input.bam
fi

# Check for BAM index file
if [ -f "${BAM}".bai ]; then
	echo "Copying input BAM index file to working directory..."
	cp "${BAM}".bai input.bam.bai
else
	echo "Input BAM index file not found. Indexing now..."
	samtools index -@ "${NCPUS}" input.bam
fi

# Get current fold json file
foldfile="splits/fold_${tid}.json" #$(ls splits/fold_*.json | sed -n "${tid}p")
echo "Using fold file: ${foldfile}"

apptainer exec --nv \
	--bind "${genomeDir}":/genome \
	--bind $(pwd):/work \
	"${sifPATH}"/chrombpnet.sif \
	chrombpnet bias pipeline \
  	-ibam /work/input.bam \
  	-d "ATAC" \
  	-g /genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  	-c /work/hg38.chrom.sizes \
  	-p /work/peaks.narrowPeak \
  	-n /work/output_fold_"${tid}"_negatives.bed \
  	-fl /work/"${foldfile}" \
  	-b "${THRESHOLD}" \
  	-o /work/output_fold_"${tid}" \
  	-fp "${NAME}"_"${THRESHOLD}"
