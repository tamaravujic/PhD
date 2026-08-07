#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=12:mem=50gb:ngpus=1:gpu_type=L40S
#PBS -N chrombpnet_pipeline
#PBS -e chrombpnet_pipeline.error
#PBS -o chrombpnet_pipeline.output

# ============================
# Script: chrombpnet_pipeline_sif.sh
# Purpose: Run ChromBPNet prep workflow on HPC
# Usage:
#   qsub -v WORKDIR=/path/to/work,MODEL=XXXX_XX_bias.h5,FOLD=fold_0 chrombpnet_pipeline_sif.sh
# Environment variables:
#   WORKDIR  : working directory with the required input files. If the train_bias_model_sif.sh script was run, WORKDIR should be set as the directory generated from the job. Must contain peaks.narrowPeak, input.bam and input.bam.bai. (required).
#   MODEL    : Name of bias model file. Should be a file name within output_fold_x/models, e.g. XXXX_XX_bias.h5 (required).
#   FOLD     : Specify training fold (fold_0, fold_1, fold_2, fold_3 or fold_4). Results folder will be called results_fold_x. Default is fold_0.
# ============================

# Function to print help
show_help() {
    echo "Usage: qsub -v WORKDIR=/path/to/work,MODEL=XXXX_XX_bias.h5,FOLD=fold_0 chrombpnet_pipeline_sif.sh"
    echo
    echo "Environment variables:
    WORKDIR  : Working directory. Use folder created by train_bias_model_sif.sh if run. Must contain peaks.narrowPeak, input.bam and input.bam.bai. (required).
    MODEL    : Name of bias model file. Should be a file name within output_fold_x/models/, e.g. XXXX_XX_bias.h5 (required).
    FOLD     : Specify training fold (fold_x). Results folder will be called results_x. Default is fold_0."
    echo
    echo "Job uses the following set paths:
    genomeDir=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38_gencode44/
    sifPATH=/rds/general/project/cebola_lab_general/live/chrombpnet"
    echo
    echo "For best practise, use a WORKDIR in the ephemeral project space and copy results back to permanent live storage after job completion."
    echo "ESSENTIAL: the results folder must NOT already exist in WORKDIR otherwise the job will exit. For failed attempts, delete the existing results folder before re-running."
    echo
    exit 0
}

# ----------- Job parameter setup -----------
# -------------------------------------------
# Get the ncpus from qstat -f $PBS_JOBID
NCPUS=$(qstat -f $PBS_JOBID | grep "Resource_List.ncpus" | cut -d "=" -f 2)
# Get the memory from qstat -f $PBS_JOBID
MEMORY=$(qstat -f $PBS_JOBID | grep "Resource_List.mem" | cut -d "=" -f 2)

# Set the CUDA environment variables
OMPI_MCA_opal_cuda_support=true
OMPI_MCA_pml="ucx" 
OMPI_MCA_osc="ucx"
UCX_MEMTYPE_CACHE=n
export OMP_NUM_THREADS="${NCPUS}"

# ----------- Check job parameters -----------
# --------------------------------------------

# Check if user requested help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
fi

# Set paths to reference files
genomeDir=/rds/general/project/cebola_lab_general/live/reference-genomes/GRCh38_gencode44/
sifPATH=/rds/general/project/cebola_lab_general/live/chrombpnet

# Set the working directory
if [ ! -n "$WORKDIR" ] || [ ! -d "$WORKDIR" ]; then
	echo "WORKDIR not set or does not exist. Exiting."
	exit 1
fi

# Change to the working directory
echo "Changing to WORKDIR... $WORKDIR"
cd "${WORKDIR}"

# Check if the required input files are present
if [ ! -f peaks.narrowPeak ] || [ ! -f input.bam ] || [ ! -f input.bam.bai ]; then
    echo "Required input files not found in WORKDIR. Confirm that peaks.narrowPeak, input.bam and input.bam.bai are present and named correctly. Exiting."
    exit 1
fi

# Confirm that MODEL variable is set
if [ ! -n "${MODEL}" ] || [ ! -f "output_${FOLD}/models/${MODEL}" ]; then
	echo "MODEL file name not set or file does not exist. Name provided: ${MODEL}. Check that output_${FOLD}/models/${MODEL} exists. Exiting."
	exit 1
fi
echo "Using MODEL: output_${FOLD}/models/${MODEL}"

# Check that FOLD variable is set
if [ ! -n "${FOLD}" ] ; then
	echo "WARNING. FOLD not specified, using default: fold_0. Exiting."
	FOLD="fold_0"
fi

# Exit if no GPUs are available
if [ $(apptainer exec --nv "${sifPATH}"/chrombpnet.sif python -c "import tensorflow as tf; print(len(tf.config.list_physical_devices('GPU')))") -eq 0 ]; then
    echo "No GPUs available. Exiting."
    exit 1
fi
# -------------------------------------------

if [ ! -f hg38.chrom.sizes ]; then
	echo "hg38.chrom.sizes not found in WORKDIR. Creating hg38.chrom.sizes file using cat "${genomeDir}"/hg38.chrom.sizes | grep -ve chrUn -ve chrEBV -ve random -ve chrM..."
	# Create hg38.chrom.sizes file
	cat "${genomeDir}"/hg38.chrom.sizes | grep -ve chrUn -ve chrEBV -ve random -ve chrM > hg38.chrom.sizes
fi

# If a pre-trained bias model is being used, then prepare the splits and non-peak regions if not already present
#if [ ! -f output_negatives.bed ]; then
#	echo "output_negatives.bed not found in ${WORKDIR}. Running prep nonpeaks steps..."

	# If the chrombpnet prep splits and chrombpnet prep nonpeaks jobs were not run, run them now.	
#	echo "Copying splits for training..."
#	cp --force -r /rds/general/project/cebola_lab_general/live/chrombpnet/splits .

#	echo "Generating non-peak regions..."
#	apptainer exec --nv \
#	--bind "${genomeDir}":/genome \
#	--bind $(pwd):/work \
#	"${sifPATH}"/chrombpnet.sif \
#	chrombpnet prep nonpeaks \
#	-g /genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#	-p /work/peaks.narrowPeak \
#	-c /work/hg38.chrom.sizes \
#	-fl /work/splits/fold_0.json \
#	-br /genome/hg38-blacklist.v2.bed \
#	-o /work/output
#fi

# If result folder already exists, exit to avoid overwriting
RESULTS="results_${FOLD}"

if [ -d "${RESULTS}" ]; then
    echo "Results directory already exists. Skipping ChromBPNet pipeline run to avoid overwriting existing results."
    exit 0
fi

echo "Running ChromBPNet pipeline..."
# Run the chrombpnet pipeline
apptainer exec --nv \
	--bind "${genomeDir}":/genome \
	--bind $(pwd):/work \
	"${sifPATH}"/chrombpnet.sif \
	chrombpnet pipeline \
	-ibam /work/input.bam \
	-d "ATAC" \
	-g /genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
	-c /work/hg38.chrom.sizes \
	-p /work/peaks.narrowPeak \
	-n /work/output_"${FOLD}"_negatives.bed \
	-fl /work/splits/"${FOLD}".json \
	-b /work/output_"${FOLD}"/models/"${MODEL}" \
	-o /work/"${RESULTS}"

# --- END ---