#!/bin/bash
# s02_Illumina_PLINK
#
# Cameron Brown 2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s02_Illumina_PLINK

# PBS directives
#---------------

#PBS -N s02_Illumina_PLINK
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -m abe
#PBS -M your.email@yourinstitution.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

# ----------------------------
# INSTRUCTIONS
# ----------------------------
# 1. Run s01_Filter_Illumina before running this script.
# 2. Edit ONLY the variables in the "USER INPUTS" section below.
# 3. Set:
#    - base_folder: your main project directory
#    - container: full path to the Singularity/Apptainer container with PLINK2
# 4. This script expects the filtered VCF produced by s01_Filter_Illumina:
#       ${base_folder}/data/processed/Illumina_Preprocessing/clean_vcf/illumina_filtered.vcf.gz
# 5. The container definition (.def file) used to build the container is available on the project GitHub.
# 6. Do NOT modify anything below the "DO NOT EDIT" line unless you understand the pipeline.
# 7. Submit the script using:
#       qsub s02_Illumina_PLINK
#
# Output:
# - PLINK .bed file
# - PLINK .bim file
# - PLINK .fam file
# - Basic validation printed to stdout
# ----------------------------

# ----------------------------
# USER INPUTS (EDIT THESE ONLY)
# ----------------------------

base_folder="/path/to/project"
container="/path/to/container.sif"

# ----------------------------
# DO NOT EDIT BELOW
# ----------------------------

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Calculate number of threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Stop at runtime errors
set -e

# Main Illumina preprocessing folder
pipeline_folder="${base_folder}/data/processed/Illumina_Preprocessing"

# Subfolders
clean_vcf_folder="${pipeline_folder}/clean_vcf"
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"

# Inputs
input_vcf="${clean_vcf_folder}/illumina_filtered.vcf.gz"

# Outputs
output_prefix="${plink_files_folder}/illumina_filtered"

# Create output folders
mkdir -p "${plink_files_folder}"

# Check inputs exist
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  echo "Run s01_Filter_Illumina first, or check base_folder is correct."
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP: Illumina VCF to PLINK conversion"
date
echo "Input VCF: $(basename "${input_vcf}")"
echo "Output prefix: $(basename "${output_prefix}")"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

# Convert filtered VCF to PLINK bed/bim/fam
singularity exec --bind /mnt/beegfs "${container}" plink2 \
  --vcf "${input_vcf}" \
  --geno 0.1 \
  --mind 0.1 \
  --make-bed \
  --threads "${threads}" \
  --out "${output_prefix}"

echo "PLINK conversion complete"
date
echo ""

echo "----------------------------------------"
echo "VALIDATION"
date
echo "----------------------------------------"

echo "Checking output files..."
if [ -f "${output_prefix}.bed" ] && [ -f "${output_prefix}.bim" ] && [ -f "${output_prefix}.fam" ]; then
  echo "PASS: .bed, .bim and .fam files found"
else
  echo "FAIL: One or more PLINK output files are missing"
  exit 1
fi
echo ""

echo "Variant count (.bim):"
variant_count=$(wc -l < "${output_prefix}.bim")
echo "${variant_count}"
echo ""

echo "Sample count (.fam):"
sample_count=$(wc -l < "${output_prefix}.fam")
echo "${sample_count}"
echo ""

echo "First 5 variants (.bim):"
head -5 "${output_prefix}.bim"
echo ""

echo "First 5 samples (.fam):"
head -5 "${output_prefix}.fam"
echo ""

echo "----------------------------------------"
echo "STEP COMPLETE: Illumina PLINK conversion successful"
date
echo "----------------------------------------"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"