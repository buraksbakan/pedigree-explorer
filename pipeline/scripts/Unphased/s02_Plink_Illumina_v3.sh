#!/bin/bash
# s02_illumina_plink_v3
#
# Cameron Brown 20Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s02_illumina_plink_v3

# PBS directives
#---------------

#PBS -N s02_illumina_plink_v3
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -m abe
#PBS -M cameron.brown.944@cranfield.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Calculate number of threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Stop at runtime errors
set -e

# Folders and files
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"

# Main Illumina preprocessing folder
pipeline_folder="${base_folder}/data/processed/Illumina_Preprocessing"

# Subfolders
clean_vcf_folder="${pipeline_folder}/clean_vcf"
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"

# Inputs 
container="${base_folder}/github_repo/pedigree-explorer/pipeline/config/Simple_Container.sif"
input_vcf="${clean_vcf_folder}/illumina.clean.vcf.gz"

# Outputs
output_prefix="${plink_files_folder}/illumina_clean"

# Create output folders
mkdir -p "${plink_files_folder}"

# Check inputs exist
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
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

# Convert cleaned VCF to PLINK bed/bim/fam
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