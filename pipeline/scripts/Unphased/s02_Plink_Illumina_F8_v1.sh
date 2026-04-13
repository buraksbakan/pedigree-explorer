#!/bin/bash
# s02_IHCAPX8_Illumina_Plink_v1
#
# Cameron Brown 2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s02_IHCAPX8_Illumina_Plink_v1

# PBS directives
#---------------

#PBS -N s02_ihcapx8_illumina_plink_v1
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

cd "$PBS_O_WORKDIR"

threads="${PBS_NCPUS:-${NCPUS:-1}}"

set -e

# Base
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"

# IHCAPX8 pipeline folder
pipeline_folder="${base_folder}/data/processed/IHCAPX8_Illumina_Preprocessing"

# Subfolders
clean_vcf_folder="${pipeline_folder}/clean_vcf"
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"

# Inputs
container="${base_folder}/github_repo/pedigree-explorer/pipeline/config/Simple_Container.sif"
input_vcf="${clean_vcf_folder}/IHCAPX8.clean.vcf.gz"

# Outputs
output_prefix="${plink_files_folder}/IHCAPX8_clean"

mkdir -p "${plink_files_folder}"

# Checks
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP: IHCAPX8 Illumina ? PLINK"
date
echo "Input: $(basename "${input_vcf}")"
echo "Output: $(basename "${output_prefix}")"
echo "Threads: ${threads}"
echo "----------------------------------------"
echo ""

# Convert
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

# Validation
echo "----------------------------------------"
echo "VALIDATION"
date
echo "----------------------------------------"

if [ -f "${output_prefix}.bed" ] && [ -f "${output_prefix}.bim" ] && [ -f "${output_prefix}.fam" ]; then
  echo "PASS: PLINK files created"
else
  echo "FAIL: Missing PLINK files"
  exit 1
fi

echo ""
echo "Variant count:"
wc -l "${output_prefix}.bim"

echo ""
echo "Sample count:"
wc -l "${output_prefix}.fam"

echo ""
echo "Preview (.bim):"
head -5 "${output_prefix}.bim"

echo ""
echo "Preview (.fam):"
head -5 "${output_prefix}.fam"

echo ""
echo "----------------------------------------"
echo "STEP COMPLETE"
date
echo "----------------------------------------"

rm -f "$PBS_O_WORKDIR/$PBS_JOBID"