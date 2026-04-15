#!/bin/bash
# s02_illumina_plink
#
# Cameron Brown 20Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s02_illumina_plink

# PBS directives
#---------------

#PBS -N s02_illumina_plink
#PBS -l nodes=1:ncpus=8
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

# Start message
echo "========================================"
echo "STEP: Illumina VCF to PLINK conversion started"
date
echo "========================================"
echo ""

# Scheduler / thread debug
echo "PBS_NODEFILE: ${PBS_NODEFILE}"
echo "PBS_NCPUS: ${PBS_NCPUS}"
echo "NCPUS: ${NCPUS}"
echo "Threads used: ${threads}"
echo ""

# Folders and files
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"
container="${base_folder}/containers/plink.sif"
input_vcf="${base_folder}/data/processed/illumina/clean/illumina.clean.vcf.gz"
output_folder="${base_folder}/data/processed/illumina/plink"
output_prefix="${output_folder}/illumina_clean"

# Make output folder if needed
mkdir -p "${output_folder}"

# Check inputs exist
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "Input VCF: ${input_vcf}"
echo "Output prefix: ${output_prefix}"
echo "Container: ${container}"
echo ""

# Convert cleaned VCF to PLINK bed/bim/fam
singularity exec --bind /mnt/beegfs "${container}" plink2 \
  --vcf "${input_vcf}" \
  --geno 0.1 \
  --mind 0.1 \
  --make-bed \
  --threads "${threads}" \
  --out "${output_prefix}"

echo ""
echo "PLINK conversion complete"
date
echo ""

# Validation block
echo "========================================"
echo "VALIDATION: Checking PLINK outputs"
date
echo "========================================"
echo ""

echo "Checking output files exist..."
if [ -f "${output_prefix}.bed" ] && [ -f "${output_prefix}.bim" ] && [ -f "${output_prefix}.fam" ]; then
  echo "PASS: .bed, .bim and .fam files found"
else
  echo "FAIL: One or more PLINK output files are missing"
  exit 1
fi
echo ""

echo "Counting variants from .bim file..."
variant_count=$(wc -l < "${output_prefix}.bim")
echo "Variants: ${variant_count}"
echo ""

echo "Counting samples from .fam file..."
sample_count=$(wc -l < "${output_prefix}.fam")
echo "Samples: ${sample_count}"
echo ""

echo "Preview of first 5 variants (.bim):"
head -5 "${output_prefix}.bim"
echo ""

echo "Preview of first 5 samples (.fam):"
head -5 "${output_prefix}.fam"
echo ""

echo "========================================"
echo "STEP COMPLETE: Illumina PLINK conversion successful"
date
echo "========================================"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"